;;;;; Pivot storage



(defstruct (pivot-key
	     (:constructor %make-pivot-key))
  (j  (error "pivot-key constructor") :type fixnum)
  (i  (error "pivot-key constructor") :type fixnum))

  

(defun make-pivot-key (i j)
  (%make-pivot-key
   :j  j
   :i  i))



(defstruct (pivot-value
	     (:constructor %make-pivot-value))
  (mc (error "pivot-value constructor") :type fixnum)
  (v  (error "pivot-value constructor") :type float)
  (ci (error "pivot-value constructor") :type fixnum))



(defun make-pivot-value (mc v ci)
  (%make-pivot-value
   :mc mc
   :v v
   :ci ci))



(defstruct (pivot-store
	     (:constructor %make-pivot-store))
  (buckets    (error "pivot-store construction error") 
	          :type hash-table)
  (sorted-mcs #() :type vector))



(defun make-pivot-store ()
  (%make-pivot-store
   :buckets    (make-hash-table)
   :sorted-mcs (make-vector fixnum)))



(defun make-bucket ()
  (make-hash-table :test #'equalp))



(defun pivot-store-sort-mcs (ps)
  (setf (pivot-store-sorted-mcs ps) 
	(sort (pivot-store-sorted-mcs ps) #'<)))



(defun reset-pivot-store (ps rcs ccs)
  (maphash #'(lambda (key bucket) (clrhash bucket))
	   (pivot-store-buckets ps))
  (let ((m (length rcs)))
    (dotimes (i m)
      (let ((rc (aref rcs i)))
	(when (zerop rc)
	  (return-from reset-pivot-store))
	(dotimes (j m)
	  (let* ((cc (aref ccs j))
		 (mc (* (- rc 1) (- cc 1))))
	    (when (zerop cc)
	      (return-from reset-pivot-store))
	    (multiple-value-bind (bucket there)
		(gethash mc (pivot-store-buckets ps))
	      (unless there
		(setf (gethash mc (pivot-store-buckets ps)) 
		      (make-bucket)))))))))
  (setf (fill-pointer (pivot-store-sorted-mcs ps)) 0)
  (maphash #'(lambda (mc bucket)
	       (vector-push-extend mc (pivot-store-sorted-mcs ps)))
	   (pivot-store-buckets ps))
  (pivot-store-sort-mcs ps))
	

(defun pivot-store-remove (ps piv-key mc)
  (multiple-value-bind (old-bucket there)
      (gethash mc (pivot-store-buckets ps))
    (assert there)
    (let ((rval (remhash piv-key old-bucket)))
      (assert rval))
    (when (zerop (hash-table-count old-bucket))
      (let ((rval (remove-in-increasing-vector mc (pivot-store-sorted-mcs ps))))
	(assert rval)))))

      

(defun pivot-store-add (ps piv-key piv-val)
  (let* ((mc (pivot-value-mc piv-val))
	 (key-index (find-index (pivot-store-sorted-mcs ps) mc)))
    (when (= -1 key-index)
      ;; if nonexisting bucket, create it
      (setf (gethash mc (pivot-store-buckets ps)) (make-bucket))
      (insert-in-increasing-vector mc (pivot-store-sorted-mcs ps)))
    ;; add to bucket
    (setf piv-val (gethash piv-key (gethash mc (pivot-store-buckets ps))))))

      

(defun pivot-store-change-v (ps piv-key mc new-v)
  (multiple-value-bind (bucket there)
      (gethash mc (pivot-store-buckets ps))
    (assert there)
    (multiple-value-bind (piv-val there)
	(gethash piv-key bucket)
      (assert there)
      (setf (pivot-value-v piv-val) new-v)
      (setf piv-val (gethash piv-key bucket)))))



(defun pivot-store-change-mc (ps piv-key mc new-mc)
  (multiple-value-bind (old-bucket there)
      (gethash mc (pivot-store-buckets ps))
    (assert there)
    (multiple-value-bind (piv-val there)
	(gethash piv-key old-bucket)
      (assert there)
      (remhash piv-key old-bucket)
      (when (zerop (hash-table-count old-bucket))
	(let ((rval (remove-in-increasing-vector mc (pivot-store-sorted-mcs ps))))
	  (assert rval)))
      (setf (pivot-value-mc piv-val) new-mc)
      (pivot-store-add ps piv-key piv-val))))
	



