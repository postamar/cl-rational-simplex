;;;;; Dynamic Markowitz pivot implementation
;;;;;



;;;; Adds a pivot in the appropriate bucket
(defun pivot-add (b i j ci)
  (let* ((m (basis-matrix-size b))
	 (key (+ (* m i) j))
	 (rc (aref (basis-matrix-row-nnz b) i))
	 (cc (aref (basis-matrix-col-nnz b) j))
	 (mc (* (- rc 1) (- cc 1))))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    (multiple-value-bind (bucket-index bucket-there)
	(markowitz-tree-splay (basis-matrix-pivot-buckets b) mc)
      (unless bucket-there
	(setf bucket-index
	      (markowitz-tree-set (basis-matrix-pivot-buckets b) mc (make-pivot-bucket))))
      (let ((bucket (markowitz-tree-value (basis-matrix-pivot-buckets b) bucket-index)))
	(multiple-value-bind (pivot-index pivot-notthere)
	    (pivot-bucket-set bucket key ci)
	  (declare (ignore pivot-index))
	  (assert pivot-notthere))))))



;;;; Removes a pivot from its bucket
(defun pivot-remove (b i j)
  (let* ((m (basis-matrix-size b))
	 (key (+ (* m i) j))
	 (rc (aref (basis-matrix-row-nnz b) i))
	 (cc (aref (basis-matrix-col-nnz b) j))
	 (mc (* (- rc 1) (- cc 1))))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    (multiple-value-bind (bucket-index bucket-there)
	(markowitz-tree-splay (basis-matrix-pivot-buckets b) mc)
      (unless bucket-there
	(setf bucket-index
	      (markowitz-tree-set (basis-matrix-pivot-buckets b) mc (make-pivot-bucket))))
      (let ((bucket (markowitz-tree-value (basis-matrix-pivot-buckets b) bucket-index)))
	(multiple-value-bind (pivot-index pivot-there)
	    (pivot-bucket-remove bucket key)
	  (if pivot-there
	      (pivot-bucket-value bucket pivot-index)
	      -1))))))



;;;; Adds a new pivot
(defun pivot-add-new (b i j ci)
  (let ((row-js (aref (basis-matrix-row-js b) i))
	(col-is (aref (basis-matrix-col-is b) j))
	(rc (aref (basis-matrix-row-nnz b) i))
	(cc (aref (basis-matrix-col-nnz b) j)))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    ;;;;
    (setf (fill-pointer (basis-matrix-pivot-i-queue b)) 0
	  (fill-pointer (basis-matrix-pivot-ci-queue b)) 0
	  (fill-pointer (basis-matrix-pivot-j-queue b)) 0)
    (vector-push-extend i (basis-matrix-pivot-i-queue b))
    (vector-push-extend ci (basis-matrix-pivot-ci-queue b))
    (vector-push-extend j (basis-matrix-pivot-j-queue b))
    (dotimes (index (length col-is))
      (let* ((oi (aref col-is index))
	     (orc (aref (basis-matrix-row-nnz b) oi))
	     (ci (pivot-remove b oi j)))
	(assert (< 0 orc most-positive-fixnum))
	(when (/= -1 ci)
	  (vector-push-extend oi (basis-matrix-pivot-i-queue b))
	  (vector-push-extend ci (basis-matrix-pivot-ci-queue b))
	  (vector-push-extend j (basis-matrix-pivot-j-queue b)))))
    (dotimes (index (length row-js))
      (let* ((oj (aref row-js index))
	     (occ (aref (basis-matrix-col-nnz b) oj))
	     (ci (pivot-remove b i oj)))
	(assert (< 0 occ most-positive-fixnum))
	(when (/= -1 ci) 
	  (vector-push-extend i (basis-matrix-pivot-i-queue b))
	  (vector-push-extend ci (basis-matrix-pivot-ci-queue b))
	  (vector-push-extend oj (basis-matrix-pivot-j-queue b)))))
    (incf (aref (basis-matrix-row-nnz b) i))
    (incf (aref (basis-matrix-col-nnz b) j))
    (vector-push-extend i (aref (basis-matrix-col-is b) j))
    (vector-push-extend j (aref (basis-matrix-row-js b) i))
    (vector-push-extend ci (aref (basis-matrix-row-cis b) i))
    ;;;;
    (dotimes (index (length (basis-matrix-pivot-ci-queue b)))
      (pivot-add b
		 (aref (basis-matrix-pivot-i-queue b) index)
		 (aref (basis-matrix-pivot-j-queue b) index)
		 (aref (basis-matrix-pivot-ci-queue b) index)))))


;;;; Fill basis according to basis header
(defun fill-basis-matrix (b lp basis-header)
  (assert (= (basis-matrix-size b) (length basis-header)))
  (reset-basis-matrix b)
  (dotimes (j (basis-matrix-size b))
    (let ((col-ref (aref basis-header j))
	  (l_j     (aref (basis-matrix-l-columns b) j)))
      (let* ((col          (aref (lp-columns lp) col-ref))
	     (col-row-refs (column-row-refs col))
	     (col-vals     (column-values col))
	     (n-nz         (length col-row-refs)))
	(setf (lu-eta-matrix-coef l_j)    (column-coef col)
	      (lu-eta-matrix-col-ref l_j) col-ref)
	(dotimes (k n-nz)
	  (let ((i (aref (lp-active-row-inds lp) (aref col-row-refs k)))
		(ci (length (lu-eta-matrix-is l_j)))
		(val     (aref col-vals k)))
	    (unless (= -1 i)
	      (vector-push-extend i   (lu-eta-matrix-is  l_j))
	      (vector-push-extend val (lu-eta-matrix-vis l_j))
	      (vector-push-extend 0.0 (lu-eta-matrix-vfs l_j))
	      (incf (aref (basis-matrix-col-nnz b) j))
	      (incf (aref (basis-matrix-row-nnz b) i))
	      (vector-push-extend i   (aref (basis-matrix-col-is b) j))
	      (vector-push-extend j   (aref (basis-matrix-row-js b) i))
	      (vector-push-extend ci  (aref (basis-matrix-row-cis b) i))))))))
  (dotimes (k (basis-matrix-size b))
    (cond ((zerop (aref (basis-matrix-col-nnz b) k))
	   (basis-matrix-column-is-redundant b k)
	   (return))
	  ((zerop (aref (basis-matrix-row-nnz b) k))
	   (basis-matrix-row-is-redundant b k)
	   (return))))
  (unless (basis-matrix-is-singular b)
    (dotimes (i (basis-matrix-size b))
      (let ((row-js (aref (basis-matrix-row-js b) i))
	    (row-cis (aref (basis-matrix-row-cis b) i)))
	(assert (= (length row-js) (length row-cis)))
	(dotimes (k (length row-js))
	  (pivot-add b i (aref row-js k) (aref row-cis k)))))
    t))



;;;; Update the pivot counts and the pivot heap
(defun pivot-count-update (b pivot-i pivot-j)
  ;; empty queue
  (setf (fill-pointer (basis-matrix-pivot-i-queue b)) 0
	(fill-pointer (basis-matrix-pivot-ci-queue b)) 0
	(fill-pointer (basis-matrix-pivot-j-queue b)) 0)
  ;; remove pivot 
  (pivot-remove b pivot-i pivot-j)
  ;; remove pivot candidates in pivot column and pivot row
  ;; and put pivots with modified markowitz count in queue
  (dotimes (index-col-is (length (aref (basis-matrix-col-is b) pivot-j)))
    (let ((i (aref (aref (basis-matrix-col-is b) pivot-j) index-col-is)))
      (dotimes (index-row-js (length (aref (basis-matrix-row-js b) i)))
	(let ((j (aref (aref (basis-matrix-row-js b) i) index-row-js)))
	  (cond ((= j pivot-j))
		((= i pivot-i)
		 (pivot-remove b i j))
		(t
		 (let ((ci (pivot-remove b i j)))
		   (unless (= -1 ci)
		     (vector-push-extend i  (basis-matrix-pivot-i-queue b))
		     (vector-push-extend ci (basis-matrix-pivot-ci-queue b))
		     (vector-push-extend j  (basis-matrix-pivot-j-queue b))))))))))
  (dotimes (index-row-js (length (aref (basis-matrix-row-js b) pivot-i)))
    (let ((j (aref (aref (basis-matrix-row-js b) pivot-i) index-row-js)))
      (dotimes (index-col-is (length (aref (basis-matrix-col-is b) j)))
	(let ((i (aref (aref (basis-matrix-col-is b) j) index-col-is)))
	  (cond ((= i pivot-i))
		((= j pivot-j)
		 (pivot-remove b i j))
		(t
		 (let ((ci (pivot-remove b i j)))
		    (unless (= -1 ci)
		      (vector-push-extend i (basis-matrix-pivot-i-queue b))
		      (vector-push-extend ci (basis-matrix-pivot-ci-queue b))
		      (vector-push-extend j (basis-matrix-pivot-j-queue b))))))))))
  ;; update sparse residual matrix and non-zero counts
  (setf (aref (basis-matrix-row-nnz b) pivot-i) most-positive-fixnum
	(aref (basis-matrix-col-nnz b) pivot-j) most-positive-fixnum)
  (dotimes (index-col-is (length (aref (basis-matrix-col-is b) pivot-j)))
    (let ((i (aref (aref (basis-matrix-col-is b) pivot-j) index-col-is))
	  (pivot-index -1))
      (unless (= i pivot-i)
	(when (zerop (decf (aref (basis-matrix-row-nnz b) i)))
	  (basis-matrix-row-is-redundant b i)
	  (return-from pivot-count-update))
	(dotimes (index-row-js (length (aref (basis-matrix-row-js b) i)))
	  (let ((j (aref (aref (basis-matrix-row-js b) i) index-row-js)))
	    (when (= j pivot-j)
	      (setf pivot-index index-row-js)
	      (return))))
	(let ((last-j (vector-pop (aref (basis-matrix-row-js b) i)))
	      (last-ci (vector-pop (aref (basis-matrix-row-cis b) i))))
	  (unless (= last-j pivot-j)
	    (setf (aref (aref (basis-matrix-row-js b) i) pivot-index) last-j
		  (aref (aref (basis-matrix-row-cis b) i) pivot-index) last-ci))))))
  (dotimes (index-row-js (length (aref (basis-matrix-row-js b) pivot-i)))
    (let ((j (aref (aref (basis-matrix-row-js b) pivot-i) index-row-js))
	  (pivot-index -1))
      (unless (= j pivot-j)
	(when (zerop (decf (aref (basis-matrix-col-nnz b) j)))
	  (basis-matrix-column-is-redundant b j)
	  (return-from pivot-count-update))
	(dotimes (index-col-is (length (aref (basis-matrix-col-is b) j)))
	  (let ((i (aref (aref (basis-matrix-col-is b) j) index-col-is)))
	    (when (= i pivot-i)
	      (setf pivot-index index-col-is)
	      (return))))
	(let ((last-i (vector-pop (aref (basis-matrix-col-is b) j))))
	  (unless (= last-i pivot-i)
	    (setf (aref (aref (basis-matrix-col-is b) j) pivot-index) last-i))))))
  ;; move pivots in queue to buckets
  (dotimes (k (length (basis-matrix-pivot-i-queue b)))
    (let ((i (aref (basis-matrix-pivot-i-queue b) k))
	  (j (aref (basis-matrix-pivot-j-queue b) k))
	  (ci (aref (basis-matrix-pivot-ci-queue b) k)))
      (pivot-add b i j ci))))
	       
	       



;;;; Select the best pivot in the matrix
(defun basis-matrix-select-pivot (b)
  (let ((m (basis-matrix-size b))
	(is-zero nil)
	(zero-i -1)
	(zero-j -1)
	(zero-ci -1))
    (map-markowitz-tree 
     #'(lambda (mc bucket)
	(declare (ignore mc))
	 (map-pivot-bucket
	  #'(lambda (pivot-key ci)
	      (let* ((i (floor pivot-key m))
		     (j (- pivot-key (* i m)))
		     (l (aref (basis-matrix-l-columns b) j))
		     (vi (aref (lu-eta-matrix-vis l) ci)))
		(if (zerop vi)
		    (setf is-zero t
			  zero-i i
			  zero-j j
			  zero-ci ci)
		    (return-from basis-matrix-select-pivot
		      (values nil i j ci)))))
	  bucket))
     (basis-matrix-pivot-buckets b))
    (values is-zero zero-i zero-j zero-ci)))
    


;;;; Perform a pivot in the matrix
(defun basis-matrix-perform-pivot (b pk)
  (multiple-value-bind (is-zero i j ci)
      (basis-matrix-select-pivot b)
    (cond (is-zero
	   (basis-matrix-row-is-redundant b i)
	   (values -1 -1 -1))
	  (t
	   (pivot-count-update b i j)
	   (if (basis-matrix-is-singular b)
	       (values -1 -1 -1)
	       (let ((perm-i (aref (basis-matrix-i->pi b) i))
		     (perm-j (aref (basis-matrix-j->pj b) j))
		     (swap-i (aref (basis-matrix-pi->i b) pk))
		     (swap-j (aref (basis-matrix-pj->j b) pk)))
		 (setf (aref (basis-matrix-i->pi b) i) pk
		       (aref (basis-matrix-j->pj b) j) pk
		       (aref (basis-matrix-pi->i b) pk) i
		       (aref (basis-matrix-pj->j b) pk) j
		       (aref (basis-matrix-i->pi b) swap-i) perm-i
		       (aref (basis-matrix-j->pj b) swap-j) perm-j
		       (aref (basis-matrix-pi->i b) perm-i) swap-i
		       (aref (basis-matrix-pj->j b) perm-j) swap-j)
		 (values i j ci)))))))



(defun is-in-buckets (b i j)
  (map-markowitz-tree 
     #'(lambda (mc bucket)
	 (map-pivot-bucket
	  #'(lambda (pivot-key ci)
	      (declare (ignore ci))
	      (let* ((bi (floor pivot-key (basis-matrix-size b)))
		     (bj (- pivot-key (* i (basis-matrix-size b)))))
		(when (and (= i bi) (= j bj))
		  (return-from is-in-buckets mc))))
	  bucket))
     (basis-matrix-pivot-buckets b))
  nil)


(defun check-buckets (b)
  (map-markowitz-tree 
     #'(lambda (mc bucket)
	 (map-pivot-bucket
	  #'(lambda (pivot-key ci)
	    (declare (ignore ci))
	      (let* ((i (floor pivot-key (basis-matrix-size b)))
		     (j (- pivot-key (* i (basis-matrix-size b))))
		     (rc (aref (basis-matrix-row-nnz b) i))
		     (cc (aref (basis-matrix-col-nnz b) j)))
		(assert (< 0 rc most-positive-fixnum))
		(assert (< 0 cc most-positive-fixnum))
		(assert (= mc (* (- cc 1) (- rc 1))))))
	  bucket))
     (basis-matrix-pivot-buckets b)))
	       
  
(defun print-buckets (b)
  (format t "buckets: ---~%")
  (map-markowitz-tree 
   #'(lambda (mc bucket)
       (format t "bucket ~A:~%" mc)
       (map-pivot-bucket
	#'(lambda (pivot-key ci)
	    (declare (ignore ci))
	    (let* ((i (floor pivot-key (basis-matrix-size b)))
		   (j (- pivot-key (* i (basis-matrix-size b))))
		   (rc (aref (basis-matrix-row-nnz b) i))
		   (cc (aref (basis-matrix-col-nnz b) j)))
	      (assert (< 0 rc most-positive-fixnum))
	      (assert (< 0 cc most-positive-fixnum))
	      (assert (= mc (* (- rc 1) (- cc 1))))
	      (format t "(~A , ~A) ~A ?= ~A ~%" i j mc (* (- rc 1) (- cc 1)))))
	bucket))
   (basis-matrix-pivot-buckets b))
  (format t "~%"))
