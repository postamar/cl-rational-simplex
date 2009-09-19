;;;;; Dynamic Markowitz pivot implementation
;;;;;



;;;; Adds a pivot in the appropriate bucket
(defun pivot-add (bm i j ci)
  (let* ((m (basis-matrix-size bm))
	 (key (+ (* m i) j))
	 (rc (aref (basis-matrix-row-nnz bm) i))
	 (cc (aref (basis-matrix-col-nnz bm) j))
	 (mc (* (- rc 1) (- cc 1))))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    (multiple-value-bind (bucket-index bucket-there)
	(markowitz-tree-splay (basis-matrix-pivot-buckets bm) mc)
      (unless bucket-there
	(setf bucket-index
	      (markowitz-tree-set (basis-matrix-pivot-buckets bm) mc (make-pivot-bucket))))
      (let ((bucket (markowitz-tree-value (basis-matrix-pivot-buckets bm) bucket-index)))
	(multiple-value-bind (pivot-index pivot-notthere)
	    (pivot-bucket-set bucket key ci)
	  (declare (ignore pivot-index))
	  (assert pivot-notthere))))))



(declaim (ftype (function (basis-matrix fixnum fixnum) fixnum) pivot-remove))
;;;; Removes a pivot from its bucket
(defun pivot-remove (bm i j)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (declare (fixnum i j))
  (let* ((m (basis-matrix-size bm))
	 (key (+ (the fixnum (* m i)) j))
	 (rc (aref (basis-matrix-row-nnz bm) i))
	 (cc (aref (basis-matrix-col-nnz bm) j))
	 (mc (* (- rc 1) (- cc 1))))
    (declare (fixnum m key rc cc mc))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    (multiple-value-bind (bucket-index bucket-there)
	(markowitz-tree-splay (basis-matrix-pivot-buckets bm) mc)
      (unless bucket-there
	(setf bucket-index
	      (markowitz-tree-set (basis-matrix-pivot-buckets bm) mc (make-pivot-bucket))))
      (let ((bucket (markowitz-tree-value (basis-matrix-pivot-buckets bm) bucket-index)))
	(multiple-value-bind (pivot-index pivot-there)
	    (pivot-bucket-remove bucket key)
	  (if pivot-there
	      (pivot-bucket-value bucket pivot-index)
	      -1))))))



;;;; Adds a new pivot
(defun pivot-add-new (bm i j ci)
  (let ((row-js (aref (basis-matrix-row-js bm) i))
	(col-is (aref (basis-matrix-col-is bm) j))
	(rc (aref (basis-matrix-row-nnz bm) i))
	(cc (aref (basis-matrix-col-nnz bm) j)))
    (assert (and (< 0 rc) (< 0 cc)))
    (assert (and (< 0 rc most-positive-fixnum) (< 0 cc most-positive-fixnum)))
    ;;;;
    (setf (fill-pointer (basis-matrix-pivot-i-queue bm)) 0
	  (fill-pointer (basis-matrix-pivot-ci-queue bm)) 0
	  (fill-pointer (basis-matrix-pivot-j-queue bm)) 0)
    (vector-push-extend i (basis-matrix-pivot-i-queue bm))
    (vector-push-extend ci (basis-matrix-pivot-ci-queue bm))
    (vector-push-extend j (basis-matrix-pivot-j-queue bm))
    (dotimes (index (length col-is))
      (let* ((oi (aref col-is index))
	     (orc (aref (basis-matrix-row-nnz bm) oi))
	     (ci (pivot-remove bm oi j)))
	(assert (< 0 orc most-positive-fixnum))
	(when (/= -1 ci)
	  (vector-push-extend oi (basis-matrix-pivot-i-queue bm))
	  (vector-push-extend ci (basis-matrix-pivot-ci-queue bm))
	  (vector-push-extend j (basis-matrix-pivot-j-queue bm)))))
    (dotimes (index (length row-js))
      (let* ((oj (aref row-js index))
	     (occ (aref (basis-matrix-col-nnz bm) oj))
	     (ci (pivot-remove bm i oj)))
	(assert (< 0 occ most-positive-fixnum))
	(when (/= -1 ci) 
	  (vector-push-extend i (basis-matrix-pivot-i-queue bm))
	  (vector-push-extend ci (basis-matrix-pivot-ci-queue bm))
	  (vector-push-extend oj (basis-matrix-pivot-j-queue bm)))))
    (incf (aref (basis-matrix-row-nnz bm) i))
    (incf (aref (basis-matrix-col-nnz bm) j))
    (vector-push-extend i (aref (basis-matrix-col-is bm) j))
    (vector-push-extend j (aref (basis-matrix-row-js bm) i))
    (vector-push-extend ci (aref (basis-matrix-row-cis bm) i))
    ;;;;
    (dotimes (index (length (basis-matrix-pivot-ci-queue bm)))
      (pivot-add bm
		 (aref (basis-matrix-pivot-i-queue bm) index)
		 (aref (basis-matrix-pivot-j-queue bm) index)
		 (aref (basis-matrix-pivot-ci-queue bm) index)))))


;;;; Fill basis according to basis header
(defun fill-basis-matrix (bm lp basis-header)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (declare ((vector fixnum) basis-header))
  (assert (= (basis-matrix-size bm) (length basis-header)))
  (reset-basis-matrix bm)
  (dotimes (j (basis-matrix-size bm))
    (let ((col-ref (aref basis-header j))
	  (u_j     (aref (basis-matrix-u-columns bm) j)))
      (let* ((col          (aref (lp-columns lp) col-ref))
	     (col-row-refs (column-row-refs col))
	     (col-vals     (column-values col))
	     (n-nz         (length col-row-refs)))
	(setf (hsv-coef u_j)    (column-coef col))
	(dotimes (k n-nz)
	  (let ((i (aref (lp-active-row-inds lp) (aref col-row-refs k)))
		(ci (hsv-length u_j))
		(val     (aref col-vals k)))
	    (unless (= -1 i)
	      (vector-push-extend i   (hsv-is  u_j))
	      (vector-push-extend val (hsv-vis u_j))
	      (incf (aref (basis-matrix-col-nnz bm) j))
	      (incf (aref (basis-matrix-row-nnz bm) i))
	      (vector-push-extend i   (aref (basis-matrix-col-is bm) j))
	      (vector-push-extend j   (aref (basis-matrix-row-js bm) i))
	      (vector-push-extend ci  (aref (basis-matrix-row-cis bm) i))))))))
  (dotimes (k (basis-matrix-size bm))
    (cond ((zerop (aref (basis-matrix-col-nnz bm) k))
	   (basis-matrix-column-is-redundant bm k)
	   (return))
	  ((zerop (aref (basis-matrix-row-nnz bm) k))
	   (basis-matrix-row-is-redundant bm k)
	   (return))))
  (unless (basis-matrix-is-singular bm)
    (dotimes (i (basis-matrix-size bm))
      (let ((row-js (aref (basis-matrix-row-js bm) i))
	    (row-cis (aref (basis-matrix-row-cis bm) i)))
	(assert (= (length row-js) (length row-cis)))
	(dotimes (k (length row-js))
	  (pivot-add bm i (aref row-js k) (aref row-cis k)))))
    t))



;;;; Update the pivot counts and the pivot heap
(defun pivot-count-update (bm pivot-i pivot-j)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (declare (fixnum pivot-i pivot-j))
  ;; empty queue
  (setf (fill-pointer (basis-matrix-pivot-i-queue bm)) 0
	(fill-pointer (basis-matrix-pivot-ci-queue bm)) 0
	(fill-pointer (basis-matrix-pivot-j-queue bm)) 0)
  ;; remove pivot 
  (pivot-remove bm pivot-i pivot-j)
  ;; remove pivot candidates in pivot column and pivot row
  ;; and put pivots with modified markowitz count in queue
  (dotimes (index-col-is (length (aref (basis-matrix-col-is bm) pivot-j)))
    (let ((i (aref (aref (basis-matrix-col-is bm) pivot-j) index-col-is)))
      (dotimes (index-row-js (length (aref (basis-matrix-row-js bm) i)))
	(let ((j (aref (aref (basis-matrix-row-js bm) i) index-row-js)))
	  (declare (fixnum i j))
	  (cond ((= j pivot-j))
		((= i pivot-i)
		 (pivot-remove bm i j))
		(t
		 (let ((ci (pivot-remove bm i j)))
		   (unless (= -1 ci)
		     (vector-push-extend i  (basis-matrix-pivot-i-queue bm))
		     (vector-push-extend ci (basis-matrix-pivot-ci-queue bm))
		     (vector-push-extend j  (basis-matrix-pivot-j-queue bm))))))))))
  (dotimes (index-row-js (length (aref (basis-matrix-row-js bm) pivot-i)))
    (let ((j (aref (aref (basis-matrix-row-js bm) pivot-i) index-row-js)))
      (dotimes (index-col-is (length (aref (basis-matrix-col-is bm) j)))
	(let ((i (aref (aref (basis-matrix-col-is bm) j) index-col-is)))
	  (declare (fixnum i j))
	  (cond ((= i pivot-i))
		((= j pivot-j)
		 (pivot-remove bm i j))
		(t
		 (let ((ci (pivot-remove bm i j)))
		    (unless (= -1 ci)
		      (vector-push-extend i (basis-matrix-pivot-i-queue bm))
		      (vector-push-extend ci (basis-matrix-pivot-ci-queue bm))
		      (vector-push-extend j (basis-matrix-pivot-j-queue bm))))))))))
  ;; update sparse residual matrix and non-zero counts
  (setf (aref (basis-matrix-row-nnz bm) pivot-i) most-positive-fixnum
	(aref (basis-matrix-col-nnz bm) pivot-j) most-positive-fixnum)
  (dotimes (index-col-is (length (aref (basis-matrix-col-is bm) pivot-j)))
    (let ((i (aref (aref (basis-matrix-col-is bm) pivot-j) index-col-is))
	  (pivot-index -1))
      (declare (fixnum i))
      (unless (= i pivot-i)
	(when (zerop (decf (aref (basis-matrix-row-nnz bm) i)))
	  (basis-matrix-row-is-redundant bm i)
	  (return-from pivot-count-update))
	(dotimes (index-row-js (length (aref (basis-matrix-row-js bm) i)))
	  (let ((j (aref (aref (basis-matrix-row-js bm) i) index-row-js)))
	    (declare (fixnum j))
	    (when (= j pivot-j)
	      (setf pivot-index index-row-js)
	      (return))))
	(let ((last-j (vector-pop (aref (basis-matrix-row-js bm) i)))
	      (last-ci (vector-pop (aref (basis-matrix-row-cis bm) i))))
	  (declare (fixnum last-j))
	  (unless (= last-j pivot-j)
	    (setf (aref (aref (basis-matrix-row-js bm) i) pivot-index) last-j
		  (aref (aref (basis-matrix-row-cis bm) i) pivot-index) last-ci))))))
  (dotimes (index-row-js (length (aref (basis-matrix-row-js bm) pivot-i)))
    (let ((j (aref (aref (basis-matrix-row-js bm) pivot-i) index-row-js))
	  (pivot-index -1))
      (declare (fixnum j))
      (unless (= j pivot-j)
	(when (zerop (decf (aref (basis-matrix-col-nnz bm) j)))
	  (basis-matrix-column-is-redundant bm j)
	  (return-from pivot-count-update))
	(dotimes (index-col-is (length (aref (basis-matrix-col-is bm) j)))
	  (let ((i (aref (aref (basis-matrix-col-is bm) j) index-col-is)))
	    (declare (fixnum i))
	    (when (= i pivot-i)
	      (setf pivot-index index-col-is)
	      (return))))
	(let ((last-i (vector-pop (aref (basis-matrix-col-is bm) j))))
	  (declare (fixnum last-i))
	  (unless (= last-i pivot-i)
	    (setf (aref (aref (basis-matrix-col-is bm) j) pivot-index) last-i))))))
  ;; move pivots in queue to buckets
  (dotimes (k (length (basis-matrix-pivot-i-queue bm)))
    (let ((i (aref (basis-matrix-pivot-i-queue bm) k))
	  (j (aref (basis-matrix-pivot-j-queue bm) k))
	  (ci (aref (basis-matrix-pivot-ci-queue bm) k)))
      (pivot-add bm i j ci))))
	       
	       



;;;; Select the best pivot in the matrix
(defun basis-matrix-select-pivot (bm)
  (let ((m (basis-matrix-size bm))
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
		     (u (aref (basis-matrix-u-columns bm) j))
		     (vi (aref (hsv-vis u) ci)))
		(if (zerop vi)
		    (setf is-zero t
			  zero-i i
			  zero-j j
			  zero-ci ci)
		    (return-from basis-matrix-select-pivot
		      (values nil i j ci)))))
	  bucket))
     (basis-matrix-pivot-buckets bm))
    (values is-zero zero-i zero-j zero-ci)))
    


;;;; Perform a pivot in the matrix
(defun basis-matrix-perform-pivot (bm pk)
  (multiple-value-bind (is-zero i j ci)
      (basis-matrix-select-pivot bm)
    (cond (is-zero
	   (basis-matrix-row-is-redundant bm i)
	   (values -1 -1 -1))
	  (t
	   (pivot-count-update bm i j)
	   (if (basis-matrix-is-singular bm)
	       (values -1 -1 -1)
	       (let ((perm-i (aref (basis-matrix-i->pi bm) i))
		     (perm-j (aref (basis-matrix-j->pj bm) j))
		     (swap-i (aref (basis-matrix-pi->i bm) pk))
		     (swap-j (aref (basis-matrix-pj->j bm) pk)))
		 (setf (aref (basis-matrix-i->pi bm) i) pk
		       (aref (basis-matrix-j->pj bm) j) pk
		       (aref (basis-matrix-pi->i bm) pk) i
		       (aref (basis-matrix-pj->j bm) pk) j
		       (aref (basis-matrix-i->pi bm) swap-i) perm-i
		       (aref (basis-matrix-j->pj bm) swap-j) perm-j
		       (aref (basis-matrix-pi->i bm) perm-i) swap-i
		       (aref (basis-matrix-pj->j bm) perm-j) swap-j)
		 (values i j ci)))))))



(defun is-in-buckets (bm i j)
  (map-markowitz-tree 
     #'(lambda (mc bucket)
	 (map-pivot-bucket
	  #'(lambda (pivot-key ci)
	      (declare (ignore ci))
	      (let* ((bi (floor pivot-key (basis-matrix-size bm)))
		     (bj (- pivot-key (* i (basis-matrix-size bm)))))
		(when (and (= i bi) (= j bj))
		  (return-from is-in-buckets mc))))
	  bucket))
     (basis-matrix-pivot-buckets bm))
  nil)


(defun check-buckets (bm)
  (map-markowitz-tree 
     #'(lambda (mc bucket)
	 (map-pivot-bucket
	  #'(lambda (pivot-key ci)
	    (declare (ignore ci))
	      (let* ((i (floor pivot-key (basis-matrix-size bm)))
		     (j (- pivot-key (* i (basis-matrix-size bm))))
		     (rc (aref (basis-matrix-row-nnz bm) i))
		     (cc (aref (basis-matrix-col-nnz bm) j)))
		(assert (< 0 rc most-positive-fixnum))
		(assert (< 0 cc most-positive-fixnum))
		(assert (= mc (* (- cc 1) (- rc 1))))))
	  bucket))
     (basis-matrix-pivot-buckets bm)))
	       
  
(defun print-buckets (bm)
  (format t "buckets: ---~%")
  (map-markowitz-tree 
   #'(lambda (mc bucket)
       (format t "bucket ~A:~%" mc)
       (map-pivot-bucket
	#'(lambda (pivot-key ci)
	    (declare (ignore ci))
	    (let* ((i (floor pivot-key (basis-matrix-size bm)))
		   (j (- pivot-key (* i (basis-matrix-size bm))))
		   (rc (aref (basis-matrix-row-nnz bm) i))
		   (cc (aref (basis-matrix-col-nnz bm) j)))
	      (assert (< 0 rc most-positive-fixnum))
	      (assert (< 0 cc most-positive-fixnum))
	      (assert (= mc (* (- rc 1) (- cc 1))))
	      (format t "(~A , ~A) ~A ?= ~A ~%" i j mc (* (- rc 1) (- cc 1)))))
	bucket))
   (basis-matrix-pivot-buckets bm))
  (format t "~%"))
