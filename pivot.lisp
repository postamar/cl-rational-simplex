(in-package :rationalsimplex)

;;;;; Dynamic Markowitz pivot implementation
;;;;;
;;;;; During factorization, pivot elements are selected
;;;;; using the Markowitz criterion: in the residual matrix,
;;;;; select the non-zero element which minimizes (nr-1).(nc-1)
;;;;; where nr and nc are the non-zero counts in the row and
;;;;; the column of the pivot element, respectively.
;;;;; The implementation is slightly complicated by the fact that 
;;;;; fill-in occurs during the factorization.




;;;; Auxilliary function for pivot-find
(defun markowitz-col (bm j col-nnz)
  (let ((min-i -1)
	(min-ci -1)
	(min-mc most-positive-fixnum))
    (dotimes (ki col-nnz (values min-i min-ci min-mc))
      (let* ((i (aref (aref (basis-matrix-col-is bm) j) ki))
	     (ri (find-index-bounded (aref (basis-matrix-row-js bm) i) 
				     (aref (basis-matrix-row-nnz bm) i)
				     j))
	     (ci (aref (aref (basis-matrix-row-cis bm) i) ri))
	     (u (aref (basis-matrix-u-columns bm) j))
	     (val (aref (hsv-vis u) ci))
	     (mc (if (zerop val) 
		     (* (basis-matrix-size bm) (basis-matrix-size bm))
		     (* (- col-nnz 1) 
			(- (aref (basis-matrix-row-nnz bm) i) 1)))))
	(when (< mc min-mc)
	  (setf min-i i
		min-ci ci
		min-mc mc))))))



;;;; Auxilliary function for pivot-find
(defun markowitz-row (bm i row-nnz)
  (let ((min-j -1)
	(min-ci -1)
	(min-mc most-positive-fixnum))
    (dotimes (kj row-nnz (values min-j min-ci min-mc))
      (let* ((j (aref (aref (basis-matrix-row-js bm) i) kj))
	     (ci (aref (aref (basis-matrix-row-cis bm) i) kj))
	     (u (aref (basis-matrix-u-columns bm) j))
	     (val (aref (hsv-vis u) ci))
	     (mc (if (zerop val) 
		     (* (basis-matrix-size bm) (basis-matrix-size bm))
		     (* (- row-nnz 1) 
			(- (aref (basis-matrix-col-nnz bm) j) 1)))))
	(when (< mc min-mc)
	  (setf min-j j
		min-ci ci
		min-mc mc))))))



;;;; Selects the best pivot element in the residual matrix
(defun pivot-find (bm)
  (let* ((m (basis-matrix-size bm))
	 (m*m (* m m))
	 (mc-min m*m)
	 (n 0)
	 (pivot-i -1)
	 (pivot-j -1)
	 (pivot-ci -1))
    
    (cond 
      ;; error checking
      ((/= 0 (aref (basis-matrix-col-bucket-sizes bm) 0))
       (error "zero weight col bucket should be empty"))
      ((/= 0 (aref (basis-matrix-row-bucket-sizes bm) 0))
       (error "zero weight row bucket should be empty"))
      ;; select singleton pivots if possible,
      ;; i.e. those in the buckets of weight 1
      ((< 0 (aref (basis-matrix-col-bucket-sizes bm) 1))
       (let ((j (aref (aref (basis-matrix-col-buckets bm) 1)
		      (- (aref (basis-matrix-col-bucket-sizes bm) 1) 1))))
	 (multiple-value-bind (i ci mc)
	     (markowitz-col bm j 1)
	   (values (= mc m*m) i j ci))))
      ((< 0 (aref (basis-matrix-row-bucket-sizes bm) 1))
       (let* ((i (aref (aref (basis-matrix-row-buckets bm) 1)
		       (- (aref (basis-matrix-row-bucket-sizes bm) 1) 1))))
	 (multiple-value-bind (j ci mc)
	     (markowitz-row bm i 1)
	   (values (= mc m*m) i j ci))))
      ;; select the best pivot element, relative to selection criterion
      (t	 
       (loop named find-pivot-loop for bucket-weight from 2 upto m
	  do (let ((col-bsize (aref (basis-matrix-col-bucket-sizes bm) bucket-weight))
		   (row-bsize (aref (basis-matrix-row-bucket-sizes bm) bucket-weight))
		   (col-bucket (aref (basis-matrix-col-buckets bm) bucket-weight))
		   (row-bucket (aref (basis-matrix-row-buckets bm) bucket-weight))
		   (col-mc-bound (* (- bucket-weight 1) (- bucket-weight 1)))
		   (row-mc-bound (* bucket-weight (- bucket-weight 1))))
	       (unless (zerop col-bsize)
		 (dotimes (kj col-bsize)
		   (let ((j (aref col-bucket kj)))
		     (multiple-value-bind (i ci mc)
			 (markowitz-col bm j bucket-weight)
			   (when (< mc mc-min)
			     (setf pivot-i i
				   pivot-j j
				   pivot-ci ci
				   mc-min mc)
			     (when (<= mc-min col-mc-bound)
			       (return-from find-pivot-loop)))
			   (when (and (<= (basis-matrix-row-col-max bm) (incf n)) 
				      (< mc-min m*m))
			     (return-from find-pivot-loop))
			   (unless (zerop row-bsize)
			     (dotimes (ki row-bsize)
			       (let ((i (aref row-bucket ki)))
				 (multiple-value-bind (j ci mc)
				     (markowitz-row bm i bucket-weight)
				   (when (< mc mc-min)
				     (setf pivot-i i
					   pivot-j j
					   pivot-ci ci
					   mc-min mc)
				     (when (<= mc-min row-mc-bound)
				       (return-from find-pivot-loop)))
				   (when (and (<= (basis-matrix-row-col-max bm) (incf n))
					      (< mc-min m*m))
				     (return-from find-pivot-loop))))))))))))
       ;; return (T if success, row index, column index, index in column-hsv)
       (values (= mc-min m*m) pivot-i pivot-j pivot-ci)))))
		   
			       


;;;; Adds a new pivot (from fill-in)
(defun pivot-add (bm i j ci)
  ;; update buckets
  (let* ((pivot-col-nnz (aref (basis-matrix-col-nnz bm) j))
	 (pivot-col-bucket (aref (basis-matrix-col-buckets bm) pivot-col-nnz))
	 (pivot-row-nnz (aref (basis-matrix-row-nnz bm) i))
	 (pivot-row-bucket (aref (basis-matrix-row-buckets bm) pivot-row-nnz))
	 (col-bucket-k (dotimes (k (aref (basis-matrix-col-bucket-sizes bm) pivot-col-nnz) 
				 (error "pivot col not found in bucket"))
			 (when (= j (aref pivot-col-bucket k))
			   (return k))))
	 (row-bucket-k (dotimes (k (aref (basis-matrix-row-bucket-sizes bm) pivot-row-nnz)
				 (error "pivot row not found in bucket"))
			 (when (= i (aref pivot-row-bucket k))
			   (return k)))))
    (setf (aref pivot-col-bucket col-bucket-k)
	  (aref pivot-col-bucket (decf (aref (basis-matrix-col-bucket-sizes bm) pivot-col-nnz))))
    (setf (aref pivot-row-bucket row-bucket-k)
	  (aref pivot-row-bucket (decf (aref (basis-matrix-row-bucket-sizes bm) pivot-row-nnz))))
    (incf pivot-col-nnz)
    (setf (aref (aref (basis-matrix-col-buckets bm) pivot-col-nnz)
		(aref (basis-matrix-col-bucket-sizes bm) pivot-col-nnz))
	  j)
    (incf (aref (basis-matrix-col-bucket-sizes bm) pivot-col-nnz))
    (incf pivot-row-nnz)
    (setf (aref (aref (basis-matrix-row-buckets bm) pivot-row-nnz)
		(aref (basis-matrix-row-bucket-sizes bm) pivot-row-nnz))
	  i)
    (incf (aref (basis-matrix-row-bucket-sizes bm) pivot-row-nnz)))
  ;; update residual matrix
  (let ((col-is (aref (basis-matrix-col-is bm) j))
	(k (aref (basis-matrix-col-nnz bm) j)))
    (setf (aref col-is k) i)
    (incf (aref (basis-matrix-col-nnz bm) j))
    (loop 
       (when (or (zerop k)
		 (< (aref col-is (- k 1)) (aref col-is k)))
	 (return))
       (rotatef (aref col-is (- k 1)) (aref col-is k))
       (decf k)))
  (let ((row-js (aref (basis-matrix-row-js bm) i))
	(row-cis (aref (basis-matrix-row-cis bm) i))
	(k (aref (basis-matrix-row-nnz bm) i)))
    (setf (aref row-js k) j
	  (aref row-cis k) ci)
    (incf (aref (basis-matrix-row-nnz bm) i))
    (loop
       (when (or (zerop k)
		 (< (aref row-js (- k 1)) (aref row-js k)))
	 (return))
       (rotatef (aref row-js (- k 1)) (aref row-js k))
       (rotatef (aref row-cis (- k 1)) (aref row-cis k))
       (decf k))))


			       
;;;; Updates pivot buckets and residual matrix after pivot selection and removal
(defun pivot-count-update (bm pivot-i pivot-j)
  (let ((pivot-row-nnz (aref (basis-matrix-row-nnz bm) pivot-i))
	(pivot-col-nnz (aref (basis-matrix-col-nnz bm) pivot-j))
	(pivot-col-is (aref (basis-matrix-col-is bm) pivot-j))
	(pivot-row-js (aref (basis-matrix-row-js bm) pivot-i)))
    ;; update on pivot row
    (dotimes (kj pivot-row-nnz)
      (let* ((j (aref pivot-row-js kj))
	     (col-nnz (aref (basis-matrix-col-nnz bm) j))
	     (new-col-nnz (- col-nnz 1))
	     (current-bucket (aref (basis-matrix-col-buckets bm) col-nnz))
	     (current-bucket-k (dotimes (k (aref (basis-matrix-col-bucket-sizes bm) col-nnz) 
					 (error "col ~A not found in bucket ~A" j col-nnz))
				 (when (= j (aref current-bucket k))
				   (return k)))))
	;; remove column from bucket
	(decf (aref (basis-matrix-col-bucket-sizes bm) col-nnz))
	(setf (aref current-bucket current-bucket-k)
	      (aref current-bucket (aref (basis-matrix-col-bucket-sizes bm) col-nnz)))
	;; add column in new bucket
	(unless (= j pivot-j)
	  (when (zerop col-nnz)
	    (basis-matrix-column-is-redundant bm j)
	    (return-from pivot-count-update))
	  (setf (aref (aref (basis-matrix-col-buckets bm) new-col-nnz)
		      (aref (basis-matrix-col-bucket-sizes bm) new-col-nnz))
		j)
	  (incf (aref (basis-matrix-col-bucket-sizes bm) new-col-nnz))
	  ;; update residual matrix
	  (let* ((col-is (aref (basis-matrix-col-is bm) j))
		 (pivot-ki (find-index-bounded col-is col-nnz pivot-i)))
	    (decf (aref (basis-matrix-col-nnz bm) j))
	    (assert (/= -1 pivot-ki))
	    (loop for ki from pivot-ki below (aref (basis-matrix-col-nnz bm) j)
	       do (setf (aref col-is ki) (aref col-is (+ ki 1))))))))
    ;; update on pivot column
    (dotimes (ki pivot-col-nnz)
      (let* ((i (aref pivot-col-is ki))
	     (row-nnz (aref (basis-matrix-row-nnz bm) i))
	     (new-row-nnz (- row-nnz 1))
	     (current-bucket (aref (basis-matrix-row-buckets bm) row-nnz))
	     (current-bucket-k (dotimes (k (aref (basis-matrix-row-bucket-sizes bm) row-nnz)
					 (error "row not found in bucket"))
				 (when (= i (aref current-bucket k))
				   (return k)))))
	;; remove row from bucket
	(decf (aref (basis-matrix-row-bucket-sizes bm) row-nnz))
	(setf (aref current-bucket current-bucket-k)
	      (aref current-bucket (aref (basis-matrix-row-bucket-sizes bm) row-nnz)))
	;; add row in new bucket
	(unless (= i pivot-i)
	  (when (zerop row-nnz)
	    (basis-matrix-row-is-redundant bm i)
	    (return-from pivot-count-update))
	  (setf (aref (aref (basis-matrix-row-buckets bm) new-row-nnz)
		      (aref (basis-matrix-row-bucket-sizes bm) new-row-nnz))
		i)
	  (incf (aref (basis-matrix-row-bucket-sizes bm) new-row-nnz))
	  ;; update residual matrix
	  (let* ((row-js (aref (basis-matrix-row-js bm) i))
		 (row-cis (aref (basis-matrix-row-cis bm) i))
		 (pivot-kj (find-index-bounded row-js row-nnz pivot-j)))
	    (decf (aref (basis-matrix-row-nnz bm) i))
	    (loop for kj from pivot-kj below (aref (basis-matrix-row-nnz bm) i)
	       do (setf (aref row-js kj) (aref row-js (+ kj 1))
			(aref row-cis kj) (aref row-cis (+ kj 1))))))))
    ;; remove pivot row and col from residual matrix
    (setf (aref (basis-matrix-row-nnz bm) pivot-i) 0
	  (aref (basis-matrix-col-nnz bm) pivot-j) 0)))

		 

;;;; Perform a pivot in the matrix
(defun basis-matrix-perform-pivot (bm pk)
  (multiple-value-bind (is-zero i j ci)
      (pivot-find bm)
    (let ((pivot-row-nnz (aref (basis-matrix-row-nnz bm) i)))
      (cond (is-zero
	     (basis-matrix-row-is-redundant bm i)
	     (values -1 -1 -1 -1))
	    (t
	     (pivot-count-update bm i j)
	     (if (basis-matrix-is-singular bm)
		 (values -1 -1 -1 -1)
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
		   (values i j ci pivot-row-nnz))))))))
			 

	       
  

