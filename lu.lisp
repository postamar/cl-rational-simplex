;;;;; LU factorization of the basis
;;;;;



;;;; makes new l-eta matrix, if necessary
(defun lu-split-column (bm pivot-i pivot-j k pivot-ci)
  (let* ((i->pi (basis-matrix-i->pi bm))
	 (pi->i (basis-matrix-pi->i bm))
	 (u (aref (basis-matrix-u-columns bm) pivot-j))
	 (u-is (lu-eta-matrix-is u))
	 (u-vis (lu-eta-matrix-vis u))
	 (u-vfs (lu-eta-matrix-vfs u))
	 (us (aref (basis-matrix-u-file bm) k))
	 (us-is (lu-eta-matrix-is us))
	 (us-vis (lu-eta-matrix-vis us))
	 (us-vfs (lu-eta-matrix-vfs us))
	 (l nil))
    ;; set coef and pivot in u eta matrix 
    (setf (lu-eta-matrix-coef us) (lu-eta-matrix-coef u))
    (vector-push-extend (aref i->pi (aref u-is pivot-ci)) us-is)
    (vector-push-extend (aref u-vis pivot-ci) us-vis)
    (vector-push-extend 0.0 us-vfs)
    ;; split column in l and u eta matrices
    (dotimes (ci (length u-is))
      (let* ((i (aref u-is ci))
	     (ip (aref i->pi i))
	     (vi (aref u-vis ci)))
	(cond ((zerop vi))
	      ((= ip k)
	       (assert (= ci pivot-ci)))
	      ((< ip k)
	       (vector-push-extend ip us-is)
	       (vector-push-extend vi us-vis)
	       (vector-push-extend 0.0 us-vfs))
	      ((not l)
	       (setf l (aref (basis-matrix-l-file bm) (basis-matrix-n-l-file bm)))
	       (incf (basis-matrix-n-l-file bm))
	       (setf (lu-eta-matrix-j l) pivot-i)
	       (setf (lu-eta-matrix-coef l) (/ 1 (aref u-vis pivot-ci)))
	       (vector-push-extend pivot-i (lu-eta-matrix-is l))
	       (vector-push-extend (aref u-vis pivot-ci) (lu-eta-matrix-vis l))
	       (vector-push-extend 0.0 (lu-eta-matrix-vfs l))
	       (vector-push-extend i (lu-eta-matrix-is l))
	       (vector-push-extend (- vi) (lu-eta-matrix-vis l))
	       (vector-push-extend 0.0 (lu-eta-matrix-vfs l)))
	      (t
	       (vector-push-extend i (lu-eta-matrix-is l))
	       (vector-push-extend (- vi) (lu-eta-matrix-vis l))
	       (vector-push-extend 0.0 (lu-eta-matrix-vfs l))))))
    ;; normalize eta matrices
    (lu-eta-normalize us)
    (when l
      (lu-eta-normalize l)
      (lu-eta-sort-indices-increasing l))
    ;; refill U column
    (reset-lu-eta-matrix u)
    (setf (lu-eta-matrix-coef u) (lu-eta-matrix-coef us))
    (dotimes (k (length (lu-eta-matrix-is us)))
      (vector-push (aref pi->i (aref (lu-eta-matrix-is us) k)) u-is)
      (vector-push (aref (lu-eta-matrix-vis us) k) u-vis)
      (vector-push (aref (lu-eta-matrix-vfs us) k) u-vfs))
    ;; sort u eta matrix
    (lu-eta-sort-indices-increasing us)
    (if l 
	(- (basis-matrix-n-l-file bm) 1)
	-1)))


;;;;
(defun lu-prepare-update (b pivot-l k)
  (let ((refs (basis-matrix-refs b))
	(m (basis-matrix-size b)))
    (loop for qi from k below m
       do (setf (aref refs qi) -1))
    (let* ((l-is (lu-eta-matrix-is pivot-l))
	   (l-n-nz (length l-is)))
      (dotimes (index l-n-nz)
	(setf (aref refs (aref (basis-matrix-i->pi b) (aref l-is index))) 
	      index)))))
	      

  


;;;;
(defun lu-update-right (bm pivot-l pivot-pj j index-j-pivot-i)
  (let* ((u (aref (basis-matrix-u-columns bm) j))
	 (m (basis-matrix-size bm))
	 (pivot-l-vis (lu-eta-matrix-vis pivot-l))
	 (flags (basis-matrix-flags bm))
	 (refs  (basis-matrix-refs bm))
	 (pivot-coef (lu-eta-matrix-coef pivot-l))
	 (pivot-n (numerator pivot-coef))
	 (pivot-d (denominator pivot-coef))
	 (u-is (lu-eta-matrix-is u))
	 (u-vis (lu-eta-matrix-vis u))
	 (u-vfs (lu-eta-matrix-vfs u))
	 (u-n-nz (length u-is)))
    (loop for ip from pivot-pj below m do
	 (setf (aref flags ip) nil))
    ;; update coef on j
    (divf (lu-eta-matrix-coef u) pivot-d) 
    ;; do rows with nonzeros
    (dotimes (index u-n-nz)
      (let* ((i (aref u-is index))
	     (ip (aref (basis-matrix-i->pi bm) i)))
	;; update rows pi /= pivot-pj (is nonzero)
	(unless (= pivot-pj ip)
	  (mulf (aref u-vis index) pivot-d))
	;; update rows pi > pivot-pj with nonzeros
	(when (< pivot-pj ip)
	  (setf (aref flags ip) t)
	  (let ((index-i-pivot-j (aref refs ip)))
	    (unless (= -1 index-i-pivot-j)
	      (incf (aref u-vis index) 
		    (* pivot-n
		       (aref u-vis index-j-pivot-i)
		       (aref pivot-l-vis index-i-pivot-j))))))))
    ;; do fill-in on rows pi > pivot-pj
    (loop for ip from (+ pivot-pj 1) below m
       do (let ((i (aref (basis-matrix-pi->i bm) ip))
		(ci (length u-is))
		(index-i-pivot-j (aref refs ip)))
	    (unless (or (= -1 index-i-pivot-j)
			(aref flags ip))
	      ;; perform fill-in
	      (let ((val (* pivot-n
			    (aref u-vis index-j-pivot-i)
			    (aref pivot-l-vis index-i-pivot-j))))
		(unless (zerop val)
		  (vector-push-extend i   u-is)
		  (vector-push-extend val u-vis)
		  (vector-push-extend 0.0 u-vfs)
		  (pivot-add-new bm i j ci))))))
    ;; update row pivot-pj (is nonzero)
    (mulf (aref u-vis index-j-pivot-i) pivot-d)
    ;; normalize j
    (lu-eta-normalize u)))



;;;; LU decomposition
;;;; returns t on success, nil if basis matrix is singular
(defun basis-matrix-lu-decomposition (bm)
  (dotimes (k (basis-matrix-size bm) t)
    ;; perform pivot
    (multiple-value-bind (pivot-i pivot-j pivot-ci)
	(basis-matrix-perform-pivot bm k) 
      (when (= -1 pivot-j)
	(return))
      (assert (not (basis-matrix-is-singular bm)))
      (assert (= pivot-i (aref (basis-matrix-pi->i bm) k)))
      (assert (= pivot-j (aref (basis-matrix-pj->j bm) k)))
      ;; make l eta matrix, if necessary
      (let ((l-file-index (lu-split-column bm pivot-i pivot-j k pivot-ci)))
	(unless (= -1 l-file-index)
	  ;; update remaining columns
	  (let* ((pivot-l       (aref (basis-matrix-l-file bm)  l-file-index))
		 (pivot-row-js  (aref (basis-matrix-row-js bm)  pivot-i))
		 (pivot-row-cis (aref (basis-matrix-row-cis bm) pivot-i))
		 (pivot-col-is  (aref (basis-matrix-col-is bm)  pivot-j)))
	    (unless (= 1 (length pivot-row-js))
	      (lu-prepare-update bm pivot-l k)
	      (dotimes (index (length pivot-row-js))
		(let ((j (aref pivot-row-js index))
		      (ci (aref pivot-row-cis index)))
		  (unless (= j pivot-j)
		    (lu-update-right bm pivot-l k j ci)))))
	    ;; empty pivot row and pivot column
	    (setf (fill-pointer pivot-row-js)    0
		  (fill-pointer pivot-row-cis)   0
		  (fill-pointer pivot-col-is)    0)))))))
	    
	    
	    
	

#|

	;; split column if necessary
	(when (> k (aref (basis-matrix-spikes b) pivot-j))
	  (lu-split-column b pivot-j k pivot-ci)
	  (setf pivot-ci 0))
	;; update l
	(lu-update-u pivot-u pivot-i pivot-ci)
	;; update columns to the right of the pivot if necessary
	(unless (= 1 (length pivot-row-js))
	  (lu-prepare-update b pivot-l k)
	  (dotimes (index (length pivot-row-js))
	    (let ((j (aref pivot-row-js index))
		  (ci (aref pivot-row-cis index)))
	      (unless (= j pivot-j)
		(assert (= pivot-i (aref (lu-eta-matrix-is 
					  (aref (basis-matrix-l-columns b) j)) ci)))
		(lu-update-right-l b pivot-l k j ci)))))
	;; update spike counters
	(dotimes (index (length pivot-row-js))
	  (let ((j (aref pivot-row-js index)))
	    (minf (aref (basis-matrix-spikes b) j) k)))
	    
|#

