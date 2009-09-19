(in-package :rationalsimplex)

;;;;; LU factorization of the basis
;;;;; 
;;;;; Factorizes B from scratch, by gaussian elimination



;;;; Makes new L-eta matrix, if necessary
(defun lu-split-column (bm pivot-i pivot-j k pivot-ci)
  (let* ((i->pi (basis-matrix-i->pi bm))
	 (fill-in-counter (aref (basis-matrix-fill-ins bm) pivot-j))
	 (u (aref (basis-matrix-u-columns bm) pivot-j))
	 (u-is (hsv-is u))
	 (u-vis (hsv-vis u))
	 (u-seq (aref (basis-matrix-u-seqs bm) pivot-j))
	 (l nil)
	 (lf nil))
    ;; re-order fill-in
    (dotimes (fill-in-k fill-in-counter)
      (let ((ci (+ (- (hsv-length u) fill-in-counter) fill-in-k)))
	(loop
	   (cond ((zerop ci)
		  (return))
		 ((< (aref u-is ci) (aref u-is (- ci 1)))
		  (rotatef (aref u-is ci) (aref u-is (- ci 1)))
		  (rotatef (aref u-vis ci) (aref u-vis (- ci 1)))
		  (cond ((= ci pivot-ci)
			 (decf pivot-ci))
			((= (- ci 1) pivot-ci)
			 (incf pivot-ci)))
		  (decf ci))
		 (t
		  (return))))))
    ;; build l
    (dotimes (ci (hsv-length u))
      (let* ((i (aref u-is ci))
	     (ip (aref i->pi i))
	     (vi (aref u-vis ci)))
	(cond ((zerop vi))
	      ((< ip k))
	      ((= ip k)
	       (assert (= ci pivot-ci)))
	      ((not l)
	       (setf l  (aref (basis-matrix-l-file bm)  (basis-matrix-n-l-file bm))
		     lf (aref (basis-matrix-lf-file bm) (basis-matrix-n-l-file bm)))
	       (setf (aref (basis-matrix-l-pivot-file bm) (basis-matrix-n-l-file bm)) pivot-i)
	       (incf (basis-matrix-n-l-file bm))
	       (incf (basis-matrix-n-l-factor-file bm))
	       (setf (hsv-coef l) (/ 1 (aref u-vis pivot-ci)))
	       (hsv-add i (- vi) l)
	       (setf (aref u-vis ci) 0))
	      (t
	       (hsv-add i (- vi) l)
	       (setf (aref u-vis ci) 0)))))
    ;; rebuild l and u
    (when l
      (hsv-normalize l)
      (copy-hsv-into-hsv-float l lf))
    (hsv-remove-zeros u)
    (hsv-normalize u)
    (copy-hsv-into-hsv-float u (aref (basis-matrix-uf-columns bm) pivot-j))
    (dotimes (ci (hsv-length u))
      (setf (aref u-seq ci) ci))
    (sort-increasing-bounded u-seq (hsv-length u) 
			     #'(lambda (k) (aref i->pi (aref (hsv-is u) k))))
    ;; return l index 
    (if l 
	(- (basis-matrix-n-l-file bm) 1)
	-1)))
    


;;;; Auxilliary function for basis-matrix-lu-factorization
(defun lu-prepare-update (b pivot-l k)
  (let ((refs (basis-matrix-refs b))
	(m (basis-matrix-size b)))
    (loop for qi from k below m
       do (setf (aref refs qi) -1))
    (let* ((l-is (hsv-is pivot-l))
	   (l-n-nz (hsv-length pivot-l)))
      (dotimes (index l-n-nz)
	(setf (aref refs (aref (basis-matrix-i->pi b) (aref l-is index))) 
	      index)))))
	      

  


;;;; Updates column in residual matrix after pivot
(defun lu-update-right (bm pivot-l pivot-pj j index-j-pivot-i)
  (let* ((u (aref (basis-matrix-u-columns bm) j))
	 (m (basis-matrix-size bm))
	 (flags (basis-matrix-flags bm))
	 (refs  (basis-matrix-refs bm))
	 (pivot-coef (hsv-coef pivot-l))
	 (pivot-n (numerator pivot-coef))
	 (pivot-d (denominator pivot-coef)))
    (loop for ip from pivot-pj below m do
	 (setf (aref flags ip) nil))
    ;; update coef on j
    (divf (hsv-coef u) pivot-d) 
    ;; do rows with nonzeros
    (dotimes (index (hsv-length u))
      (let* ((i (aref (hsv-is u) index))
	     (ip (aref (basis-matrix-i->pi bm) i)))
	;; update rows pi /= pivot-pj (is nonzero)
	(unless (= pivot-pj ip)
	  (mulf (aref (hsv-vis u) index) pivot-d))
	;; update rows pi > pivot-pj with nonzeros
	(when (< pivot-pj ip)
	  (setf (aref flags ip) t)
	  (let ((index-i-pivot-j (aref refs ip)))
	    (unless (= -1 index-i-pivot-j)
	      (incf (aref (hsv-vis u) index) 
		    (* pivot-n
		       (aref (hsv-vis u) index-j-pivot-i)
		       (aref (hsv-vis pivot-l) index-i-pivot-j))))))))
    ;; do fill-in on rows pi > pivot-pj
    (loop for ip from (+ pivot-pj 1) below m
       do (let ((i (aref (basis-matrix-pi->i bm) ip))
		(index-i-pivot-j (aref refs ip)))
	    (unless (or (= -1 index-i-pivot-j)
			(aref flags ip))
	      ;; perform fill-in
	      (let ((val (* pivot-n
			    (aref (hsv-vis u) index-j-pivot-i)
			    (aref (hsv-vis pivot-l) index-i-pivot-j))))
		(unless (zerop val)
		  (pivot-add bm i j (hsv-length u))
		  (hsv-add i val u)
		  (incf (aref (basis-matrix-fill-ins bm) j)))))))
    ;; update row pivot-pj (is nonzero)
    (mulf (aref (hsv-vis u) index-j-pivot-i) pivot-d)
    ;; normalize j
    (hsv-normalize u)))




;;;; LU factorization
;;;; returns t on success, nil if basis matrix is singular
(defun basis-matrix-lu-factorization (bm)
  (dotimes (k (basis-matrix-size bm) t)
    ;; select and perform pivot
    (multiple-value-bind (pivot-i pivot-j pivot-ci pivot-row-nnz)
	(basis-matrix-perform-pivot bm k) 
      (when (= -1 pivot-j)
	(return)) ; pivot unsuccessful
      (assert (not (basis-matrix-is-singular bm)))
      (assert (= pivot-i (aref (basis-matrix-pi->i bm) k)))
      (assert (= pivot-j (aref (basis-matrix-pj->j bm) k)))
      ;; make L eta matrix, if necessary
      (let ((l-file-index (lu-split-column bm pivot-i pivot-j k pivot-ci)))
	(unless (= -1 l-file-index)
	  ;; update remaining columns
	  (let* ((pivot-l       (aref (basis-matrix-l-file bm)  l-file-index))
		 (pivot-row-js  (aref (basis-matrix-row-js bm)  pivot-i))
		 (pivot-row-cis (aref (basis-matrix-row-cis bm) pivot-i)))
	    (unless (= 1 pivot-row-nnz)
	      (lu-prepare-update bm pivot-l k)
	      (dotimes (index pivot-row-nnz)
		(let ((j (aref pivot-row-js index))
		      (ci (aref pivot-row-cis index)))
		  (unless (= j pivot-j)
		    (lu-update-right bm pivot-l k j ci)))))))))))
	    

