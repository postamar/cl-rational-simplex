;;;;; LU factorization of the basis
;;;;;



;;;; makes new l-eta matrix, if necessary
(defun lu-split-column (bm pivot-i pivot-j k pivot-ci)
  (let* ((i->pi (basis-matrix-i->pi bm))
	 (fill-in-counter (aref (basis-matrix-fill-ins bm) pivot-j))
	 (u-seq (aref (basis-matrix-u-seqs bm) pivot-j))
	 (u (aref (basis-matrix-u-columns bm) pivot-j))
	 (u-is (hsv-is u))
	 (u-vis (hsv-vis u))
	 (l nil)
	 (lf nil))
    ;; re-order fill-in
    (dotimes (fill-in-k fill-in-counter)
      (let ((ci (+ (- (length u-is) fill-in-counter) fill-in-k)))
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
	       (incf (basis-matrix-n-l-file bm))
	       (setf (hsv-j l) pivot-i)
	       (setf (hsv-coef l) (/ 1 (aref u-vis pivot-ci)))
	       (vector-push-extend pivot-i (hsv-is l))
	       (vector-push-extend (aref u-vis pivot-ci) (hsv-vis l))
	       (vector-push-extend i (hsv-is l))
	       (vector-push-extend (- vi) (hsv-vis l))
	       (setf (aref u-vis ci) 0))
	      (t
	       (vector-push-extend i (hsv-is l))
	       (vector-push-extend (- vi) (hsv-vis l))
	       (setf (aref u-vis ci) 0)))))
    ;; rebuild l and u
    (when l
      (hsv-normalize l)
      (copy-hsv-into-hsv-float l lf))
    (hsv-remove-zeros u)
    (hsv-normalize u)
    (copy-hsv-into-hsv-float u (aref (basis-matrix-uf-columns bm) pivot-j))
    ;; build permuted sequence
    (dotimes (ci (hsv-length u))
      (vector-push-extend ci u-seq))
    (sort u-seq #'< :key #'(lambda (k) (aref i->pi (aref u-is k))))
    ;; return l index 
    (if l 
	(- (basis-matrix-n-l-file bm) 1)
	-1)))
    


;;;;
(defun lu-prepare-update (b pivot-l k)
  (let ((refs (basis-matrix-refs b))
	(m (basis-matrix-size b)))
    (loop for qi from k below m
       do (setf (aref refs qi) -1))
    (let* ((l-is (hsv-is pivot-l))
	   (l-n-nz (length l-is)))
      (dotimes (index l-n-nz)
	(setf (aref refs (aref (basis-matrix-i->pi b) (aref l-is index))) 
	      index)))))
	      

  


;;;;
(defun lu-update-right (bm pivot-l pivot-pj j index-j-pivot-i)
  (let* ((u (aref (basis-matrix-u-columns bm) j))
	 (m (basis-matrix-size bm))
	 (pivot-l-vis (hsv-vis pivot-l))
	 (flags (basis-matrix-flags bm))
	 (refs  (basis-matrix-refs bm))
	 (pivot-coef (hsv-coef pivot-l))
	 (pivot-n (numerator pivot-coef))
	 (pivot-d (denominator pivot-coef))
	 (u-is (hsv-is u))
	 (u-vis (hsv-vis u))
	 (u-n-nz (length u-is)))
    (loop for ip from pivot-pj below m do
	 (setf (aref flags ip) nil))
    ;; update coef on j
    (divf (hsv-coef u) pivot-d) 
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
		  (incf (aref (basis-matrix-fill-ins bm) j))
		  (pivot-add-new bm i j ci))))))
    ;; update row pivot-pj (is nonzero)
    (mulf (aref u-vis index-j-pivot-i) pivot-d)
    ;; normalize j
    (hsv-normalize u)))




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
	    
	    
	    
	


