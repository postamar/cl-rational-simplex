;;;;; LU factorization of the basis
;;;;;


;;;; Makes upper triangular eta matrices
(defun make-upper-triangular-etas (b) 
  (let ((nu (length (basis-matrix-u-columns b)))
	(i->pi (basis-matrix-i->pi b))
	(j->pj (basis-matrix-j->pj b)))
    (dotimes (k nu t)
      (let* ((u (aref (basis-matrix-u-columns b) k))
	     (u-is (lu-eta-matrix-is u))
	     (u-vis (lu-eta-matrix-vis u))
	     (u-vfs (lu-eta-matrix-vfs u))
	     (n (length u-is))
	     (us (make-lu-eta-matrix))
	     (us-is (lu-eta-matrix-is us))
	     (us-vis (lu-eta-matrix-vis us))
	     (us-vfs (lu-eta-matrix-vfs us)))
	(setf (lu-eta-matrix-coef us) (lu-eta-matrix-coef u)
	      (lu-eta-matrix-j us) (aref j->pj (lu-eta-matrix-j u)))
	(dotimes (i n)
	  (vector-push-extend (aref i->pi (aref u-is i)) us-is)
	  (vector-push-extend (aref u-vis i) us-vis)
	  (vector-push-extend (aref u-vfs i) us-vfs))
	(lu-eta-sort-indices-increasing us)
	(vector-push-extend us (basis-matrix-u-file b))))))
      


;;;; Splits column in lower and upper parts
(defun lu-split-column (b j pj pivot-ci)
  (let* ((l (aref (basis-matrix-l-columns b) j))
	 (u (make-lu-eta-matrix))
	 (l-is (lu-eta-matrix-is l))
	 (l-vis (lu-eta-matrix-vis l))
	 (l-vfs (lu-eta-matrix-vfs l))
	 (counter 0)
	 (l-n-nz (length l-is)))
    (lu-eta-select-pivot l pivot-ci)
    (setf (lu-eta-matrix-j u) j
	  (lu-eta-matrix-coef u) (lu-eta-matrix-coef l))
    (vector-push-extend (aref (basis-matrix-pi->i b) pj) (lu-eta-matrix-is u))
    (vector-push-extend (/ 1 (lu-eta-matrix-coef l)) (lu-eta-matrix-vis u))
    (vector-push-extend 1.0                          (lu-eta-matrix-vfs u))
    (dotimes (index l-n-nz)
      (let ((i  (aref l-is index))
	    (vi (aref l-vis index)))
	(if (zerop vi)
	    (progn (setf (aref (basis-matrix-refs b) counter) index) 
		   (incf counter))
	    (when (< (aref (basis-matrix-i->pi b) i) pj)
	      (setf (aref (basis-matrix-refs b) counter) index) 
	      (incf counter)
	      (vector-push-extend i (lu-eta-matrix-is u))
	      (vector-push-extend (aref l-vis index) (lu-eta-matrix-vis u))
	      (vector-push-extend (aref l-vfs index) (lu-eta-matrix-vfs u))))))
    (lu-eta-normalize u)
    (vector-push-extend u (basis-matrix-u-columns b))
    (unless (zerop counter)
      (loop 
	 (let ((index (aref (basis-matrix-refs b) (decf counter))))
	   (lu-eta-remove l index))
	 (when (zerop counter)
	   (return))))))


;;;; Makes lower ETA matrix
(defun lu-update-l (l pivot-i pivot-ci)
  (let ((l-is (lu-eta-matrix-is l))
	(l-vis (lu-eta-matrix-vis l))
	(l-vfs (lu-eta-matrix-vfs l))
	(l-n-nz (length (lu-eta-matrix-is l)))
	(n (numerator (lu-eta-matrix-coef l)))
	(d (denominator (lu-eta-matrix-coef l))))
    (assert (< pivot-ci l-n-nz))
    (assert (= pivot-i (aref l-is pivot-ci)))
    (lu-eta-select-pivot l pivot-ci)
    ;; update values
    (let ((vgcd d))
      (loop for index from 1 below l-n-nz do
	   (let ((vi (aref l-vis index)))
	     (setf vgcd (gcd vgcd vi))))
      (loop for index from 1 below l-n-nz do
	   (mulf (aref l-vis index) (/ (- n) vgcd)))
      (setf (lu-eta-matrix-coef l) (/ vgcd (* n (aref l-vis 0)))
	    (aref l-vis 0)         (/ d vgcd))
      (loop for index from 0 below l-n-nz do
	   (setf (aref l-vfs index) 
		 (float (* (lu-eta-matrix-coef l) (aref l-vis index))))))))



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
(defun lu-update-right-l (b pivot-l pivot-pj j index-j-pivot-i)
  (let ((l (aref (basis-matrix-l-columns b) j))
	(m (basis-matrix-size b))
	(pivot-l-vis (lu-eta-matrix-vis pivot-l))
	(flags (basis-matrix-flags b))
	(refs  (basis-matrix-refs b)))
    (loop for ip from pivot-pj below m do
	 (setf (aref flags ip) nil))
    (let* ((pivot-coef (lu-eta-matrix-coef pivot-l))
	   (pivot-n (numerator pivot-coef))
	   (pivot-d (denominator pivot-coef))
	   (l-is (lu-eta-matrix-is l))
	   (l-vis (lu-eta-matrix-vis l))
	   (l-vfs (lu-eta-matrix-vfs l))
	   (l-n-nz (length l-is)))
      ;; update coef on j
      (divf (lu-eta-matrix-coef l) pivot-d) 
      ;; update row pivot-pj (is nonzero)
      (mulf (aref l-vis index-j-pivot-i) pivot-n)
      ;; do rows with nonzeros
      (dotimes (index l-n-nz)
	(let* ((i (aref l-is index))
	       (ip (aref (basis-matrix-i->pi b) i)))
	  ;; update rows pi /= pivot-pj (is nonzero)
	  (unless (= pivot-pj ip)
	    (mulf (aref l-vis index) pivot-d))
	    ;; update rows pi > pivot-pj with nonzeros
	  (when (< pivot-pj ip)
	    (setf (aref flags ip) t)
	    (let ((index-i-pivot-j (aref refs ip)))
	      (unless (= -1 index-i-pivot-j)
		(incf (aref l-vis index) 
		      (* (aref l-vis index-j-pivot-i)
			 (aref pivot-l-vis index-i-pivot-j))))))))
      ;; do fill-in on rows pi > pivot-pj
      (loop for ip from (+ pivot-pj 1) below m
	 do (let ((i (aref (basis-matrix-pi->i b) ip))
		  (ci (length l-is))
		  (index-i-pivot-j (aref refs ip)))
	      (unless (or (= -1 index-i-pivot-j)
			  (aref flags ip))
		;; perform fill-in
		(let ((val (* (aref l-vis index-j-pivot-i)
			      (aref pivot-l-vis index-i-pivot-j))))
		  (unless (zerop val)
		    (vector-push-extend i   l-is)
		    (vector-push-extend val l-vis)
		    (vector-push-extend 0.0 l-vfs)
		    (pivot-add-new b i j ci))))))
		    
      ;; update row pivot-pj (is nonzero)
      (mulf (aref l-vis index-j-pivot-i) (aref pivot-l-vis 0))
      ;; normalize j
      (lu-eta-normalize l)
      )))




;;;; LU decomposition
;;;; returns t on success, nil if basis matrix is singular
(defun basis-matrix-lu-decomposition (b)
  (dotimes (k (basis-matrix-size b) (make-upper-triangular-etas b))
    ;; perform pivot
    (multiple-value-bind (pivot-i pivot-j pivot-ci)
	(basis-matrix-perform-pivot b k) 
      (when (= -1 pivot-j)
	(return))
      (assert (not (basis-matrix-is-singular b)))
      (assert (= pivot-i (aref (basis-matrix-pi->i b) k)))
      (assert (= pivot-j (aref (basis-matrix-pj->j b) k)))
      (let ((pivot-l (aref (basis-matrix-l-columns b) pivot-j))
	    (pivot-row-js    (aref (basis-matrix-row-js b)    pivot-i))
	    (pivot-row-cis   (aref (basis-matrix-row-cis b)   pivot-i))
	    (pivot-col-is    (aref (basis-matrix-col-is b)    pivot-j)))
	;; split column if necessary
	(when (> k (aref (basis-matrix-spikes b) pivot-j))
	  (lu-split-column b pivot-j k pivot-ci)
	  (setf pivot-ci 0))
	;; update l
	(lu-update-l pivot-l pivot-i pivot-ci)
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
	;; empty pivot row and pivot column
	(setf (fill-pointer pivot-row-js)    0
	      (fill-pointer pivot-row-cis)   0
	      (fill-pointer pivot-col-is)    0)))))
	    


