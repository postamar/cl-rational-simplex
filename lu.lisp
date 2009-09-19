;;;;; LU factorization of the basis
;;;;;



;;;; Splits column in lower and upper parts
(defun lu-split-column (b j pj pivot-ci)
  (let* ((l (aref (basis-l-columns b) j))
	 (u (make-lu-eta-matrix))
	 (l-is (lu-eta-matrix-is l))
	 (l-vis (lu-eta-matrix-vis l))
	 (l-vfs (lu-eta-matrix-vfs l))
	 (counter 0)
	 (l-n-nz (length l-is)))
    (lu-eta-select-pivot l pivot-ci)
    (setf (lu-eta-matrix-j u) j
	  (lu-eta-matrix-coef u) (lu-eta-matrix-coef l))
    (vector-push-extend (aref (basis-pi->i b) pj)    (lu-eta-matrix-is u))
    (vector-push-extend (/ 1 (lu-eta-matrix-coef l)) (lu-eta-matrix-vis u))
    (vector-push-extend 1.0                          (lu-eta-matrix-vfs u))
    (dotimes (index l-n-nz)
      (let ((i  (aref l-is index))
	    (vi (aref l-vis index)))
	(if (zerop vi)
	    (progn (setf (aref (basis-refs b) counter) index) 
		   (incf counter))
	    (when (< (aref (basis-i->pi b) i) pj)
	      (setf (aref (basis-refs b) counter) index) 
	      (incf counter)
	      (vector-push-extend i (lu-eta-matrix-is u))
	      (vector-push-extend (aref l-vis index) (lu-eta-matrix-vis u))
	      (vector-push-extend (aref l-vfs index) (lu-eta-matrix-vfs u))))))
    (lu-eta-normalize u)
    (vector-push-extend u (basis-u-columns b))
    (unless (zerop counter)
      (loop 
	 (let ((index (aref (basis-refs b) (decf counter))))
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
  (let ((refs (basis-refs b))
	(m (basis-size b)))
    (loop for qi from k below m
       do (setf (aref refs qi) -1))
    (let* ((l-is (lu-eta-matrix-is pivot-l))
	   (l-n-nz (length l-is)))
      (dotimes (index l-n-nz)
	(setf (aref refs (aref (basis-i->pi b) (aref l-is index))) 
	      index)))))
	      

  
;;;;
(defun lu-update-right-l (b pivot-l pivot-pj j index-j-pivot-i)
  (let ((l (aref (basis-l-columns b) j))
	(m (basis-size b))
	(pivot-l-vis (lu-eta-matrix-vis pivot-l))
	(flags (basis-flags b))
	(refs  (basis-refs b)))
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
	       (ip (aref (basis-i->pi b) i)))
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
	 do (let ((i (aref (basis-pi->i b) ip))
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
		    (incf (aref (basis-row-nnz b) i))
		    (vector-push-extend j    (aref (basis-row-js b) i))
		    (vector-push-extend ci   (aref (basis-row-cis b) i))
		    (incf (aref (basis-col-nnz b) j))
		    (vector-push-extend i    (aref (basis-col-is b) j))
		    (check-pivot-heap-elements b)
		    (let ((href (pivot-heap-append b j i ci)))
		      (vector-push-extend href (aref (basis-row-hrefs b) i))
		      (vector-push-extend href (aref (basis-col-hrefs b) j))))))))
      ;; update row pivot-pj (is nonzero)
      (mulf (aref l-vis index-j-pivot-i) (aref pivot-l-vis 0))
      ;; normalize j
      (lu-eta-normalize l)
      ;; update value changes in heap
      (let* ((hrefs   (aref (basis-col-hrefs b) j))
	     (n-hrefs (length hrefs)))
	(dotimes (index n-hrefs)
	  (pivot-heap-update b (aref hrefs index))))
      )))




;;;; LU decomposition
;;;; returns t on success, nil if basis matrix is singular
(defun basis-lu-decomposition (b)
  (dotimes (k (basis-size b) t)
    (check-pivot-heap b)
    (format t " start ~A~%" k)
    (print-pivot-heap b)
    ;; perform pivot
    (multiple-value-bind (pivot-j pivot-i pivot-ci)
	(basis-perform-pivot b k) 
      (when (= -1 pivot-j)
	(return))
      (assert (not (basis-is-singular b)))
      (assert (= pivot-i (aref (basis-pi->i b) k)))
      (assert (= pivot-j (aref (basis-pj->j b) k)))
      (let ((pivot-l (aref (basis-l-columns b) pivot-j))
	    (pivot-row-hrefs (aref (basis-row-hrefs b) pivot-i))
	    (pivot-row-js    (aref (basis-row-js b)    pivot-i))
	    (pivot-row-cis   (aref (basis-row-cis b)   pivot-i))
	    (pivot-col-hrefs (aref (basis-col-hrefs b) pivot-j))
	    (pivot-col-is    (aref (basis-col-is b)    pivot-j)))
	(print-pivot-heap b)
	(check-pivot-heap b)
	(format t " mid~A~%" k)
	;; split column if necessary
	(when (> k (aref (basis-spikes b) pivot-j))
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
		(assert (= pivot-i (aref (lu-eta-matrix-is (aref (basis-l-columns b) j)) ci)))
		(lu-update-right-l b pivot-l k j ci))))) 
	(check-pivot-heap b)
	(format t " end~A~%" k)
	;; update spike counters
	(dotimes (index (length pivot-row-js))
	  (let ((j (aref pivot-row-js index)))
	    (minf (aref (basis-spikes b) j) k)))
	;; empty pivot row and pivot column
	(setf (fill-pointer pivot-row-hrefs) 0
	      (fill-pointer pivot-row-js)    0
	      (fill-pointer pivot-row-cis)   0
	      (fill-pointer pivot-col-hrefs) 0
	      (fill-pointer pivot-col-is)    0)))))
	    



;;;; Verifies the LU decomposition
(defun lu-check (b ob)
  (assert (= (basis-size b) (basis-size ob)))
  (let* ((m (basis-size b))
	 (j->pj (make-array m :initial-element 0 :element-type 'fixnum))
	 (pj->j (make-array m :initial-element 0 :element-type 'fixnum))
	 (pa (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (ua (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (la (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (da (make-array (list m m) :initial-element 0 :element-type 'rational)))
    ;; fill the arrays
    (dotimes (j m)
      (setf (aref j->pj j) j
	    (aref pj->j j) j))
    (dotimes (j m)
      (let ((ol (aref (basis-l-columns ob) j))
	    (l (aref (basis-l-columns b) j)))
	(dotimes (index (length (lu-eta-matrix-is l)))
	  (let ((i (aref (lu-eta-matrix-is l) index)))
	    (setf (aref la (aref (basis-i->pi b) i) (aref (basis-j->pj b) j))
		  (* (lu-eta-matrix-coef l) (aref (lu-eta-matrix-vis l) index)))))
	(dotimes (index (length (lu-eta-matrix-is ol)))
	  (let ((i (aref (lu-eta-matrix-is ol) index)))
	    (setf (aref pa (aref (basis-i->pi b) i) j)
		  (* (lu-eta-matrix-coef ol) (aref (lu-eta-matrix-vis ol) index)))
	    (setf (aref da (aref (basis-i->pi b) i) j)
		  (* (lu-eta-matrix-coef ol) (aref (lu-eta-matrix-vis ol) index)))))))
    (dotimes (index (length (basis-u-columns b)))
      (let* ((u (aref (basis-u-columns b) index))
	     (j (lu-eta-matrix-j u)))
	(dotimes (r (length (lu-eta-matrix-is u)))
	  (let ((i (aref (lu-eta-matrix-is u) r)))
	    (setf (aref ua (aref (basis-i->pi b) i) (aref (basis-j->pj b) j))
		  (* (lu-eta-matrix-coef u) (aref (lu-eta-matrix-vis u) r)))))))
   ; (print-2d-array da)
    ;; go ahead
    (dotimes (k m)
      ;; permutate columns 
      (let ((j1 (aref (basis-lu-ppivots1 b) k))
	    (j2 (aref (basis-lu-ppivots2 b) k)))
	(let ((pj1 (aref j->pj j1))
	      (pj2 (aref j->pj j2)))
	  (unless (= j1 j2)
	    (dotimes (i m)
	      (rotatef (aref pa i pj1) (aref pa i pj2))
	      (rotatef (aref j->pj j1) (aref j->pj j2))
	      (rotatef (aref pj->j pj1) (aref pj->j pj2))
	      (rotatef (aref da i pj1) (aref da i pj2))))))
      ;; left-multiply by pivot
      ;; rows i > k
      (loop for i from (+ k 1) below m do
	   (dotimes (j m)
	     (incf (aref da i j)
		   (* (aref la i k)
		      (aref da k j)))))
      ;; row k
      (dotimes (j m)
	(mulf (aref da k j) (aref la k k))))
    ;; add diagonal 1s
    (dotimes (k m)
      (setf (aref ua k k) 1))
    ;; print arrays
    ;(print-2d-array da)
    ;(print-2d-array pa)
    ;(print-2d-array ua)
    ;(print-2d-array la)
    ;; check for equality
    (dotimes (i m t)
      (dotimes (j m)
	(assert (= (aref da i j) (aref ua i j)))))))
	  
