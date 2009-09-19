;;;;; LU factorization of the basis
;;;;;



;;;; Permutates column k with column p >= k if |p| > ppivot-coef |k|
;;;; returns -1 is matrix is singular (this should not happen)
(defun lu-partial-pivoting (b pk in-block block-pj-end)
  (let ((current-abs-val 0)
	(largest-abs-val 0)
	(k (aref (basis-pj->j b) pk))
	(i (aref (basis-pi->i b) pk))
	(pivot-index -1)
	(pivot-j -1))
    ;; check current column
    (let* ((l (aref (basis-l-columns b) k))
	   (l-is (lu-eta-matrix-is l))
	   (l-n-nz (length l-is)))
      (dotimes (index l-n-nz)
	(when (= pk (aref (basis-i->pi b) (aref l-is index)))
	  (let ((abs-val (abs (lu-eta-get-ratio l index))))
	    (unless (zerop abs-val)
	      (setf current-abs-val abs-val
		    pivot-index index
		    pivot-j k)))
	  (return))))
    ;; check the others in the block
    (when in-block
      (loop for pj from (+ pk 1) upto block-pj-end do
	   (let ((j (aref (basis-pj->j b) pj)))
	     (when (<= pk (aref (basis-spikes b) j))
	       (let* ((block-l (aref (basis-l-columns b) j))
		      (block-l-is (lu-eta-matrix-is block-l))
		      (block-l-n-nz (length block-l-is))
		      (col-index -1)
		      (abs-val -1))
		 (dotimes (index block-l-n-nz)
		   (when (= i (aref block-l-is index))
		     (setf col-index index)
		     (setf abs-val (abs (lu-eta-get-ratio block-l index)))
		     (return)))
		 (when (and (/= -1 col-index)
			  (or (= -1 pivot-j)
			      (< largest-abs-val abs-val)))
		   (setf pivot-j j
			 pivot-index col-index
			 largest-abs-val abs-val)))))))
    ;; perform the partial pivoting test
    (unless (= -1 pivot-j)
      (unless (or (= k pivot-j)
		  (= 0.0 (basis-ppivot-coef b))
		  (> (* current-abs-val (basis-ppivot-coef b)) largest-abs-val))
	(setf (aref (basis-lu-ppivots b) pk) (aref (basis-j->pj b) pivot-j))
	(basis-permutate-columns b pivot-j k))
      t)))
		   


;;;; Stores the locations of bump blocks
(defun lu-identify-blocks (b)
  (let ((block-s -1)
	(m (basis-size b)))
    ;; go from last to first column in permuted matrix
    (loop for pj from (- m 1) downto 0 do
	 (let* ((j (aref (basis-pj->j b) pj))
		(s (aref (basis-spikes b) j)))
	   ;; store current column permutation matrix entries
	   (setf (aref (basis-lu-ppivots b) pj) pj)
	   (setf (aref (basis-col-counts b) j) pj)
	   (setf (aref (basis-row-counts b) pj) j)
	   ;; find bump blocks
	   (when (= -1 s)
	     (setf s (setf (aref (basis-spikes b) j) pj)))
	   (cond ((= pj block-s)
		  ;; start of block
		  (setf block-s -1)
		  (vector-push-extend pj (basis-lu-blocks b)))
		 ((and (= -1 block-s)
		       (/= s pj))
		  ;; end of block
		  (setf block-s s)
		  (vector-push-extend pj (basis-lu-blocks b))))))))
	


;;;; Splits column in lower and upper parts
(defun lu-split-column (b j pj)
  (let* ((l (aref (basis-l-columns b) j))
	 (u (make-lu-eta-matrix))
	 (l-is (lu-eta-matrix-is l))
	 (l-vis (lu-eta-matrix-vis l))
	 (l-vfs (lu-eta-matrix-vfs l))
	 (index-remove (basis-pip-spikes b))
	 (l-n-nz (length l-is)))
    (setf (lu-eta-matrix-j u) j
	  (lu-eta-matrix-coef u) (lu-eta-matrix-coef l))
    (vector-push-extend (aref (basis-pi->i b) pj)    (lu-eta-matrix-is u))
    (vector-push-extend (/ 1 (lu-eta-matrix-coef l)) (lu-eta-matrix-vis u))
    (vector-push-extend 1.0                          (lu-eta-matrix-vfs u))
    (setf (fill-pointer index-remove) 0)
    (dotimes (index l-n-nz)
      (let ((i  (aref l-is index))
	    (vi (aref l-vis index)))
	(if (zerop vi)
	    (vector-push-extend index index-remove)
	    (when (< (aref (basis-i->pi b) i) pj)
	      (vector-push-extend index index-remove)
	      (vector-push-extend i (lu-eta-matrix-is u))
	      (vector-push-extend (aref l-vis index) (lu-eta-matrix-vis u))
	      (vector-push-extend (aref l-vfs index) (lu-eta-matrix-vfs u))))))
    (lu-eta-normalize u)
    (vector-push-extend u (basis-u-columns b))
    (loop
       (when (zerop (length index-remove))
	 (return))
       (let ((index (vector-pop index-remove)))
	  (lu-eta-remove l index)))))



;;;; Makes lower ETA matrix
(defun lu-update-l (l i)
  (let ((l-is (lu-eta-matrix-is l))
	(l-vis (lu-eta-matrix-vis l))
	(l-vfs (lu-eta-matrix-vfs l))
	(l-n-nz (length (lu-eta-matrix-is l)))
	(pivot-index -1)
	(n (numerator (lu-eta-matrix-coef l)))
	(d (denominator (lu-eta-matrix-coef l))))
    ;; find pivot element
    (dotimes (index l-n-nz)
      (when (= i (aref l-is index))
	(setf pivot-index index)
	(return)))
    (assert (/= -1 pivot-index))
    (lu-eta-select-pivot l pivot-index)
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
(defun lu-prepare-pivot (b pivot-l pivot-i)
  (let ((refs (basis-refs b))
	(m (basis-size b)))
    (loop for index from pivot-i below m do
	 (setf (aref refs index) -1))
    (let* ((l-is (lu-eta-matrix-is pivot-l))
	   (l-n-nz (length l-is)))
      (dotimes (index l-n-nz)
	(setf (aref refs (aref (basis-i->pi b) (aref l-is index))) 
	      index)))))
	      

  
;;;;
(defun lu-perform-pivot (b pivot-l pivot-pj j)
  (let ((l (aref (basis-l-columns b) j))
	(m (basis-size b))
	(pivot-l-vis (lu-eta-matrix-vis pivot-l))
	(index-j-pivot-i -1)
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
      (dotimes (index l-n-nz)
	(when (= pivot-pj (aref (basis-i->pi b) (aref l-is index)))
	  (setf index-j-pivot-i index)
	  (return)))
      (unless (= -1 index-j-pivot-i)
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
	(loop for ip from (+ pivot-pj 1) below m do
	     (let ((i (aref (basis-pi->i b) ip))
		   (index-i-pivot-j (aref refs ip)))
	       (unless (or (= -1 index-i-pivot-j)
			   (aref flags ip))
		 ;; perform fill-in
		 (let ((val (* (aref l-vis index-j-pivot-i)
			       (aref pivot-l-vis index-i-pivot-j))))
		   (unless (zerop val)
		     (vector-push-extend i   l-is)
		     (vector-push-extend val l-vis)
		     (vector-push-extend 0.0 l-vfs))))))
	;; update row pivot-pj (is nonzero)
	(mulf (aref l-vis index-j-pivot-i) (aref pivot-l-vis 0))
	;; normalize j
	(lu-eta-normalize l)))))


    
;;;; LU decomposition
;;;; returns t on success, nil if basis matrix is singular
(defun lu-decomposition (b)
  (lu-identify-blocks b)
  (let ((block-pj-start -1)
	(block-pj-end -1)
	(in-block t)
	(success t))
    (dotimes (pj (basis-size b))
      ;; see if we are in a bump block
      (when (and in-block
		 (< block-pj-end pj))
	(setf in-block nil)
	(unless (zerop (length (basis-lu-blocks b)))
	  (setf block-pj-start (vector-pop (basis-lu-blocks b)))
	  (setf block-pj-end   (vector-pop (basis-lu-blocks b)))))
      (when (= pj block-pj-start)
	(setf in-block t))
      ;; partial pivoting and singularity check
      (unless (lu-partial-pivoting b pj in-block block-pj-end)
	(setf success nil)
	(return))
      ;; decompose the matrix
      ;; process current column
      (let* ((j (aref (basis-pj->j b) pj))
	     (l (aref (basis-l-columns b) j)))
	;; if the column is a spike, split it in two
	(unless (= pj (aref (basis-spikes b) j))
	  (lu-split-column b j pj))
	;; update l
	(lu-update-l l (aref (basis-pi->i b) pj))
	;; propagate l in bump block
	(when in-block
	  (lu-prepare-pivot b l pj)
	  (loop for block-pj from (+ pj 1) upto block-pj-end do
	       (let ((block-j (aref (basis-pj->j b) block-pj)))
		 (when (>= pj (aref (basis-spikes b) block-j))
		   (lu-perform-pivot b l pj block-j)))))))
    ;; reset original column permutation matrix
    (dotimes (k (basis-size b) success)
      (setf (aref (basis-j->pj b) k) (aref (basis-col-counts b) k)
	    (aref (basis-pj->j b) k) (aref (basis-row-counts b) k)))))

		   




;;;; Verifies the LU decomposition
(defun lu-check (b ob)
  (assert (= (basis-size b) (basis-size ob)))
  (let* ((m (basis-size b))
	 (ua (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (la (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (da (make-array (list m m) :initial-element 0 :element-type 'rational)))
    ;; fill the arrays
    (dotimes (j m)
      (let ((ol (aref (basis-l-columns ob) j))
	    (l (aref (basis-l-columns b) j)))
	(dotimes (index (length (lu-eta-matrix-is l)))
	  (let ((i (aref (lu-eta-matrix-is l) index)))
	    (setf (aref la (aref (basis-i->pi b) i) (aref (basis-j->pj b) j))
		  (* (lu-eta-matrix-coef l) (aref (lu-eta-matrix-vis l) index)))))
	(dotimes (index (length (lu-eta-matrix-is ol)))
	  (let ((i (aref (lu-eta-matrix-is ol) index)))
	    (setf (aref da (aref (basis-i->pi b) i) (aref (basis-j->pj b) j))
		  (* (lu-eta-matrix-coef ol) (aref (lu-eta-matrix-vis ol) index)))))))
    (dotimes (index (length (basis-u-columns b)))
      (let* ((u (aref (basis-u-columns b) index))
	     (j (lu-eta-matrix-j u)))
	(dotimes (r (length (lu-eta-matrix-is u)))
	  (let ((i (aref (lu-eta-matrix-is u) r)))
	    (setf (aref ua (aref (basis-i->pi b) i) (aref (basis-j->pj b) j))
		  (* (lu-eta-matrix-coef u) (aref (lu-eta-matrix-vis u) r)))))))
;    (print-2d-array da)
    ;; go ahead
    (dotimes (k m)
      ;; permutate columns 
      (let ((pp (aref (basis-lu-ppivots b) k)))
	(unless (= pp k)
	  (dotimes (i m)
	    (rotatef (aref la i pp) (aref la i k))
	    (rotatef (aref ua i pp) (aref ua i k))
	    (rotatef (aref da i pp) (aref da i k)))))
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
;    (print-2d-array da)
;    (print-2d-array ua)
;    (print-2d-array la)
    ;; check for equality
    (dotimes (i m t)
      (dotimes (j m)
	(assert (= (aref da i j) (aref ua i j)))))))
	  
    
      
		      
    
    
      
      
	    
      
      
    
