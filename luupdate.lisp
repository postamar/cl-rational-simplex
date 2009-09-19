

;;;;
(defun lu-update-replace-column (bm j lp col-ref)
  (let* ((u (aref (basis-matrix-u-columns bm) j))
	 (m (basis-matrix-size bm))
	 (col (column-hsv (adjvector-column-ref (lp-columns lp) col-ref)))
	 (u-row-seq (aref (basis-matrix-u-row-seqs bm) j))
	 (i->pi (basis-matrix-i->pi bm))
	 (j->pj (basis-matrix-j->pj bm))
	 (pivot-i (aref (hsv-is u) 
			(aref (basis-matrix-u-row-seqs bm) (- (hsv-length u) 1))))
	 (pivot-ci -1))
    ;; remove column from u row representation
    ;; add to update values row
    (dotimes (ci (hsv-length col))
      (let* ((i (aref (hsv-is u) ci))
	     (col-seq (aref (basis-matrix-u-col-seqs bm) i))
	     (ci-seq (aref (basis-matrix-u-ci-seqs bm) i))
	     (rowlength (aref (basis-matrix-u-row-lengths bm) i))
	     (ri (find-index-bounded col-seq rowlength j)))
	(unless (= -1 ri)
	  (loop for k from (+ ri 1) below rowlength
	     do (setf (aref col-seq (- k 1)) (aref col-seq k)
		      (aref ci-seq (- k 1)) (aref col-seq k)))
	  (decf (aref (basis-matrix-u-row-lengths bm) i)))))
    ;; fill column
    (reset-hsv u)
    (setf (hsv-coef u) (hsv-coef col))
    (dotimes (ci (hsv-length col))
      (let ((i (adjvector-fixnum-ref (lp-active-row-inds lp) (aref (hsv-is col) ci))))
	(hsv-add i (aref (hsv-vis col) ci) u)
	(setf (aref u-row-seq ci) ci)))
    (hsv-sort-indices-increasing u)
    (sort-increasing-bounded u-row-seq (hsv-length u) 
			     #'(lambda (k) (aref i->pi (aref (hsv-is u) k))))
    ;; add to u row representation
    (dotimes (ci (hsv-length u))
      (let* ((i (aref (hsv-is u) ci))
	     (col-seq (aref (basis-matrix-u-col-seqs bm) i))
	     (ci-seq (aref (basis-matrix-u-ci-seqs bm) i))
	     (rowlength (aref (basis-matrix-u-row-lengths bm) i))
	     (ri (find-closest-index-bounded col-seq rowlength j)))
	(assert (< (aref col-seq ri) j (aref col-seq (+ ri 1))))
	(when (= i pivot-i)
	  (setf pivot-ci ci))
	(loop for k from (aref (basis-matrix-u-row-lengths bm) i) above ri
	   do (setf (aref col-seq k) (aref col-seq (- k 1))
		    (aref ci-seq k) (aref col-seq (- k 1))))
	(setf (aref col-seq ri) j
	      (aref ci-seq ri) ci)
	(incf (aref (basis-matrix-u-row-lengths bm) i))))
    ;; fill update values row and update u values
    (let ((rowlength (aref (basis-matrix-u-row-lengths bm) pivot-i))
	  (mus (basis-matrix-update-row-vals bm))
	  (col-seq (aref (basis-matrix-u-col-seqs bm) pivot-i))
	  (ci-seq (aref (basis-matrix-u-ci-seqs bm) pivot-i)))
      (loop for k from (aref j->pj j) below m
	 do (setf (aref mus k) 0))
      (dotimes (ri rowlength)
	(let* ((uj (aref col-seq ri))
	       (ci (aref ci-seq ri))
	       (u (aref (basis-matrix-u-columns bm) uj)))
	  (assert (= pivot-i (aref (hsv-is u) ci)))
	  (setf (aref mus (aref j->pj uj))
		(* (hsv-coef u) (aref (hsv-vis u) ci)))
	  (setf (aref (hsv-vis u) ci) 0))))
    ;; return
    pivot-ci))



;;;;
(defun lu-update-permute-spike (bm pivot-j pivot-ci)
  (let* ((spike-col (aref (basis-matrix-u-columns bm) pivot-j))
	 (m (basis-matrix-size bm))
	 (pivot-i (aref (hsv-is spike-col) pivot-ci))
	 (lastk (- m 1))
	 (i->pi (basis-matrix-i->pi bm))
	 (j->pj (basis-matrix-j->pj bm))
	 (pi->i (basis-matrix-pi->i bm))
	 (pj->j (basis-matrix-pj->j bm))
	 (pivot-k (aref j->pj pivot-j)))
    (loop for k from pivot-k below lastk
       do (let ((next-i (aref pi->i (+ k 1)))
		(next-j (aref pj->j (+ k 1))))
	    (setf (aref pi->i k) next-i
		  (aref i->pi next-i) k
		  (aref pj->j k) next-j
		  (aref j->pj next-j) k)))
    (setf (aref pi->i lastk) pivot-i
	  (aref i->pi pivot-i) lastk
	  (aref pj->j lastk) pivot-j
	  (aref j->pj pivot-j) lastk)
    pivot-k))



;;;; 
(defun lu-update-compute-values (bm pivot-k)
  (let* ((m (basis-matrix-size bm))
	 (pj->j (basis-matrix-pj->j bm))
	 (i->pi (basis-matrix-i->pi bm))
	 (cdenom 0)
	 (mus (basis-matrix-update-row-vals bm)))
    ;; calculate the update row values
    (loop for pj from pivot-k below m
       do (let* ((rhs (aref mus pj))
		 (residue 0)
		 (ukk 0)
		 (j (aref pj->j pj))
		 (u (aref (basis-matrix-u-columns bm) j)))
	    ;; get current coef and residue
	    (loop for ci from (- (hsv-length u) 1) downto 0
	       do (let* ((i (aref (hsv-is u) ci))
			 (ip (aref i->pi i)))
		    (cond ((< ip pivot-k)
			   (return))
			  ((> pj ip)
			   (incf residue (* (aref (hsv-vis u) ci) (aref mus ip))))
			  ((= pj ip)
			   (setf ukk (aref (hsv-vis u) ci))))))
	    ;; compute the current value
	    (if (= pj (- m 1))
		(progn 
		  (assert (= rhs ukk))
		  (setf (aref mus pj)
			(+ residue ukk)))
		(progn
		  (setf (aref mus pj) (/ (- (/ rhs (hsv-coef u)) residue) ukk))
		  (let ((denomj (denominator (aref mus pj))))
		    (if (zerop cdenom)
			(setf cdenom denomj)
			(mulf cdenom (/ denomj (gcd cdenom denomj)))))))))
    cdenom))



;;;;
(defun lu-update-build-eta (bm pivot-k pivot-j pivot-ci cdenom)
  (let* ((m (basis-matrix-size bm))
	 (l (aref (basis-matrix-l-file bm) (basis-matrix-n-l-file bm)))
	 (pj->j (basis-matrix-pj->j bm))
	 (pivot-col (aref (basis-matrix-u-columns bm) pivot-j))
	 (mus (basis-matrix-update-row-vals bm)))
    ;; build eta column
    (reset-hsv l)
    (setf (hsv-coef l) (/ 1 cdenom))
    (hsv-add (aref pj->j (- m 1)) cdenom l)
    (loop for k from pivot-k below (- m 1)
       unless (zerop (aref mus k))
       do (hsv-add (aref pj->j k) (* cdenom (aref mus k)) l))
    (when (< 1 (hsv-length l))
      (hsv-sort-indices-increasing l)
      (hsv-normalize l)
      (setf (aref (basis-matrix-l-pivot-file bm) (basis-matrix-n-l-file bm)) 
	    (aref pj->j (- m 1)))
      (incf (basis-matrix-n-l-file bm)))
    ;; update new column in u
    (let ((pivot-denom (denominator (aref mus (- m 1)))))
      (dotimes (ci (hsv-length pivot-col))
	(if (= ci pivot-ci)
	    (setf (aref (hsv-vis pivot-col) ci) (numerator (aref mus (- m 1))))
	    (mulf (aref (hsv-vis pivot-col) ci) pivot-denom))))
    (hsv-normalize pivot-col)))
    
    
    
;;;;
(defun lu-update (bm j lp col-ref)
  (let ((pivot-ci (lu-update-replace-column bm j lp col-ref)))
    (let ((pivot-k (lu-update-permute-spike bm j pivot-ci)))
      (let ((cdenom (lu-update-compute-values bm pivot-k)))
	(lu-update-build-eta bm pivot-k j pivot-ci cdenom)))))

  
    
