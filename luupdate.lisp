

;;;;
(defun lu-update-replace-column (bm pivot-j spike)
  (let* ((u (aref (basis-matrix-u-columns bm) pivot-j))
	 (m (basis-matrix-size bm))
	 (mus (basis-matrix-update-row-vals bm))
	 (u-seq (aref (basis-matrix-u-seqs bm) pivot-j))
	 (lastk (- m 1))
	 (i->pi (basis-matrix-i->pi bm))
	 (j->pj (basis-matrix-j->pj bm))
	 (pi->i (basis-matrix-pi->i bm))
	 (pj->j (basis-matrix-pj->j bm))
	 (pivot-i (aref (hsv-is u) (aref u-seq (- (hsv-length u) 1))))
	 (pivot-k (aref j->pj pivot-j))
	 (pivot-ci (if (zerop (hsv-length spike)) 
		       -1
		       (hsv-find pivot-i spike))))
    ;; fill column
    (reset-hsv u)
    (cond ((zerop (hsv-length spike))
	   (setf pivot-ci 0)
	   (hsv-add pivot-i 0 u))
	  ((/= -1 pivot-ci)
	   (copy-hsv-into-hsv spike u))
	  ((< (aref (hsv-is spike) (- (hsv-length spike) 1)) pivot-i)
	   (copy-hsv-into-hsv spike u)
	   (setf pivot-ci (hsv-length u))
	   (hsv-add pivot-i 0 u))
	  ((< pivot-i (aref (hsv-is spike) 0))
	   (hsv-add pivot-i 0 u)
	   (setf pivot-ci 0
		 (hsv-coef u) (hsv-coef spike))
	   (dotimes (ci (hsv-length spike))
	     (hsv-add (aref (hsv-is spike) ci) (aref (hsv-vis spike) ci) u)))
	  (t 
	   (setf (hsv-coef u) (hsv-coef spike))
	   (dotimes (ci (hsv-length spike))
	     (let ((i (aref (hsv-is spike) ci)))
	       (hsv-add i (aref (hsv-vis spike) ci) u)
	       (when (and (< i pivot-i)
			  (< pivot-i (aref (hsv-is spike) (+ ci 1))))
		 (setf pivot-ci (+ ci 1))
		 (hsv-add pivot-i 0 u))))))
    (assert (/= -1 pivot-ci))
    (assert (= pivot-i (aref (hsv-is u) pivot-ci)))
    (dotimes (ci (hsv-length u))
      (setf (aref u-seq ci) ci))
    ;; fill update values row and update u values
    (loop for k from 0 below m
       do (setf (aref mus k) 0))
    (loop for k from (aref j->pj pivot-j) below lastk
       do (let* ((pj (+ k 1))
		 (j (aref pj->j pj))
		 (uj (aref (basis-matrix-u-columns bm) j))
		 (ci (hsv-find pivot-i uj)))
	    (setf (aref mus k) 
		  (if (= -1 ci)
		      0
		      (- (* (hsv-coef uj) (aref (hsv-vis uj) ci)))))
	    (unless (= -1 ci)
	      (setf (aref (hsv-vis uj) ci) 0))))
    (setf (aref mus lastk) (* (hsv-coef u) (aref (hsv-vis u) pivot-ci)))
    ;; update permutation matrices
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
    ;; update sequence for new column
    (sort-increasing-bounded (aref (basis-matrix-u-seqs bm) pivot-j) 
			     (hsv-length u)
			     #'(lambda (k) (aref i->pi (aref (hsv-is u) k))))
    ;; return
    (values pivot-ci pivot-k)))



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
		 (u-row-seq (aref (basis-matrix-u-seqs bm) j))
		 (u (aref (basis-matrix-u-columns bm) j)))
	    ;; get current coef and residue
	    (loop for k from (- (hsv-length u) 1) downto 0
	       do (let* ((ci (aref u-row-seq k))
			 (i (aref (hsv-is u) ci))
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
		  (assert (= rhs (* (hsv-coef u) ukk)))
		  (when (zerop ukk)
		    (when (zerop residue)
		      (print u))
		    (assert (not (zerop residue))))
		  (setf (aref mus pj)
			(+ residue ukk)))
		(progn
		  (setf (aref mus pj) (/ (- (/ rhs (hsv-coef u)) residue) ukk))
		  (unless (zerop (aref mus pj))
		    (let ((denomj (denominator (aref mus pj))))
		      (if (zerop cdenom)
			  (setf cdenom denomj)
			  (mulf cdenom (/ denomj (gcd cdenom denomj))))))))))
    cdenom))



;;;;
(defun lu-update-build-eta (bm pivot-k pivot-j pivot-ci cdenom)
  (let* ((m (basis-matrix-size bm))
	 (l (aref (basis-matrix-l-file bm) (basis-matrix-n-l-file bm)))
	 (pi->i (basis-matrix-pi->i bm))
	 (pivot-col (aref (basis-matrix-u-columns bm) pivot-j))
	 (mus (basis-matrix-update-row-vals bm)))
    ;; build eta column
    (if (zerop cdenom)
	(loop for k from pivot-k below (- m 1)
	   do (assert (zerop (aref mus k))))
	(progn 
	  (assert (loop for k from pivot-k below (- m 1)
		     unless (zerop (aref mus k))
		     do (return t)))
	  (loop for k from pivot-k below (- m 1)
	     unless (zerop (aref mus k))
	     do (assert (= (denominator (aref mus k)) 
			   (gcd (denominator (aref mus k)) cdenom))))))
    (unless (zerop cdenom)
      (reset-hsv l)
      (setf (hsv-coef l) (/ 1 cdenom))
      (loop for k from pivot-k below (- m 1)
	 unless (zerop (aref mus k))
	 do (hsv-add (aref pi->i k) (* cdenom (aref mus k)) l))
      (hsv-normalize l)
      (when (< 1 (hsv-length l))
	(hsv-sort-indices-increasing l))
      (setf (aref (basis-matrix-l-pivot-file bm) (basis-matrix-n-l-file bm)) 
	    (aref pi->i (- m 1)))
      (incf (basis-matrix-n-l-file bm)))
    ;; update new column in u
    (let ((pivot-denom (denominator (aref mus (- m 1)))))
      (divf (hsv-coef pivot-col) pivot-denom)
      (dotimes (ci (hsv-length pivot-col))
	(if (= ci pivot-ci)
	    (setf (aref (hsv-vis pivot-col) ci) (numerator (aref mus (- m 1))))
	    (mulf (aref (hsv-vis pivot-col) ci) pivot-denom))))
    (hsv-normalize pivot-col)))
    
   
    
;;;;
(defun lu-update (bm j spike)
  (multiple-value-bind (pivot-ci pivot-k)
      (lu-update-replace-column bm j spike)
    (let ((cdenom (lu-update-compute-values bm pivot-k)))
      (lu-update-build-eta bm pivot-k j pivot-ci cdenom))))
  
    
