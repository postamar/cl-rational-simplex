
;;;;
(defun simplex-compute-reduced-costs (sd)
  (let* ((b (simplex-basis sd))
	 (bm (basis-matrix b))
	 (m (basis-matrix-size bm))
	 (lp (simplex-lp sd))
	 (bh (basis-header b))
	 (rcosts (basis-reduced-costs b))
	 (flags (basis-column-flags b))
	 (n (length rcosts))
	 (cdenom 1))
    ;; load basic variable costs into input vector
    (dotimes (k m)
      (let* ((col (adjvector-column-ref (lp-columns lp) (aref bh k)))
	     (cd (denominator (column-c col))))
	(mulf cdenom (/ cd (gcd cd cdenom)))))
    (setf (hsv-coef (simplex-hsv sd)) (/ 1 cdenom)
	  (hsv-length (simplex-hsv sd)) 0)
    (dotimes (k m)
      (let* ((col (adjvector-column-ref (lp-columns lp) (aref bh k)))
	     (c (column-c col)))
	(assert (= (denominator c) (gcd (denominator c) cdenom)))
	(assert (integerp  (* (numerator c) (/ cdenom (denominator c)))))
	(hsv-add k (* (numerator c) (/ cdenom (denominator c))) (simplex-hsv sd))))
    ;; compute multipliers
    (btran (simplex-btran sd) (simplex-hsv sd))
    ;; compute reduced costs
    (dotimes (j n)
      (setf (aref rcosts j)
	    (if (eq (aref flags j) 'basic)
		0
		(lp-get-cost lp j))))
    (compute-pivot-row sd)
    (dotimes (k (simplex-pivot-row-length sd))
      (let ((j (aref (simplex-pivot-row-col-refs sd) k)))
	(unless (eq (aref flags j) 'basic)
	  (decf (aref rcosts j) 
		(aref (simplex-pivot-row-values sd) k)))))
    (setf (simplex-pivot-row-length sd) 0)))




;;;;
(defun simplex-compute-primal-values (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b))
	 (tr (simplex-ftran sd))
	 (m (basis-matrix-size bm))
	 (flags (basis-column-flags b))
	 (rhs (make-array m :initial-element 0 :element-type 'rational))
	 (cdenom 1))
    ;; compute right-hand side in solving Bb = rhs for b
    (dotimes (i m)
      (let* ((v 0)
	     (row-ref (adjvector-fixnum-ref (lp-active-row-refs lp) i))
	     (row (adjvector-row-ref (lp-rows lp) row-ref)))
	(dotimes (k (adjvector-fixnum-fill-pointer (row-col-refs row)))
	  (let* ((col-ref (adjvector-fixnum-ref (row-col-refs row) k))
		 (flag (aref flags col-ref))
		 (col (adjvector-column-ref (lp-columns lp) col-ref))
		 (a (rational-in-column col (adjvector-fixnum-ref (row-col-indices row) k))))
	    (cond ((eq flag 'nonbasic-lower-bound)
		   (assert (column-has-l col))
		   (decf v (* a (column-l col))))
		  ((eq flag 'nonbasic-upper-bound)
		   (assert (column-has-u col))
		   (decf v (* a (column-u col)))))))
	(setf (aref rhs i) v)
	(mulf cdenom (/ (denominator v) (gcd (denominator v) cdenom)))))
    ;; set rhs for solving Bb = rhs for b
    (setf (hsv-coef (simplex-hsv sd)) (/ 1 cdenom)
	  (hsv-length (simplex-hsv sd)) 0)
    (dotimes (i m)
      (setf (aref (basis-primal-values b) i) 0)
      (let ((v (aref rhs i)))
	(unless (zerop v)
	  (hsv-add i (* (numerator v) (/ cdenom (denominator v))) (simplex-hsv sd)))))
    (ftran tr (simplex-hsv sd))
    (check-ftran sd tr (simplex-hsv sd))
    ;; fill primal value array
    (dotimes (k (hsv-length (tran-hsv tr)))
      (let ((i (aref (hsv-is (tran-hsv tr)) k))
	    (x (* (hsv-coef (tran-hsv tr)) (aref (hsv-vis (tran-hsv tr)) k))))
	(setf (aref (basis-primal-values b) i) x)))))


;;;;    
(defun simplex-prepare-phase2 (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b))
	 (m (basis-matrix-size bm))
	 (flags (basis-column-flags b))
	 (z 0))
    (setf (basis-in-phase1 b) nil)
    ;; compute reduced costs
    (simplex-compute-reduced-costs sd)
    ;; set nonbasic variable flags and compute objective value
    (dotimes (j (adjvector-column-fill-pointer (lp-columns lp)))
      (let* ((col (adjvector-column-ref (lp-columns lp) j))
	     (c (* (- (lp-obj-sense lp)) (column-c col)))
	     (d (aref (basis-reduced-costs b) j))
	     (flag (aref flags j)))
	(cond ((eq flag 'basic))
	      ((eq flag 'nonbasic-lower-bound)
	       (assert (column-has-l col))
	       (incf z (* c (column-l col))))
	      ((eq flag 'nonbasic-upper-bound)
	       (assert (column-has-u col))
	       (incf z (* c (column-u col))))
	      ;; from now on, only inactive or boxed
	      ((< 0 d)
	       (assert (column-has-l col))
	       (setf (aref flags j) 'nonbasic-lower-bound)
	       (incf z (* c (column-l col))))
	      ((< d 0)
	       (assert (column-has-u col))
	       (setf (aref flags j) 'nonbasic-upper-bound)
	       (incf z (* c (column-u col))))
	      ((and (column-has-l col) (not (column-has-u col)))
	       (setf (aref flags j) 'nonbasic-lower-bound)
	       (incf z (* c (column-l col))))
	      ((and (not (column-has-l col)) (column-has-u col))
	       (setf (aref flags j) 'nonbasic-upper-bound)
	       (incf z (* c (column-u col))))
	      ((<= 0 c)
	       (setf (aref flags j) 'nonbasic-lower-bound)
	       (incf z (* c (column-l col))))
	      (t 
	       (setf (aref flags j) 'nonbasic-upper-bound)
	       (incf z (* c (column-u col)))))))
    ;; compute primal values
    (simplex-compute-primal-values sd)
    ;; finish compute objective value
    (dotimes (k m)
      (incf z (* (- (lp-obj-sense lp)) 
		 (column-c (adjvector-column-ref (lp-columns lp) (aref (basis-header b) k)))
		 (aref (basis-primal-values b) k))))
    (setf (basis-obj-value b) z)))



;;;;
(defun dual-simplex (sd &key (min-z) (max-z) (cutoff) (phase1only nil))
  (symbol-macrolet 
      ((b (simplex-basis sd))
       (z (basis-obj-value (simplex-basis sd))))
    (let ((n-iterations 0))
      ;; phase 1
      (loop
	 (print (float z))
	 (cond 
	   ((= z 0)
	    (return))
	   ((and cutoff (<= cutoff n-iterations))
	    (return-from dual-simplex 'iteration-cutoff)))
	 (incf n-iterations)
	 (let ((status (simplex-iteration sd)))
	   (cond ((eq status 'dual-feasible))
		 ((eq status 'optimal)
		  (if (< z 0)
		      (return-from dual-simplex 'unbounded)
		      (return)))
		 (t
		  (return-from dual-simplex status)))))
      ;; phase 2
      (check-dual-feasability sd)
      (when phase1only
	(return-from dual-simplex 'dual-feasible))
      (print 'phase2)
      (simplex-prepare-phase2 sd)
      (check-dual-feasability sd)
      (check-primal-values sd)
      (loop 
	 (print (float z))
	 (cond 
	   ((and min-z (<= z min-z))
	    (return-from dual-simplex 'reached-min-z))
	   ((and max-z (<= max-z z))
	    (return-from dual-simplex 'reached-max-z))
	   ((and cutoff (<= cutoff n-iterations))
	    (return-from dual-simplex 'iteration-cutoff)))
	 (incf n-iterations)
	 (let ((status (simplex-iteration sd)))
	   (check-primal-update-phase2 sd)
	   (check-primal-values sd)
	   (unless (eq status 'dual-feasible)
	     (return-from dual-simplex status)))))))



