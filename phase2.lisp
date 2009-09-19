
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
      (let* ((col (aref (lp-columns lp) (aref bh k)))
	     (cd (denominator (column-c col))))
	(mulf cdenom (/ cd (gcd cd cdenom)))))
    (setf (simplex-vector-coef sd) (/ 1 cdenom)
	  (fill-pointer (simplex-vector-indices sd)) 0
	  (fill-pointer (simplex-vector-values sd)) 0)
    (dotimes (k m)
      (let* ((col (aref (lp-columns lp) (aref bh k)))
	     (c (column-c col)))
	(assert (= (denominator c) (gcd (denominator c) cdenom)))
	(vector-push-extend k (simplex-vector-indices sd))
	(vector-push-extend (* (numerator c) (/ cdenom (denominator c)))
			    (simplex-vector-values sd))))
    ;; compute multipliers
    (simplex-btran sd)
    (copy-vector-from-tran sd)
    ;; compute reduced costs
    (assert (= n (length (simplex-pivot-row sd))))
    (dotimes (j n)
      (setf (aref (simplex-pivot-row sd) j)
	    (- (column-c (aref (lp-columns lp) j)))))
    (compute-pivot-row sd)
    (dotimes (j n)
      (let ((flag (aref flags j)))
	(setf (aref rcosts j) 
	      (if (eq flag 'basic)
		  0
		  (- (aref (simplex-pivot-row sd) j))))
	(setf (aref (simplex-pivot-row sd) j) 0)))))



;;;;
(defun simplex-compute-primal-values (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b))
	 (tr (simplex-tran sd))
	 (m (basis-matrix-size bm))
	 (flags (basis-column-flags b))
	 (rhs (make-nvector m 0 rational))
	 (cdenom 1))
    ;; compute right-hand side in solving Bb = rhs for b
    (dotimes (i m)
      (let* ((v 0)
	     (row-ref (aref (lp-active-row-refs lp) i))
	     (row (aref (lp-rows lp) row-ref)))
	(dotimes (k (length (row-col-refs row)))
	  (let* ((col-ref (aref (row-col-refs row) k))
		 (flag (aref flags col-ref))
		 (col (aref (lp-columns lp) col-ref))
		 (a (rational-in-column col (aref (row-col-indices row) k))))
	    (cond ((eq flag 'nonbasic-lower-bound)
		   (assert (column-has-l col))
		   (decf v (* a (column-l col))))
		  ((eq flag 'nonbasic-upper-bound)
		   (assert (column-has-u col))
		   (decf v (* a (column-u col)))))))
	(setf (aref rhs i) v)
	(mulf cdenom (/ (denominator v) (gcd (denominator v) cdenom)))))
    ;; set rhs for solving Bb = rhs for b
    (setf (simplex-vector-coef sd) (/ 1 cdenom)
	  (fill-pointer (simplex-vector-indices sd)) 0
	  (fill-pointer (simplex-vector-values sd)) 0)
    (dotimes (i m)
      (setf (aref (basis-primal-values b) i) 0)
      (let ((v (aref rhs i)))
	(unless (zerop v)
	  (vector-push-extend i (simplex-vector-indices sd))
	  (vector-push-extend (* (numerator v) (/ cdenom (denominator v)))
			      (simplex-vector-values sd)))))
    (simplex-ftran sd)
    ;; fill primal value array
    (dotimes (k (length (tran-indices tr)))
      (let ((i (aref (tran-indices tr) k))
	    (x (* (tran-coef tr) (aref (tran-values tr) k))))
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
    (dotimes (j (length (lp-columns lp)))
      (let* ((col (aref (lp-columns lp) j))
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
		 (column-c (aref (lp-columns lp) (aref (basis-header b) k)))
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



