(in-package :rationalsimplex)

;;;;; Dual simplex algorithm implementation
;;;;;



;;;; Returns T if basis needs to be refactorized from scratch
(defun simplex-basis-matrix-refactorize-p (sd)
  (zerop (mod (stats-total-iters (simplex-stats sd))
	      (basis-matrix-refactorization-period (basis-matrix (simplex-basis sd))))))
  


;;;; Refactorizes basis from scratch
(defun simplex-basis-matrix-refactorize (sd)
  (let* ((b (simplex-basis sd))
	 (bm (simplex-alt-basis-matrix sd))
	 (bh (basis-header b))
	 (leaving-ref (aref bh (simplex-pivot-row-index sd)))
	 (entering-ref (simplex-pivot-col-ref sd)))
    (fill-basis-matrix bm (simplex-lp sd) (basis-header b) leaving-ref entering-ref)
    (unless (basis-matrix-lu-factorization bm)
      (error "basis redundancy"))
    (when *checks*
      (let ((new-bh (copy-seq bh)))
	(setf (aref new-bh (simplex-pivot-row-index sd)) entering-ref)
	(check-u-seqs bm)
	(check-lu (simplex-lp sd) bm new-bh)))))


;;;; Computes spike for basis factorization update
(defun simplex-basis-matrix-update-compute-spike (sd)
  (let ((spike (tran-hsv (simplex-spike-ftran sd))))
    (set-column-as-simplex-spike-vector sd (simplex-pivot-col-ref sd))
    (ftran-l (simplex-spike-ftran sd) (simplex-spike-hsv sd))
    (hsv-remove-zeros spike)
    (hsv-normalize spike)))


;;;; Updates basis factorization
(defun simplex-basis-matrix-update (sd)
  (let* ((b (simplex-basis sd))
	 (bm (basis-matrix b))
	 (bh (basis-header b))
	 (spike (tran-hsv (simplex-spike-ftran sd))))
    (lu-update bm (simplex-pivot-row-index sd) spike)
    (when *checks*
      (let ((new-bh (copy-seq bh)))
	(setf (aref new-bh (simplex-pivot-row-index sd)) (simplex-pivot-col-ref sd))
	(check-u-seqs bm)
	(check-lu (simplex-lp sd) bm new-bh)))))


;;;; Final basis matrix update operation and thread waiting,
;;;; Ready for next iteration after this.
(defun simplex-finalize-basis-matrix-update (sd)
  (if (simplex-basis-matrix-refactorize-p sd)
      (let ((bm (simplex-alt-basis-matrix sd)))
	(rotatef (simplex-alt-basis-matrix sd) (basis-matrix (simplex-basis sd)))
	(setf (tran-bm (simplex-ftran sd)) bm
	      (tran-bm (simplex-flip-ftran sd)) bm
	      (tran-bm (simplex-dse-ftran sd)) bm
	      (tran-bm (simplex-spike-ftran sd)) bm
	      (tran-bm (simplex-btran sd)) bm))
      (simplex-basis-matrix-update sd)))


;;;; 
(defun simplex-iteration-dse-ftran (sd)
  (ftran (simplex-dse-ftran sd) (tran-hsv (simplex-btran sd)))
  (check-ftran (simplex-basis sd) (simplex-lp sd) 
	       (simplex-dse-ftran sd) (tran-hsv (simplex-btran sd))))


;;;;
(defun simplex-iteration-btran (sd)
  (btran (simplex-btran sd) (simplex-hsv sd))
  (check-btran (simplex-basis sd) (simplex-lp sd) 
	       (simplex-btran sd) (simplex-hsv sd)))


;;;;
(defun simplex-iteration-ftran (sd)
  (ftran (simplex-ftran sd) (simplex-hsv sd))
  (check-ftran (simplex-basis sd) (simplex-lp sd) 
	       (simplex-ftran sd) (simplex-hsv sd)))


;;;;; Performs an iteration of the dual simplex algorithm
(defun simplex-iteration (sd)
  ;; reset
  (reset-simplex sd)
  ;; pricing
  (choose-exiting-basis-index sd)
  (when (= -1 (simplex-pivot-row-index sd))
    (return-from simplex-iteration 'optimal))
  ;; btran
  (hsv-add (simplex-pivot-row-index sd) 1 (simplex-hsv sd))
  (simplex-iteration-btran sd)
  ;; dse ftran
  (thread-launch *dse-ftran-thread* 
		 #'(lambda () (simplex-iteration-dse-ftran sd)))
  ;; pivot row
  (compute-pivot-row sd)
  (check-pivot-row sd)
  ;; ratio test
  (when (= -1 (choose-entering-basis-index sd))
    (return-from simplex-iteration 'infeasible))
  ;;
  (when (simplex-basis-matrix-refactorize-p sd)
    ;; launch basis matrix refactorization 
    (thread-launch *basis-matrix-factor-thread* 
		   #'(lambda () (simplex-basis-matrix-refactorize sd))))
  ;; ftrans
  (set-column-as-simplex-vector sd (simplex-pivot-col-ref sd))
  (simplex-iteration-ftran sd)
  ;;
  (unless (simplex-basis-matrix-refactorize-p sd)
    ;; launch basis matrix factorization update
    (thread-launch *basis-matrix-factor-thread* 
		   #'(lambda () (simplex-basis-matrix-update-compute-spike sd))))
  ;; update DSE weights
  (thread-result *dse-ftran-thread*)
  (thread-launch *dse-weight-update-thread* 
		 #'(lambda () (simplex-basis-update-dse sd)))
  ;; basis change and update
  (simplex-basis-update sd)
  ;; complete the iteration
  (thread-result *basis-matrix-factor-thread*)
  (simplex-finalize-basis-matrix-update sd)
  (thread-result *dse-weight-update-thread*)
  (check-infeas-vector (simplex-basis sd) (simplex-lp sd))
  (check-reduced-costs sd)
  (check-dual-feasability sd)
  (check-dse-weights sd)
  'dual-feasible)



;;;; Computes reduced costs from scratch given a basis 
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
	     (c (column-c col))
	     (cd (denominator c)))
	(unless (zerop c)
	  (mulf cdenom (/ cd (gcd cd cdenom))))))
    (setf (hsv-coef (simplex-hsv sd)) (/ 1 cdenom)
	  (hsv-length (simplex-hsv sd)) 0)
    (dotimes (k m)
      (let* ((col (adjvector-column-ref (lp-columns lp) (aref bh k)))
	     (c (column-c col)))
	(assert (= (denominator c) (gcd (denominator c) cdenom)))
	(assert (integerp  (* (numerator c) (/ cdenom (denominator c)))))
	(unless (zerop c)
	  (hsv-add k (* (numerator c) (/ cdenom (denominator c))) (simplex-hsv sd)))))
    ;; compute multipliers
    (btran (simplex-btran sd) (simplex-hsv sd))
    (check-btran b lp (simplex-btran sd) (simplex-hsv sd))
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




;;;; Computes primal values from scratch given a basis
(defun simplex-compute-primal-values (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b))
	 (tr (simplex-ftran sd))
	 (m (basis-matrix-size bm))
	 (flags (basis-column-flags b))
	 (infeas (basis-primal-infeas b))
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
    (check-ftran b lp tr (simplex-hsv sd))
    ;; fill primal value array and infeasability vector
    (dotimes (k (hsv-length (tran-hsv tr)))
      (let ((i (aref (hsv-is (tran-hsv tr)) k))
	    (x (* (hsv-coef (tran-hsv tr)) (aref (hsv-vis (tran-hsv tr)) k))))
	(setf (aref (basis-primal-values b) i) x)))
    (reset-splay-tree-fixnum-rational infeas)
    (dotimes (i m)
      (let ((col (adjvector-column-ref (lp-columns lp) (aref (basis-header b) i)))
	    (x (aref (basis-primal-values b) i)))
	(cond ((and (column-has-l col) (< x (column-l col)))
	       (splay-tree-fixnum-rational-set infeas i 
					       (* (- (column-l col) x)
						  (- (column-l col) x))))
	      ((and (column-has-u col) (< (column-u col) x))
	       (splay-tree-fixnum-rational-set infeas i 
					       (* (- x (column-u col))
						  (- x (column-u col))))))))))
	



;;;; Given a dual-feasible basis, prepares phase 2 of the simplex algorithm
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



;;;; Executes the dual simplex algorithm
;;;; allows various interruptions as parameters
(defun dual-simplex (sd &key  
		     (min-z) 
		     (max-z)
		     (z-print-freq 1)
		     (max-total-time)
		     (max-phase-time)
		     (max-total-iters) 
		     (max-phase-iters))
  (let ((start-time (get-internal-real-time))
	(phase2-start-time 0))
    (symbol-macrolet 
	((b (simplex-basis sd))
	 (st (simplex-stats sd))
	 (z (basis-obj-value (simplex-basis sd))))
      (macrolet 
	  ((print-z ()
	     `(format t "~&    ~12,D        ~16,5F~%"
		      (stats-total-iters st)
		      (coerce z 'double-float)))
	   (exit (status) 
	     `(progn 
		(setf (stats-total-duration st) 
		      (/ (float (- (get-internal-real-time) start-time))
			 (float internal-time-units-per-second)))
		(if (basis-in-phase1 b) 
		    (setf (stats-phase1-duration st) (stats-total-duration st))
		    (setf (stats-phase2-duration st) 
			  (/ (float (- (get-internal-real-time) phase2-start-time))
			     (float internal-time-units-per-second))))
		(return-from dual-simplex (prog1 ,status (print-z))))))
	(check-infeas-vector b (simplex-lp sd))
	;; phase 1
	(format t "~&      Iterations       Phase 1 objective")
	(print-z)
	(loop
	   (cond 
	     ((= z 0)
	      (return))
	     ((and max-phase-iters (<= max-phase-iters (stats-phase1-iters st)))
	      (exit 'phase1-iteration-count-cutoff))
	     ((and max-total-iters (<= max-total-iters (stats-total-iters st)))
	      (exit 'total-iteration-count-cutoff))
	     ((and max-phase-time 
		   (< (* max-phase-time internal-time-units-per-second)
		      (- (get-internal-real-time) start-time)))
	      (exit 'phase1-duration-cutoff))
	     ((and max-total-time 
		   (< (* max-total-time internal-time-units-per-second)
		      (- (get-internal-real-time) start-time)))
	      (exit 'total-duration-cutoff)))
	   (incf (stats-phase1-iters st))
	   (incf (stats-total-iters st))
	   (when (zerop (mod (stats-total-iters st) z-print-freq))
	     (print-z))
	   (let ((status (simplex-iteration sd)))
	     (cond ((eq status 'dual-feasible))
		   ((eq status 'optimal)
		    (if (< z 0)
			(exit 'unbounded)
			(return)))
		   (t
		    (exit status)))))
	;; phase 2
	(setf (stats-phase1-duration st) 
	      (/ (float (- (get-internal-real-time) start-time))
		 (float internal-time-units-per-second)))
	(print-z)
	(check-dual-feasability sd)
	(format t "~&Dual-feasible basis found, going to phase 2.")
	(format t "~&      Iterations       Phase 2 objective (~A)" 
		(lp-obj-name (simplex-lp sd)))
	(setf phase2-start-time (get-internal-real-time))
	(simplex-prepare-phase2 sd)
	(print-z)	
	(check-infeas-vector b (simplex-lp sd))
	(check-dual-feasability sd)
	(check-primal-values sd)
	(loop 
	   (cond 
	     ((and min-z (<= z min-z))
	      (exit 'reached-min-z))
	     ((and max-z (<= max-z z))
	      (exit 'reached-max-z))
	     ((and max-phase-iters (<= max-phase-iters (stats-phase2-iters st)))
	      (exit 'phase2-iteration-count-cutoff))
	     ((and max-total-iters (<= max-total-iters (stats-total-iters st)))
	      (exit 'total-iteration-count-cutoff))
	     ((and max-phase-time 
		   (< (* max-phase-time internal-time-units-per-second)
		      (- (get-internal-real-time) phase2-start-time)))
	      (exit 'phase2-duration-cutoff))
	     ((and max-total-time 
		   (< (* max-total-time internal-time-units-per-second)
		      (- (get-internal-real-time) start-time)))
	      (exit 'total-duration-cutoff)))
	   (incf (stats-phase2-iters (simplex-stats sd)))
	   (incf (stats-total-iters (simplex-stats sd)))
	   (when (zerop (mod (stats-total-iters st) z-print-freq))
	     (print-z))
	   (let ((status (simplex-iteration sd)))
	     (check-primal-update-phase2 sd)
	     (check-primal-values sd)
	     (unless (eq status 'dual-feasible)
	       (exit status))))))))
		 




;;;;; DEBUGGING

;;;;
(defun check-primal-update-phase2 (sd)
  (when *checks*
    (let ((orig-primal (copy-seq (basis-primal-values (simplex-basis sd)))))
      ;; get multipliers
      (simplex-compute-primal-values sd)
      (let ((test 
	     (dotimes (k (length orig-primal) t)
	       (unless (= (aref orig-primal k)
			  (aref (basis-primal-values (simplex-basis sd)) k))
		 (print (cons k (aref (basis-header (simplex-basis sd)) k)))
		 (return nil)))))
	(unless test
	  (error "phase2 primal update"))))))
  
  

;;;;
(defun check-reduced-costs (sd)
  (when *checks*
    (let ((orig-rcosts (copy-seq (basis-reduced-costs (simplex-basis sd)))))
      ;; get multipliers
      (simplex-compute-reduced-costs sd)
      (let ((test 
	     (dotimes (j (length orig-rcosts) t)
	       (let ((flag (aref (basis-column-flags (simplex-basis sd)) j)))
		 (when (or (eq flag 'nonbasic-lower-bound)
			   (eq flag 'nonbasic-upper-bound))
		   (unless (= (aref orig-rcosts j)
			      (aref (basis-reduced-costs (simplex-basis sd)) j))
		     (return nil)))))))
	(unless test
	  (let ((diff (copy-seq orig-rcosts)))
	     (dotimes (j (length orig-rcosts))
	       (let ((flag (aref (basis-column-flags (simplex-basis sd)) j)))
		 (if (or (eq flag 'nonbasic-lower-bound)
			 (eq flag 'nonbasic-upper-bound))
		     (decf (aref diff j) (aref (basis-reduced-costs (simplex-basis sd)) j))
		     (setf (aref diff j) 0))))
	     (print '--)
	     (print orig-rcosts)
	     (print '--should-be)
	     (print (basis-reduced-costs (simplex-basis sd)))
	     (print '---diff)
	     (print diff)
	  (error "reduced costs")))))))



  


