(in-package :rationalsimplex)

;;;;; Simplex object and basis updating functions
;;;;; at the end of a dual simplex iteration
;;;;;



;;;; Determines non-basic primal value flag for exiting variable
(defun basis-exiting-variable-flag (in-phase1 exit-val exit-rcost exit-col)
  (if in-phase1
      (cond ((and (not (column-has-l exit-col))
		  (not (column-has-u exit-col)))
	     (error "exiting variable is free in phase 2"))
	    ((and (= -1 exit-val)
		  (not (column-has-l exit-col)))
	     'nonbasic-lower-bound)
	    ((and (= 0 exit-val)
		  (column-has-l exit-col)
		  (<= 0 exit-rcost))
	     'nonbasic-lower-bound)
	    ((and (= 0 exit-val)
		  (column-has-u exit-col)
		  (>= 0 exit-rcost))
	     'nonbasic-upper-bound)
	    ((and (= 1 exit-val)
		  (not (column-has-u exit-col)))
	     'nonbasic-upper-bound)
	    (t
	     (error "bad exiting variable or reduced cost value in phase 1")))
      (cond ((and (column-has-l exit-col)
		  (= (column-l exit-col) exit-val)
		  (<= 0 exit-rcost))
	     'nonbasic-lower-bound)
	    ((and (column-has-u exit-col)
		  (= (column-u exit-col) exit-val)
		  (>= 0 exit-rcost))
	     'nonbasic-upper-bound)
	    ((and (not (column-has-l exit-col))
		  (not (column-has-u exit-col)))
	     (error "exiting variable is free in phase 2"))
	    (t
	     (error "bad exiting variable value or reduced cost in phase 2")))))



;;;; Updates primal values, primal infeasabilities, and dual-steepest-edge weights
(defun simplex-basis-update-primal (sd)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (let* ((b (simplex-basis sd))
	 (dse-tr (simplex-dse-ftran sd))
	 (alphaq (tran-hsv (simplex-ftran sd)))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-row-weight (aref (basis-dse-weights b) exit-row))
	 (primal-values (basis-primal-values b))
	 (lpcols (lp-columns (simplex-lp sd)))
	 (alpha-index (hsv-find exit-row alphaq)))
    (declare (fixnum alpha-index))
    (assert (/= -1 alpha-index))
    (assert (/= 0 (aref (hsv-vis alphaq) alpha-index)))
    ;; compute primal step
    (setf (simplex-primal-step sd)
	  (/ (simplex-delta sd)
	     (* (aref (hsv-vis alphaq) alpha-index)
		(hsv-coef alphaq))))
    ;; update values and weights
    (dotimes (k (hsv-length alphaq))
      (let ((i (aref (hsv-is alphaq) k)))
	;; update basic primal values
	(decf (aref primal-values i) 
	      (* (simplex-primal-step sd)
		 (hsv-coef alphaq)
		 (aref (hsv-vis alphaq) k)))
	;; update primal infeasabilities
	(basis-update-primal-infeasability 
	 (basis-primal-infeas b)
	 i
	 (aref primal-values i) 
	 (adjvector-column-ref lpcols (aref (basis-header b) i))
	 (basis-in-phase1 b))
	;; update dse weights
	(if (= k alpha-index)
	    (let ((d (* (hsv-coef alphaq)
			(aref (hsv-vis alphaq) alpha-index))))
	      (assert (= exit-row-weight (aref (basis-dse-weights b) exit-row)))
	      (divf (aref (basis-dse-weights b) exit-row) (* d d)))
	    (let ((ratio (/ (aref (hsv-vis alphaq) k)
			    (aref (hsv-vis alphaq) alpha-index)))
		  (tau-index (hsv-find i (tran-hsv dse-tr))))
	      (declare (fixnum tau-index))
	      (unless (= -1 tau-index)
		(decf (aref (basis-dse-weights b) i)
		      (* 2 ratio 
			 (hsv-coef (tran-hsv dse-tr))
			 (aref (hsv-vis (tran-hsv dse-tr)) tau-index))))
	      (incf (aref (basis-dse-weights b) i)
		    (* ratio ratio exit-row-weight))))))))



;;;; Updates the basis matrix
;;;; Either rebuilding or just updating the LU factorization
(defun simplex-basis-matrix-update (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b)))
    (if (zerop (mod (stats-total-iters (simplex-stats sd))
		    (basis-matrix-refactorization-period bm)))
	(progn 
	  (fill-basis-matrix bm lp (basis-header b))
	  (unless (basis-matrix-lu-factorization bm)
	    (error "basis redundancy")))
	(let ((spike (tran-hsv (simplex-ftran sd))))
	  (set-column-as-simplex-vector sd (simplex-pivot-col-ref sd))
	  (ftran-l (simplex-ftran sd) (simplex-hsv sd))
	  (hsv-remove-zeros spike)
	  (hsv-normalize spike)
	  (lu-update bm (simplex-pivot-row-index sd) spike)))
    (check-u-seqs bm)
    (check-lu lp bm (basis-header b))))




;;;; Temporarily uses pivot row arrays to compute weighted column sum
(defun simplex-compute-flip-ftran-rhs (sd)
  (let* ((lp (simplex-lp sd))
	 (inphase1 (basis-in-phase1 (simplex-basis sd)))
	 (n-cols (simplex-n-flips sd))
	 (last-heap-index (- n-cols 1))
	 (heap (simplex-add-heap sd))
	 (rhs-is (simplex-pivot-row-col-refs sd))
	 (rhs-values (simplex-pivot-row-values sd))
	 (rhs-length -1)
	 (flags (basis-column-flags (simplex-basis sd)))
	 (prev-row-ref -1))
    (labels ((get-row-ref (k)
	       (let* ((counter (aref (simplex-add-counters sd) k))
		      (col-ref (aref (simplex-flip-col-refs sd) k))
		      (col-hsv (column-hsv (adjvector-column-ref (lp-columns lp) col-ref))))
		 (if (= counter (hsv-length col-hsv))
		     most-positive-fixnum
		     (aref (hsv-is col-hsv) counter))))
	     (sift-down (root end)
	       (loop
		  (let ((child (+ 1 (* 2 root))))
		    (when (> child end)
		      (return))
		    (when (and (< child end)
			       (> (get-row-ref (aref heap child))
				  (get-row-ref (aref heap (+ child 1)))))
		      (incf child))
		    (when (<= (get-row-ref (aref heap root)) 
			      (get-row-ref (aref heap child)))
		      (return))
		    (rotatef (aref heap root) (aref heap child))
		    (setf root child)))))
      ;; initialization
      (assert (not (zerop n-cols)))
      (reset-hsv (simplex-hsv sd))
      (dotimes (k n-cols)
	(let* ((col-ref (aref (simplex-flip-col-refs sd) k))
	       (col (adjvector-column-ref (lp-columns lp) col-ref))
	       (diff (if (eq (aref flags col-ref) 'nonbasic-lower-bound)
			 (column-u-minus-l inphase1 col)
			 (column-l-minus-u inphase1 col))))
	  (setf (aref (simplex-add-counters sd) k) 0
		(aref (simplex-flip-col-coefs sd) k) (* diff (hsv-coef (column-hsv col)))
		(aref (simplex-add-heap sd) k) k)))
      ;; heapify
      (loop for heapify-start from (floor (- last-heap-index 1) 2) downto 0
	 do (sift-down heapify-start last-heap-index))
      ;; build weighted column sum
      (loop
	 ;; pop element with smallest row index off heap
	 (let* ((smallest-k (aref heap 0))
		(row-ref (get-row-ref smallest-k)))
	   (cond ((= row-ref most-positive-fixnum)
		  (incf rhs-length)
		  (return))
		 ((> prev-row-ref row-ref)
		  (error "bad row index order"))
		 ((= prev-row-ref row-ref))
		 (t
		  (incf rhs-length)
		  (setf prev-row-ref row-ref
			(aref rhs-is rhs-length) row-ref
			(aref rhs-values rhs-length) 0)))
	   ;; update rhs value
	   (incf (aref rhs-values rhs-length)
		 (* (aref (simplex-flip-col-coefs sd) smallest-k)
		    (aref (hsv-vis (column-hsv (adjvector-column-ref 
						(lp-columns lp)
						(aref (simplex-flip-col-refs sd) smallest-k))))	
			  (aref (simplex-add-counters sd) smallest-k))))
	   ;; update counter and heap
	   (incf (aref (simplex-add-counters sd) smallest-k))
	   (sift-down 0 last-heap-index)))
      ;; normalize rhs
      (let ((cdenom (common-denominator rhs-values rhs-length)))
	(setf (hsv-coef (simplex-hsv sd)) (/ 1 cdenom))
	(dotimes (k rhs-length)
	  (assert (integerp (* cdenom (aref rhs-values k))))
	  (hsv-add (adjvector-fixnum-ref (lp-active-row-inds lp) (aref rhs-is k))
		   (* cdenom (aref rhs-values k)) (simplex-hsv sd)))
	(hsv-sort-indices-increasing (simplex-hsv sd))
	(hsv-normalize (simplex-hsv sd))
	(check-flip-rhs sd)))))



;;;; Ties it all together.
;;;; Updates reduced costs, computes primal step, updates basis,
;;;; as well as the entering and exiting variable primal values.
;;;; Performs bound flips, updates objective, primal values and infeasabilites.
(defun simplex-basis-update (sd)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (let* ((b (simplex-basis sd))
	 (header (basis-header b))
	 (lp (simplex-lp sd))
	 (q (simplex-pivot-col-ref sd))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-col-ref (aref header exit-row))
	 (primal-values (basis-primal-values b))
	 (primal-infeas (basis-primal-infeas b))
	 (flags (basis-column-flags b))
	 (rcosts (basis-reduced-costs b)))
    ;; update reduced costs
    (setf (aref rcosts exit-col-ref) (- (simplex-dual-step sd)))
    (unless (zerop (simplex-dual-step sd))
      (dotimes (k (simplex-pivot-row-length sd))
	(let* ((j (aref (simplex-pivot-row-col-refs sd) k))
	       (flag (aref flags j)))
	  (when (or (eq flag 'nonbasic-lower-bound)
		    (eq flag 'nonbasic-upper-bound))
	    (decf (aref rcosts j)
		  (* (simplex-dual-step sd) 
		     (aref (simplex-pivot-row-values sd) k)))))))
    ;; in case of bound flips
    (unless (zerop (simplex-n-flips sd))
      (simplex-compute-flip-ftran-rhs sd)
      (ftran (simplex-flip-ftran sd) (simplex-hsv sd))
      (check-ftran b lp (simplex-flip-ftran sd) (simplex-hsv sd))
      ;; update objective, primal basic values and primal infeasabilities
      (let ((delta-x (tran-hsv (simplex-flip-ftran sd))))
	(dotimes (k (hsv-length delta-x))
	  (let* ((i (aref (hsv-is delta-x) k))
		 (xk (* (hsv-coef delta-x) (aref (hsv-vis delta-x) k)))
		 (col (adjvector-column-ref (lp-columns lp) (aref (basis-header b) i))))
	    (decf (aref primal-values i) xk)
	    (decf (basis-obj-value b) (* xk (column-c col)))
	    (basis-update-primal-infeasability primal-infeas i
					       (aref primal-values i) col 
					       (basis-in-phase1 b)))))
      ;; update primal nonbasic values
      (dotimes (k (simplex-n-flips sd))
	(let ((j (aref (simplex-flip-col-refs sd) k)))
	    (setf (aref flags j)
		  (if (eq (aref flags j) 'nonbasic-lower-bound)
		      'nonbasic-upper-bound
		      'nonbasic-lower-bound)))))
    ;; compute primal step, update primal basic values, update dse weights
    (simplex-basis-update-primal sd)
    ;; update exiting variable
    (setf (aref flags exit-col-ref)
	  (basis-exiting-variable-flag (basis-in-phase1 b) 
				       (aref primal-values exit-row)
				       (aref rcosts exit-col-ref)
				       (adjvector-column-ref (lp-columns lp) exit-col-ref)))
    ;; update entering variable
    (setf (aref header exit-row) q)
    (setf (aref primal-values exit-row) 
	  (+ (simplex-primal-step sd)
	     (nonbasic-value (basis-in-phase1 b)
			     (adjvector-column-ref (lp-columns lp) q)
			     (aref flags q))))
    (basis-update-primal-infeasability primal-infeas 
				       exit-row
				       (aref primal-values exit-row) 
				       (adjvector-column-ref (lp-columns lp) q)
				       (basis-in-phase1 b))
    (setf (aref flags q) 'basic)
    ;; update objective
    (incf (basis-obj-value b)
	  (* (simplex-dual-step sd)
	     (simplex-delta sd)))
    ;; update basis matrix
    (simplex-basis-matrix-update sd)))

