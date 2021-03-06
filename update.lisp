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



;;;; Updates dual-steepest-edge weights
(defun simplex-basis-update-dse (sd)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (let* ((b (simplex-basis sd))
	 (dse-tr (simplex-dse-ftran sd))
	 (alphaq (tran-hsv (simplex-ftran sd)))
	 (exit-row (simplex-pivot-row-index sd))
	 (alpha-index (hsv-find exit-row alphaq))
	 (alphaqr (aref (hsv-vis alphaq) alpha-index))
	 (m (length (basis-header b)))
	 (tau (tran-hsv dse-tr)) 
	 (flags (simplex-dse-update-flags sd)))
    ;; reset flags
    (bit-xor flags flags t)
    ;; build new dse-coef and new dse-vis
    (let* ((alphaq-n (numerator (hsv-coef alphaq)))
	   (alphaq-d (denominator (hsv-coef alphaq)))
	   (tau-n (numerator (hsv-coef tau)))
	   (tau-d (denominator (hsv-coef tau)))
	   (dse-n (numerator (basis-dse-coef b)))
	   (dse-d (denominator (basis-dse-coef b)))
	   (a (* tau-d dse-n))
	   (d (* alphaqr alphaqr))
	   (c (* alphaq-n alphaq-n))
	   (f (aref (basis-dse-weight-vis b) exit-row))
	   (g (* 2 alphaqr tau-n dse-d))
	   (h (gcd d f))
	   (a.h (* a h))
	   (j (gcd a.h g))
	   (c.j (* c j))
	   (l (* alphaq-d alphaq-d a f))
	   (gcd-c.j-l (gcd c.j l))
	   (c/ (/ c.j gcd-c.j-l))
	   (a/ (/ a.h j))
	   (g/ (/ g j))
	   (d/ (/ d h))
	   (f/ (/ f h)))
      (declare (integer alphaq-n alphaq-d tau-n tau-d dse-n dse-d))
      (declare (integer c/ a/ g/ d/ f/))
      ;; update coefficient
      (setf (basis-dse-coef b) (/ gcd-c.j-l (* dse-d d c tau-d)))
      ;; update dse weights for nonzero elements in alphaq vector
      (dotimes (k (hsv-length alphaq))
	(let ((i (aref (hsv-is alphaq) k))
	      (alphaqi (aref (hsv-vis alphaq) k)))
	  (unless (zerop alphaqi)
	    (setf (sbit flags i) 1)
	    (setf (aref (basis-dse-weight-vis b) i)
		  (if (= i exit-row)
		      (/ l gcd-c.j-l)
		      (let ((taui-ref (hsv-find i tau)))
			(declare (fixnum taui-ref))
			(if (= -1 taui-ref)
			    (* c/ a/ (+ (* d/ (aref (basis-dse-weight-vis b) i)) 
					(* f/ alphaqi alphaqi)))
			    (* c/ (- (* a/ (+ (* d/ (aref (basis-dse-weight-vis b) i)) 
					      (* f/ alphaqi alphaqi)))
				     (* g/ alphaqi (aref (hsv-vis tau) taui-ref)))))))))))
      ;; update all other dse weights
      (dotimes (i m)
	(when (zerop (sbit flags i))
	  (setf (aref (basis-dse-weight-vis b) i) 
		(* c/ a/ d/ (aref (basis-dse-weight-vis b) i))))))
    ;; normalize signs and values
    (let ((coef-signum (signum (basis-dse-coef b)))
	  (common-factor (aref (basis-dse-weight-vis b) 0)))
      (declare (integer common-factor))
      (dotimes (i m)
	(when (= 1 (setf common-factor
			 (gcd common-factor (aref (basis-dse-weight-vis b) i))))
	  (return)))
      (when (> 0 coef-signum)
	(setf common-factor (- common-factor)))
      (unless (= 1 common-factor)
	(setf (basis-dse-coef b) (* (basis-dse-coef b) common-factor))
	(dotimes (i m)
	  (let ((r (/ (aref (basis-dse-weight-vis b) i) common-factor)))
	    (declare (integer r))
	    (setf (aref (basis-dse-weight-vis b) i) r)))))))







;;;; Updates primal values, primal infeasabilities, and dual-steepest-edge weights
(defun simplex-basis-update-primal (sd)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (let* ((b (simplex-basis sd))
	 (alphaq (tran-hsv (simplex-ftran sd)))
	 (exit-row (simplex-pivot-row-index sd))
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
    ;; update values
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
	 (basis-in-phase1 b))))))



;;;; Temporarily uses pivot row arrays to compute weighted column sum
(defun simplex-compute-flip-ftran (sd)
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
	(check-flip-rhs sd))
      ;; perform ftran
      (ftran (simplex-flip-ftran sd) (simplex-hsv sd))
      (check-ftran (simplex-basis sd) lp (simplex-flip-ftran sd) (simplex-hsv sd)))))
      



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
      (simplex-compute-flip-ftran sd)
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
	     (simplex-delta sd)))))

