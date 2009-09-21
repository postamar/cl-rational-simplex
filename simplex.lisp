(in-package :rationalsimplex)

;;;;; Simplex data structure implementation
;;;;; An object contains all data necessary during a dual simplex iteration
;;;;;



;;;; Statistics data structure
(defstruct stats
  (total-duration 0.0 :type float)
  (phase1-duration 0.0 :type float)
  (phase2-duration 0.0 :type float)
  (total-iters 0 :type fixnum)
  (phase1-iters 0 :type fixnum)
  (phase2-iters 0 :type fixnum))



;;;; Simplex data structure
(symbol-macrolet 
    ((err (error "simplex constructor")))
  (defstruct (simplex
	       (:constructor %make-simplex)
	       (:print-function print-simplex))
    (lp                  err :type lp)
    (stats               err :type stats)
    (basis               err :type basis)
    (alt-basis-matrix    err :type basis-matrix)
    (btran               err :type tran)
    (ftran               err :type tran)
    (flip-ftran          err :type tran)
    (dse-ftran           err :type tran)
    (spike-ftran         err :type tran)
    (hsv                 err :type hsv)
    (spike-hsv           err :type hsv)
    (pivot-row-col-refs  err :type (simple-array fixnum 1))
    (pivot-row-values    err :type (simple-array rational 1))
    (pivot-row-length    0   :type fixnum)
    (breakpoints         err :type simple-bit-vector)
    (n-breakpoints       0   :type fixnum)
    (pivot-row-index     -1  :type fixnum)
    (pivot-col-ref       -1  :type fixnum)
    (delta               0   :type rational)
    (dual-step           0   :type rational)
    (primal-step         0   :type rational)
    (add-heap            err :type (simple-array fixnum 1))
    (add-counters        err :type (simple-array fixnum 1))
    (add-refs            err :type (simple-array fixnum 1))
    (flip-col-refs       err :type (simple-array fixnum 1))
    (flip-col-coefs      err :type (simple-array rational 1))
    (n-flips             0   :type fixnum)
    (dse-update-flags    err :type simple-bit-vector)))



;;;; Printer function
(defun print-simplex (sd stream depth)
  (declare (ignore depth))
  (format stream "#SIMPLEX:~% N-BREAKPOINTS  ~A~%  PIVOT ROW INDEX  ~A~%  PIVOT COL REF  ~A~%  EXIT VAR DELTA  ~A~%  DUAL STEP  ~A~%  PRIMAL STEP  ~A~%PIVOT ROW   ~A~%RHS ~A~%"
	  (simplex-n-breakpoints sd)
	  (simplex-pivot-row-index sd)
	  (simplex-pivot-col-ref sd)
	  (simplex-delta sd)
	  (simplex-dual-step sd)
	  (simplex-primal-step sd)
	  (let ((v (make-array (length (simplex-pivot-row-col-refs sd))
			       :element-type 'rational)))
	    (loop for k from 0 below (simplex-pivot-row-length sd)
	       do (setf (aref v (aref (simplex-pivot-row-col-refs sd) k))
			(aref (simplex-pivot-row-values sd) k)))
	    v)
	  (simplex-hsv sd)))



;;;; Constructor
(defun make-simplex (lp basis)
  (let ((n (adjvector-column-fill-pointer (lp-columns lp)))
	(m (length (basis-header basis)))
	(refac-period (basis-matrix-refactorization-period (basis-matrix basis))))
    (%make-simplex
     :lp                 lp
     :stats              (make-stats)
     :basis              basis
     :alt-basis-matrix   (make-basis-matrix :lp lp :refactorization-period refac-period)
     :btran              (make-tran (basis-matrix basis))
     :ftran              (make-tran (basis-matrix basis))
     :flip-ftran         (make-tran (basis-matrix basis))
     :dse-ftran          (make-tran (basis-matrix basis))
     :spike-ftran        (make-tran (basis-matrix basis))
     :hsv                (make-hsv)
     :spike-hsv          (make-hsv)
     :pivot-row-col-refs (make-array n :initial-element -1 :element-type 'fixnum)
     :pivot-row-values   (make-array n :initial-element 0 :element-type 'rational)
     :breakpoints        (make-array n :initial-element 0 :element-type 'bit)
     :add-heap           (make-array n :initial-element -1 :element-type 'fixnum)
     :add-counters       (make-array n :initial-element -1 :element-type 'fixnum)
     :add-refs           (make-array n :initial-element -1 :element-type 'fixnum)
     :flip-col-refs      (make-array n :initial-element -1 :element-type 'fixnum)
     :flip-col-coefs     (make-array n :initial-element 0 :element-type 'rational)
     :dse-update-flags   (make-array m :element-type 'bit))))

    

;;;; Resets the data structure to start new iteration
(defun reset-simplex (sd)
  (reset-hsv (simplex-hsv sd))
  (bit-xor (simplex-breakpoints sd) (simplex-breakpoints sd) t)
  (setf (simplex-n-breakpoints sd) 0
	(simplex-pivot-row-index sd) -1
	(simplex-pivot-col-ref sd) -1
	(simplex-delta sd) 0
	(simplex-dual-step sd) 0
	(simplex-primal-step sd) 0
	(simplex-n-flips sd) 0))
 


;;;; Sets a column of the LP (usually the entering column)
;;;; in the simplex object's intermediate result vector
(defun set-column-as-simplex-vector (sd q)
  (let* ((lp (simplex-lp sd))
	 (col (adjvector-column-ref (lp-columns lp) q)))
    (reset-hsv (simplex-hsv sd))
    (setf (hsv-coef (simplex-hsv sd)) (hsv-coef (column-hsv col)))
    (dotimes (k (hsv-length (column-hsv col)))
      (let ((ind (adjvector-fixnum-ref (lp-active-row-inds lp)
				       (aref (hsv-is (column-hsv col)) k))))
	(unless (= -1 ind)
	  (hsv-add ind (aref (hsv-vis (column-hsv col)) k) (simplex-hsv sd)))))))


;;;; Sets a column of the LP (usually the entering column)
;;;; in the simplex object's intermediate result vector for spikes
(defun set-column-as-simplex-spike-vector (sd q)
  (let* ((lp (simplex-lp sd))
	 (col (adjvector-column-ref (lp-columns lp) q)))
    (reset-hsv (simplex-spike-hsv sd))
    (setf (hsv-coef (simplex-spike-hsv sd)) (hsv-coef (column-hsv col)))
    (dotimes (k (hsv-length (column-hsv col)))
      (let ((ind (adjvector-fixnum-ref (lp-active-row-inds lp)
				       (aref (hsv-is (column-hsv col)) k))))
	(unless (= -1 ind)
	  (hsv-add ind (aref (hsv-vis (column-hsv col)) k) (simplex-spike-hsv sd)))))))


;;;;; DEBUGGING

;;;;
(defun check-dual-feasability (sd)
  (when *checks*
    (let* ((b (simplex-basis sd))
	   (rcosts (basis-reduced-costs b))
	   (flags (basis-column-flags b)))
      (assert (= (length rcosts) (length flags)))
      (let ((test 
	     (dotimes (j (length flags) t)
	       (let ((flag (aref flags j))
		     (d (aref rcosts j)))
		 (when (eq flag 'nonbasic-lower-bound)
		   (unless (<= 0 d)
		     (print flag)
		     (print d)
		     (return nil)))
		 (when (eq flag 'nonbasic-upper-bound)
		   (unless (<= d 0)
		     (print flag)
		     (print d)
		     (return nil)))))))
	(unless test
	  (error "current basis not dual feasible"))))))
  
  
;;;;
(defun check-primal-feasability (sd)
  (when *checks*
    (let* ((b (simplex-basis sd))
	   (bh (basis-header b))
	   (lp (simplex-lp sd)))
      (check-primal-values sd)
      (dotimes (k (length bh))
	(let* ((col (adjvector-column-ref (lp-columns lp) (aref bh k)))
	       (xk (aref (basis-primal-values b) k)))
	  (when (column-has-l col)
	    (assert (<= (column-l col) xk)))
	  (when (column-has-u col)
	    (assert (<= xk (column-u col)))))))))
  
  

;;;;
(defun check-primal-values (sd)
  (when *checks*
    (let* ((b (simplex-basis sd))
	   (bh (basis-header b))
	   (lp (simplex-lp sd))
	   (n (adjvector-column-fill-pointer (lp-columns lp)))
	   (x (make-array n :initial-element 0 :element-type 'rational)))
      ;; fill x
      (dotimes (k (length bh))
	(let* ((col (adjvector-column-ref (lp-columns lp) (aref bh k)))
	       (xk (aref (basis-primal-values b) k)))
	  (setf (aref x (column-ref col)) xk)))
      (dotimes (j n)
	(let ((col (adjvector-column-ref (lp-columns lp) j))
	      (flag (aref (basis-column-flags b) j)))
	  (cond ((eq flag 'basic))
		((eq flag 'nonbasic-lower-bound)
		 (assert (column-has-l col))
		 (setf (aref x j) (column-l col)))
		((eq flag 'nonbasic-upper-bound)
		 (assert (column-has-u col))
		 (setf (aref x j) (column-u col)))
		(t
		 (error "not in phase 2? ~A" flag)))))
      ;; check if rows add up
      (dotimes (row-ref (adjvector-row-fill-pointer (lp-rows lp)))
	(let* ((row (adjvector-row-ref (lp-rows lp) row-ref))
	       (total 0))
	  (when (row-is-active row)
	    (assert (= (adjvector-fixnum-fill-pointer (row-col-refs row)) 
		       (adjvector-fixnum-fill-pointer (row-col-indices row))))
	    (dotimes (k (adjvector-fixnum-fill-pointer (row-col-refs row)))
	      (let* ((j (adjvector-fixnum-ref (row-col-refs row) k))
		     (col (adjvector-column-ref (lp-columns lp) j))
		     (aij (rational-in-column col (adjvector-fixnum-ref (row-col-indices row) k))))
		(incf total (* aij (aref x j)))))
	    (unless (zerop total)
	      (print total)
	      (print x)
	      (print row)))
	  (assert (zerop total)))))))
  
  

  

;;;;
(defun check-phase1-objective (sd)
  (when *checks*
    (let* ((lp (simplex-lp sd))
	   (b (simplex-basis sd))
	   (flags (basis-column-flags b))
	   (n (adjvector-column-fill-pointer (lp-columns lp)))
	   (bh (basis-header b))
	   (z 0))
      (dotimes (j n)
	(let ((col (adjvector-column-ref (lp-columns lp) j))
	      (flag (aref flags j)))
	  (cond ((zerop (column-c col)))
		((and (eq flag 'nonbasic-upper-bound)
		      (not (column-has-u col)))
		 (incf z (* (- (lp-obj-sense lp)) (column-c col))))
		((and (eq flag 'nonbasic-lower-bound)
		      (not (column-has-l col)))
		 (decf z (* (- (lp-obj-sense lp)) (column-c col)))))))
      (dotimes (k (length bh))
	(let* ((j (aref bh k))
	       (v (aref (basis-primal-values b) k))
	       (col (adjvector-column-ref (lp-columns lp) j)))
	  (assert (eq 'basic (aref flags j)))
	  (incf z (* (- (lp-obj-sense lp)) (column-c col) v))))
      (assert (= z (basis-obj-value b))))))
  


;;;;  
(defun check-pivot-row (sd)
  (when *checks*
    (let* ((lp (simplex-lp sd))
	   (b (simplex-basis sd))
	   (m (basis-matrix-size (basis-matrix b)))
	   (flags (basis-column-flags b))
	   (n (length flags))
	   (rho (tran-hsv (simplex-btran sd)))
	   (vc (make-array m :initial-element 0 :element-type 'rational))
	   (alpha (make-array n :initial-element 0 :element-type 'rational)))
      (dotimes (k (simplex-pivot-row-length sd))
	(setf (aref alpha (aref (simplex-pivot-row-col-refs sd) k))
	      (aref (simplex-pivot-row-values sd) k)))
      (dotimes (j n)
	(let* ((col (adjvector-column-ref (lp-columns lp) j))
	       (flag (aref flags j))
	       (x 0))
	  (when (or (eq flag 'nonbasic-lower-bound)
		    (eq flag 'nonbasic-upper-bound))
	  (dotimes (i m)
	    (setf (aref vc i) 0))
	  (dotimes (k (hsv-length (column-hsv col)))
	    (let* ((row-ref (aref (hsv-is (column-hsv col)) k))
		   (row-ind (adjvector-fixnum-ref (lp-active-row-inds lp) row-ref)))
	      (unless (= -1 row-ind)
		(setf (aref vc row-ind)
		      (aref (hsv-vis (column-hsv col)) k)))))
	  (dotimes (k (hsv-length rho))
	    (let* ((row-index (aref (hsv-is rho) k))
		   (row-value (aref (hsv-vis rho) k)))
	      (incf x (* (aref vc row-index) row-value))))
	  (mulf x (hsv-coef (column-hsv col)))
	  (mulf x (hsv-coef rho))
	  (assert (= x (aref alpha j)))))))))



;;;;
(defun check-dse-weights (sd)
  (when *checks*
    (let ((tr (simplex-btran sd))
	  (beta-coef (basis-dse-coef (simplex-basis sd)))
	  (beta-vis (basis-dse-weight-vis (simplex-basis sd)))
	  (m (basis-matrix-size (basis-matrix (simplex-basis sd))))
	  (rhs (simplex-hsv sd))
	  (test nil))
      (dotimes (i m t)
	(reset-hsv rhs)
	(hsv-add i 1 rhs)
	(btran tr rhs)
	(let ((weight 0)
	      (betai (* beta-coef (aref beta-vis i)))
	    (rho (tran-hsv tr)))
	  (dotimes (k (hsv-length rho))
	    (let ((rhok (aref (hsv-vis rho) k)))
	      (incf weight (* rhok rhok))))
	  (mulf weight (hsv-coef rho))
	  (mulf weight (hsv-coef rho))
	  (unless (= weight betai)
	    (setf test t)
	    (return))))
      (when test
	(check-lu (simplex-lp sd)
		  (basis-matrix (simplex-basis sd)) 
		  (basis-header (simplex-basis sd)))
	(dotimes (i m t)
	  (reset-hsv rhs)
	  (hsv-add i 1 rhs)
	  (btran tr rhs)
	  (check-btran (simplex-basis sd) (simplex-lp sd) tr rhs)
	  (let ((weight 0)
		(betai (* beta-coef (aref beta-vis i)))
		(rho (tran-hsv tr)))
	    (dotimes (k (hsv-length rho))
	      (let ((rhok (aref (hsv-vis rho) k)))
		(incf weight (* rhok rhok))))
	    (mulf weight (hsv-coef rho))
	    (mulf weight (hsv-coef rho))
	    (assert (= weight betai))))))))
  

	

(defun compute-dense-pivot-row (sd)
  (let ((dpr (make-array (length (simplex-pivot-row-values sd)) :element-type 'rational)))
    (dotimes (k (hsv-length (simplex-hsv sd)) dpr)
      (let* ((row-index (aref (hsv-is (simplex-hsv sd)) k))
	     (row-value (aref (hsv-vis (simplex-hsv sd)) k))
	     (row-ref (adjvector-fixnum-ref (lp-active-row-refs (simplex-lp sd)) row-index))
	     (row (adjvector-row-ref (lp-rows (simplex-lp sd)) row-ref)))
	(dotimes (row-k (adjvector-fixnum-fill-pointer (row-col-refs row)))
	  (let* ((col-ref (adjvector-fixnum-ref (row-col-refs row) row-k))
		 (col (adjvector-column-ref (lp-columns (simplex-lp sd)) col-ref))
		 (col-value (aref (hsv-vis (column-hsv col)) (adjvector-fixnum-ref (row-col-indices row) row-k)))
		 (flag (aref (basis-column-flags (simplex-basis sd)) col-ref)))
	    (unless (eq flag 'basic)
	      (incf (aref dpr col-ref)
		    (* (hsv-coef (simplex-hsv sd))
		       (hsv-coef (column-hsv col))
		       row-value
		       col-value)))))))))


  

(defun check-flip-rhs (sd)
  (when *checks*
    (let* ((lp (simplex-lp sd))
	   (inphase1 (basis-in-phase1 (simplex-basis sd)))
	   (m (basis-matrix-size (basis-matrix (simplex-basis sd))))
	   (v1 (make-array m :element-type 'rational))
	   (v2 (make-array m :element-type 'rational)))
      (dotimes (k (hsv-length (simplex-hsv sd)))
	(setf (aref v1 (aref (hsv-is (simplex-hsv sd)) k))
	      (* (hsv-coef (simplex-hsv sd)) (aref (hsv-vis (simplex-hsv sd)) k))))
      (dotimes (kj (simplex-n-flips sd))
	(let* ((j (aref (simplex-flip-col-refs sd) kj))
	       (flag (aref (basis-column-flags (simplex-basis sd)) j))
	       (col (adjvector-column-ref (lp-columns lp) j))
	       (col-coef (* (hsv-coef (column-hsv col))
			    (if (eq flag 'nonbasic-lower-bound)
				(column-u-minus-l inphase1 col)
				(column-l-minus-u inphase1 col)))))
	  (dotimes (k (hsv-length (column-hsv col)))
	    (incf (aref v2 (aref (hsv-is (column-hsv col)) k))
		  (* col-coef (aref (hsv-vis (column-hsv col)) k))))))
      (unless (equalp v1 v2)
	(print v1)
	(print 'should-be)
	(print v2))
      (assert (equalp v1 v2)))))
  




