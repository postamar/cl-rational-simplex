;;;; 

(symbol-macrolet 
    ((err (error "simplex constructor")))
  (defstruct (simplex
	       (:constructor %make-simplex)
	       (:print-function print-simplex))
    (lp                  err :type lp)
    (basis               err :type basis)
    (btran               err :type tran)
    (ftran               err :type tran)
    (flip-ftran          err :type tran)
    (dse-ftran           err :type tran)
    (hsv                 err :type hsv)
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
    (n-flips             0   :type fixnum)))



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

(defun make-simplex (lp basis)
  (let ((n (adjvector-column-fill-pointer (lp-columns lp))))
    (%make-simplex
     :lp                 lp
     :basis              basis
     :btran              (make-tran (basis-matrix basis))
     :ftran              (make-tran (basis-matrix basis))
     :flip-ftran         (make-tran (basis-matrix basis))
     :dse-ftran          (make-tran (basis-matrix basis))
     :hsv                (make-hsv)
     :pivot-row-col-refs (make-array n :initial-element -1 :element-type 'fixnum)
     :pivot-row-values   (make-array n :initial-element 0 :element-type 'rational)
     :breakpoints        (make-array n :initial-element 0 :element-type 'bit)
     :add-heap           (make-array n :initial-element -1 :element-type 'fixnum)
     :add-counters       (make-array n :initial-element -1 :element-type 'fixnum)
     :add-refs           (make-array n :initial-element -1 :element-type 'fixnum)
     :flip-col-refs      (make-array n :initial-element -1 :element-type 'fixnum)
     :flip-col-coefs     (make-array n :initial-element 0 :element-type 'rational))))
    

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
 



(defun get-delta (in-phase1 x col)
  (if in-phase1
      (cond ((and (column-has-l col)
		  (< x 0))
	     x)
	    ((and (column-has-u col)
		  (< 0 x))
	     x)
	    (t 0))
      (cond ((and (column-has-l col)
		  (< x (column-l col)))
	     (- x (column-l col)))
	    ((and (column-has-u col)
		  (< (column-u col) x))
	     (- x (column-u col)))
	    (t 0))))



;;;; Chooses exiting variable
(defun choose-exiting-basis-index (sd)
  (let* ((lp (simplex-lp sd))
	 (b (simplex-basis sd))
	 (m (basis-matrix-size (basis-matrix b)))
	 (best-k '())
	 (infty nil)
	 (best-val 0))
    (dotimes (k m)
      (let* ((weight (aref (basis-dse-weights b) k))
	     (col (adjvector-column-ref (lp-columns lp) (aref (basis-header b) k)))
	     (x (aref (basis-primal-values b) k))
	     (delta (get-delta (basis-in-phase1 b) x col)))
	(cond ((< weight 0)
	       (error "negative dse weight"))
	      ((zerop delta))
	      ((and (zerop weight) infty)
	       (push k best-k))
	      ((zerop weight)
	       (setf best-k (list k)
		     infty t))
	      (infty)
	      (t
	       (let ((val (/ (* delta delta) weight)))
		 (cond
		   ((< best-val val)
		    (setf best-val val
			  best-k (list k)))
		   ((= best-val val)
		    (push k best-k))))))))
    (if (and (not infty) (zerop best-val))
	(setf (simplex-delta sd) 0
	      (simplex-pivot-row-index sd) -1)
	(let ((choice-k (elt best-k (random (length best-k)))))
	  (setf (simplex-delta sd)
		(get-delta (basis-in-phase1 b)
		     (aref (basis-primal-values b) choice-k)
		     (adjvector-column-ref (lp-columns lp) 
					   (aref (basis-header b) choice-k)))
		(simplex-pivot-row-index sd) 
		choice-k)))))
	  


;;;; Row-wise pivot row computation
(defun compute-pivot-row (sd)
  (let* ((rho (tran-hsv (simplex-btran sd)))
	 (n-rows (hsv-length rho))
	 (last-heap-index (- n-rows 1))
	 (heap (simplex-add-heap sd))
	 (prev-col-ref -1))
    (labels ((get-col-ref (k)
	       (let* ((counter (aref (simplex-add-counters sd) k))
		      (row-ref (aref (simplex-add-refs sd) k))
		      (row (adjvector-row-ref (lp-rows (simplex-lp sd)) row-ref))
		      (row-col-refs  (row-col-refs row)))
		 (if (= counter (adjvector-fixnum-fill-pointer row-col-refs))
		     most-positive-fixnum
		     (adjvector-fixnum-ref row-col-refs counter))))
	     (sift-down (root end)
	       (loop
		  (let ((child (+ 1 (* 2 root))))
		    (when (> child end)
		      (return))
		    (when (and (< child end)
			       (> (get-col-ref (aref heap child))
				  (get-col-ref (aref heap (+ child 1)))))
		      (incf child))
		    (when (<= (get-col-ref (aref heap root)) 
			      (get-col-ref (aref heap child)))
		      (return))
		    (rotatef (aref heap root) (aref heap child))
		    (setf root child)))))
      ;; initialization
      (when (zerop n-rows)
	(setf (simplex-pivot-row-length sd) 0)
	(return-from compute-pivot-row))
      (setf (simplex-pivot-row-length sd) -1)
      (dotimes (k n-rows)
	(setf (aref (simplex-add-counters sd) k) 0
	      (aref (simplex-add-heap sd) k) k
	      (aref (simplex-add-refs sd) k) (adjvector-fixnum-ref (lp-active-row-refs (simplex-lp sd))
								   (aref (hsv-is rho) k))))
      ;; heapify
      (loop for heapify-start from (floor (- last-heap-index 1) 2) downto 0
	 do (sift-down heapify-start last-heap-index))
      ;; build pivot row
      (loop
	 ;; pop element with smallest col-ref off heap
	 (let* ((smallest-k (aref heap 0))
		(row (adjvector-row-ref (lp-rows (simplex-lp sd)) (aref (simplex-add-refs sd) smallest-k)))
		(col-ref (get-col-ref smallest-k)))
	   (cond ((= col-ref most-positive-fixnum)
		  (incf (simplex-pivot-row-length sd))
		  (return))
		 ((> prev-col-ref col-ref)
		  (print row)
		  (error "bad col-ref order"))
		 ((= col-ref prev-col-ref))
		 (t
		  (setf prev-col-ref col-ref)
		  (incf (simplex-pivot-row-length sd))
		  (setf (aref (simplex-pivot-row-col-refs sd) 
			      (simplex-pivot-row-length sd)) col-ref)
		  (setf (aref (simplex-pivot-row-values sd)
			      (simplex-pivot-row-length sd)) 0)))
	   ;; update pivot row value
	   (incf (aref (simplex-pivot-row-values sd)
		       (simplex-pivot-row-length sd))
		 (* (aref (hsv-vis (column-hsv (adjvector-column-ref (lp-columns (simplex-lp sd)) col-ref)))
			  (adjvector-fixnum-ref (row-col-indices row) (aref (simplex-add-counters sd) smallest-k)))
		    (aref (hsv-vis rho) smallest-k)))
	   ;; update counter and heap
	   (incf (aref (simplex-add-counters sd) smallest-k))
	   (sift-down 0 last-heap-index)))
      ;; multiply pivot row by coefficients
      (let ((coef1 (hsv-coef rho)))
	(dotimes (k (simplex-pivot-row-length sd))
	  (mulf (aref (simplex-pivot-row-values sd) k)
		(* (hsv-coef (column-hsv (adjvector-column-ref (lp-columns (simplex-lp sd))
							       (aref (simplex-pivot-row-col-refs sd) k))))
		   coef1)))))))
	      
      

;;;;
(defun choose-entering-basis-index (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (best-val 0)
	 (best-k -1)
	 (delta-obj 0)
	 (boxed nil)
	 (slope (abs (simplex-delta sd)))
	 (sign (signum (simplex-delta sd))))
    ;; set breakpoints
    (dotimes (k (simplex-pivot-row-length sd))
      (let* ((j (aref (simplex-pivot-row-col-refs sd) k))
	     (col (adjvector-column-ref (lp-columns (simplex-lp sd)) j))
	     (alpha (* sign (aref (simplex-pivot-row-values sd) k)))
	     (flag (aref (basis-column-flags b) j)))
	(when (or (and (not (column-has-l col))
		       (not (column-has-u col))
		       (not (basis-in-phase1 b))
		       (or (eq flag 'nonbasic-lower-bound)
			   (eq flag 'nonbasic-upper-bound)))
		  (and (eq flag 'nonbasic-lower-bound)
		       (< 0 alpha))
		  (and (eq flag 'nonbasic-upper-bound)
		       (< alpha 0)))
	  (setf (sbit (simplex-breakpoints sd) j) 1)
	  (incf (simplex-n-breakpoints sd)))))
    ;; bound-flipping ratio test
    (when (zerop (simplex-n-breakpoints sd))
      (return-from choose-entering-basis-index -1))
    (loop 
       (setf best-k -1
	     boxed nil)
       ;; pick q
       (dotimes (k (simplex-pivot-row-length sd))
	 (let ((j (aref (simplex-pivot-row-col-refs sd) k)))
	   (unless (zerop (sbit (simplex-breakpoints sd) j))
	     (let* ((col (adjvector-column-ref (lp-columns lp) j))
		    (rcost (aref (basis-reduced-costs (simplex-basis sd)) j))
		    (val (/ rcost (* sign (aref (simplex-pivot-row-values sd) k)))))
	       (when (or (= -1 best-k)
			 (< val best-val)
			 (and (= val best-val)
			      (not boxed)))
		 (setf best-k k
		       best-val val
		       boxed (or (basis-in-phase1 b)
				 (and (column-has-u col) (column-has-l col)))))))))
       ;; update breakpoints and test to exit loop 
       (let* ((q (aref (simplex-pivot-row-col-refs sd) best-k))
	      (colq (adjvector-column-ref (lp-columns lp) q)))
	 (setf (sbit (simplex-breakpoints sd) q) 0)
	 (decf (simplex-n-breakpoints sd))
	 (when (or (not boxed)
		   (zerop (simplex-n-breakpoints sd)))
	   (return))
	 (let ((next-slope (- slope 
			      (* (column-u-minus-l (basis-in-phase1 b) colq)
				 (abs (aref (simplex-pivot-row-values sd) best-k))))))
	   (when (< next-slope 0)
	     (return))
	   (unless (zerop (column-u-minus-l (basis-in-phase1 b) colq))
	     (setf (aref (simplex-flip-col-refs sd) (simplex-n-flips sd)) q)
	     (incf (simplex-n-flips sd))
	     (incf delta-obj
		   (* (column-c colq)
		      (if (eq (aref (basis-column-flags b) q) 'nonbasic-lower-bound)
			  (column-u-minus-l (basis-in-phase1 b) colq)
			  (column-l-minus-u (basis-in-phase1 b) colq))))
	   (setf slope next-slope)))))
    ;; set dual step, set delta, set and return pivot column
    (if (zerop best-val)
	(prog1
	    (setf (simplex-pivot-col-ref sd) 
		  (aref (simplex-pivot-row-col-refs sd) best-k))
	  (setf (simplex-dual-step sd) 0
		(simplex-n-flips sd) 0))
	(prog1 
	    (setf (simplex-pivot-col-ref sd) 
		  (aref (simplex-pivot-row-col-refs sd) best-k))
	  (incf (basis-obj-value b) delta-obj)
	  (setf (simplex-delta sd) (* sign slope))
	  (setf (simplex-dual-step sd) 
		(/ (aref (basis-reduced-costs b) (simplex-pivot-col-ref sd))
		   (aref (simplex-pivot-row-values sd) best-k)))))))
  



;;;;
(defun choose-entering-basis-index-no-bfrt (sd)
  (let* ((b (simplex-basis sd))
	 (best-val 0)
	 (best-cols '())
	 (sign (signum (simplex-delta sd))))
    ;; perform pricing
    (dotimes (k (simplex-pivot-row-length sd))
      (let* ((j (aref (simplex-pivot-row-col-refs sd) k))
	     (col (adjvector-column-ref (lp-columns (simplex-lp sd)) j))
	     (alpha (* sign (aref (simplex-pivot-row-values sd) k)))
	     (flag (aref (basis-column-flags (simplex-basis sd)) j)))
	(when (and (or (eq flag 'nonbasic-lower-bound)
		       (eq flag 'nonbasic-upper-bound))
		   (or (and (not (column-has-l col))
			    (not (column-has-u col))
			    (not (basis-in-phase1 b)))
		       (and (eq flag 'nonbasic-lower-bound)
			    (< 0 alpha))
		       (and (eq flag 'nonbasic-upper-bound)
			    (< alpha 0))))
	  (let* ((rcost (aref (basis-reduced-costs (simplex-basis sd)) j))
		 (val (/ rcost alpha)))
	    (cond ((= val best-val)
		   (push k best-cols))
		  ((or (< val best-val)
		       (eq '() best-cols))
		   (setf best-cols (list k)
			 best-val val))))
	  (incf (simplex-n-breakpoints sd)))))
    ;; select entering variable randomly
    (if (eq '() best-cols)
	-1
	(let* ((entering-col-k (elt best-cols (random (length best-cols))))
	       (entering-col (aref (simplex-pivot-row-col-refs sd) entering-col-k)))
	  (setf (simplex-dual-step sd) 
		(/ (aref (basis-reduced-costs b) entering-col)
		   (aref (simplex-pivot-row-values sd) entering-col-k)))
	  (setf (simplex-pivot-col-ref sd) entering-col)))))




;;;; 
(defun set-entering-column-as-simplex-vector (sd)
  (let* ((q (simplex-pivot-col-ref sd))
	 (lp (simplex-lp sd))
	 (col (adjvector-column-ref (lp-columns lp) q)))
    (reset-hsv (simplex-hsv sd))
    (setf (hsv-coef (simplex-hsv sd)) (hsv-coef (column-hsv col)))
    (dotimes (k (hsv-length (column-hsv col)))
      (let ((ind (adjvector-fixnum-ref (lp-active-row-inds lp)
				       (aref (hsv-is (column-hsv col)) k))))
	(unless (= -1 ind)
	  (hsv-add ind (aref (hsv-vis (column-hsv col)) k) (simplex-hsv sd)))))))



;;;;
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



;;;; 
(defun simplex-basis-update-primal (sd)
  (let* ((b (simplex-basis sd))
	 (dse-tr (simplex-dse-ftran sd))
	 (alphaq (tran-hsv (simplex-ftran sd)))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-row-weight (aref (basis-dse-weights b) exit-row))
	 (primal-values (basis-primal-values b))
	 (alpha-index (hsv-find exit-row alphaq)))
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
	;; update dse weights
	(if (= k alpha-index)
	    (let ((d (* (hsv-coef alphaq)
			(aref (hsv-vis alphaq) alpha-index))))
	      (assert (= exit-row-weight (aref (basis-dse-weights b) exit-row)))
	      (divf (aref (basis-dse-weights b) exit-row) (* d d)))
	    (let ((ratio (/ (aref (hsv-vis alphaq) k)
			    (aref (hsv-vis alphaq) alpha-index)))
		  (tau-index (hsv-find i (tran-hsv dse-tr))))
	      (unless (= -1 tau-index)
		(decf (aref (basis-dse-weights b) i)
		      (* 2 ratio 
			 (hsv-coef (tran-hsv dse-tr))
			 (aref (hsv-vis (tran-hsv dse-tr)) tau-index))))
	      (incf (aref (basis-dse-weights b) i)
		    (* ratio ratio exit-row-weight))))))))


;;;;
(defun simplex-basis-matrix-update (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b)))
    (unless (basis-matrix-lu-decomposition bm)
      (error "basis redundancy"))
    (check-lu lp bm (basis-header b))))
  


;;;; temporarily uses pivot row arrays to compute weighted column sum
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




;;;;	     
(defun simplex-basis-update (sd)
  (let* ((b (simplex-basis sd))
	 (header (basis-header b))
	 (lp (simplex-lp sd))
	 (q (simplex-pivot-col-ref sd))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-col-ref (aref header exit-row))
	 (primal-values (basis-primal-values b))
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
      (check-ftran sd (simplex-flip-ftran sd) (simplex-hsv sd))
      ;; update objective and primal basic values
      (let ((delta-x (tran-hsv (simplex-flip-ftran sd))))
	(dotimes (k (hsv-length delta-x))
	  (let* ((i (aref (hsv-is delta-x) k))
		 (xk (* (hsv-coef delta-x) (aref (hsv-vis delta-x) k)))
		 (cj (column-c (adjvector-column-ref (lp-columns lp) 
						     (aref (basis-header b) i)))))
	    (decf (aref primal-values i) xk)
	    (decf (basis-obj-value b) (* xk cj)))))
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
    (setf (aref flags q) 'basic)
    ;; update objective
    (incf (basis-obj-value b)
	  (* (simplex-dual-step sd)
	     (simplex-delta sd)))
    ;; update basis matrix
    (fill-basis-matrix (basis-matrix b) lp header)
    (simplex-basis-matrix-update sd)))



;;;; 
(defun simplex-iteration (sd)
  ;; reset
  (reset-simplex sd)
  ;; pricing
  (choose-exiting-basis-index sd)
  (when (= -1 (simplex-pivot-row-index sd))
    (return-from simplex-iteration 'optimal))
  ;; btran
  (hsv-add (simplex-pivot-row-index sd) 1 (simplex-hsv sd))
  (btran (simplex-btran sd) (simplex-hsv sd))
  (check-btran sd (simplex-btran sd) (simplex-hsv sd))
  ;; pivot row
  (compute-pivot-row sd)
  (check-pivot-row sd)
  ;; ratio test
  (when (= -1 (choose-entering-basis-index sd))
    (return-from simplex-iteration 'infeasible))
  ;; ftrans
  (ftran (simplex-dse-ftran sd) (tran-hsv (simplex-btran sd)))
  (check-ftran sd (simplex-dse-ftran sd) (tran-hsv (simplex-btran sd)))
  (set-entering-column-as-simplex-vector sd)
  (ftran (simplex-ftran sd) (simplex-hsv sd))
  (check-ftran sd (simplex-ftran sd) (simplex-hsv sd))
  ;; basis change and update
  (simplex-basis-update sd)
  (check-reduced-costs sd)
  (check-dual-feasability sd)
  (check-dse-weights sd)
  'dual-feasible)


