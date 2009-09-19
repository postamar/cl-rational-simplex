;;;; 


(defstruct (simplex
	     (:constructor %make-simplex))
  (lp               nil :type lp)
  (basis            nil :type basis)
  (tran             nil :type tran)
  (dse-tran         nil :type tran)
  (vector-coef      1   :type rational)
  (vector-indices   #() :type vector)
  (vector-values    #() :type vector)
  (pivot-row        #() :type vector)
  (n-breakpoints    0   :type fixnum)
  (breakpoints      #() :type bit-vector)
  (pivot-row-index  -1  :type fixnum)
  (pivot-col-ref    -1  :type fixnum)
  (delta            0   :type rational)
  (dual-step        0   :type rational)
  (primal-step      0   :type rational)
  (flip-col-indices #() :type vector))


(defun make-simplex (lp basis)
  (let ((n (length (lp-columns lp)))
	(m (length (lp-rows lp))))
    (%make-simplex
     :lp               lp
     :basis            basis
     :tran             (make-tran (basis-matrix basis))
     :dse-tran             (make-tran (basis-matrix basis))
     :vector-indices   (make-nvector m -1 fixnum)
     :vector-values    (make-nvector m 0 integer)
     :pivot-row        (make-nvector n 0 rational)
     :breakpoints      (make-nvector n 0 bit)
     :flip-col-indices (make-nvector n -1 fixnum))))
    

;;;; Resets the data structure to start new iteration
(defun reset-simplex (sd)
  (setf (simplex-vector-coef sd) 1
	(fill-pointer (simplex-vector-indices sd)) 0
	(fill-pointer (simplex-vector-values sd)) 0
	(simplex-n-breakpoints sd) 0
	(simplex-pivot-row-index sd) -1
	(simplex-pivot-col-ref sd) -1
	(simplex-delta sd) 0
	(simplex-dual-step sd) 0
	(simplex-primal-step sd) 0
	(fill-pointer (simplex-flip-col-indices sd)) 0)
  (bit-xor (simplex-breakpoints sd) (simplex-breakpoints sd) t)
  (let ((n (length (simplex-pivot-row sd))))
    (dotimes (j n)
      (setf (aref (simplex-pivot-row sd) j) 0))))
 

;;;; Copies TRAN vector into simplex data structure
(defun copy-vector-from-tran (sd)
  (let* ((tr (simplex-tran sd))
	 (n (length (tran-indices tr))))
    (setf (simplex-vector-coef sd) (tran-coef tr)
	  (fill-pointer (simplex-vector-indices sd)) n
	  (fill-pointer (simplex-vector-values sd))  n)
    (dotimes (k n)
      (setf (aref (simplex-vector-indices sd) k) (aref (tran-indices tr) k)
	    (aref (simplex-vector-values sd) k)  (aref (tran-values tr) k)))))



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
	 (m (length (basis-header b)))
	 (best-val 0))
    (dotimes (k m)
      (let* ((weight (aref (basis-dse-weights b) k))
	     (col (aref (lp-columns lp) (aref (basis-header b) k)))
	     (x (aref (basis-primal-values b) k))
	     (delta (get-delta (basis-in-phase1 b) x col))
	     (val (/ (* delta delta) weight)))
	(when (< best-val val)
	  (setf best-val val
		(simplex-delta sd) delta
		(simplex-pivot-row-index sd) k))))))


;;;;
(defun simplex-btran (sd)
  (btran (simplex-tran sd) 
	 (basis-matrix-i->pi (basis-matrix (simplex-basis sd)))
	 (basis-matrix-j->pj (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pi->i (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pj->j (basis-matrix (simplex-basis sd)))
	 (simplex-vector-coef sd)
	 (simplex-vector-indices sd)
	 (simplex-vector-values sd)))


;;;;
(defun simplex-dse-ftran (sd)
  (ftran (simplex-dse-tran sd)
	 (basis-matrix-i->pi (basis-matrix (simplex-basis sd)))
	 (basis-matrix-j->pj (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pi->i (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pj->j (basis-matrix (simplex-basis sd)))
	 (tran-coef (simplex-tran sd))
	 (tran-indices (simplex-tran sd))
	 (tran-values (simplex-tran sd))))






;;;; Row-wise pivot row computation
(defun compute-pivot-row (sd)
  (dotimes (k (length (simplex-vector-indices sd)))
    (let* ((row-index (aref (simplex-vector-indices sd) k))
	   (row-value (aref (simplex-vector-values sd) k))
	   (row-ref (aref (lp-active-row-refs (simplex-lp sd)) row-index))
	   (row (aref (lp-rows (simplex-lp sd)) row-ref)))
      (dotimes (row-k (length (row-col-refs row)))
	(let* ((col-ref (aref (row-col-refs row) row-k))
	       (col (aref (lp-columns (simplex-lp sd)) col-ref))
	       (col-value (aref (column-values col) (aref (row-col-indices row) row-k)))
	       (flag (aref (basis-column-flags (simplex-basis sd)) col-ref)))
	(unless (eq flag 'basic)
	  (incf (aref (simplex-pivot-row sd) col-ref)
		(* (simplex-vector-coef sd)
		   (column-coef col)
		   row-value
		   col-value))))))))


		       


;;;;
(defun choose-entering-basis-index (sd)
  (let* ((r (simplex-pivot-row-index sd))
	 (b (simplex-basis sd))
	 (best-val 0)
	 (n (length (simplex-pivot-row sd)))
	 (p (aref (basis-header b) r))
	 (x (aref (basis-primal-values b) r))
	 (col (aref (lp-columns (simplex-lp sd)) p))
	 (sign (if (basis-in-phase1 b)
		   (if (column-has-l col)
		       (if (< x 0) -1 1)
		       (if (< x -1) -1 1))
		   (if (and (column-has-l col) (< x (column-l col))) -1 1))))
    ;; perform pricing
    (dotimes (j n)
      (let ((col (aref (lp-columns (simplex-lp sd)) j))
	    (alpha (* sign (aref (simplex-pivot-row sd) j)))
	    (flag (aref (basis-column-flags (simplex-basis sd)) j)))
	(when (and (or (eq flag 'nonbasic-lower-bound)
		       (eq flag 'nonbasic-upper-bound))
		   (or (and (not (column-has-l col))
			    (not (column-has-u col)))
		       (and (eq flag 'nonbasic-lower-bound)
			    (< 0 alpha))
		       (and (eq flag 'nonbasic-upper-bound)
			    (< alpha 0))))
	  (let* ((rcost (aref (basis-reduced-costs (simplex-basis sd)) j))
		 (val (/ rcost alpha)))
	    (when (or (= -1 (simplex-pivot-col-ref sd))
		      (< val best-val))
	      (setf (simplex-pivot-col-ref sd) j
		    best-val val)))
	  (setf (bit (simplex-breakpoints sd) j) 1)
	  (incf (simplex-n-breakpoints sd)))))
    (unless (= -1 (simplex-pivot-col-ref sd))
      (setf (simplex-dual-step sd) 
	    (/ (aref (basis-reduced-costs b) (simplex-pivot-col-ref sd))
	       (aref (simplex-pivot-row sd) (simplex-pivot-col-ref sd)))))
    (simplex-pivot-col-ref sd)))



;;;; 
(defun set-entering-column-as-simplex-vector (sd)
  (let* ((q (simplex-pivot-col-ref sd))
	 (lp (simplex-lp sd))
	 (indices (simplex-vector-indices sd))
	 (values (simplex-vector-values sd))
	 (col (aref (lp-columns lp) q)))
    (setf (fill-pointer indices) 0
	  (fill-pointer values) 0
	  (simplex-vector-coef sd) (column-coef col))
    (dotimes (k (length (column-row-refs col)))
      (let ((ind (aref (lp-active-row-inds lp) (aref (column-row-refs col) k))))
	(unless (= -1 ind)
	  (vector-push ind indices)
	  (vector-push (aref (column-values col) k) values))))))




;;;;
(defun simplex-ftran (sd)
  (ftran (simplex-tran sd) 
	 (basis-matrix-i->pi (basis-matrix (simplex-basis sd)))
	 (basis-matrix-j->pj (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pi->i (basis-matrix (simplex-basis sd)))
	 (basis-matrix-pj->j (basis-matrix (simplex-basis sd)))
	 (simplex-vector-coef sd)
	 (simplex-vector-indices sd)
	 (simplex-vector-values sd)))



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
	 (dse-tr (simplex-dse-tran sd))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-row-weight (aref (basis-dse-weights b) exit-row))
	 (primal-values (basis-primal-values b))
	 (alpha-index (find-index (simplex-vector-indices sd) exit-row)))
    (assert (/= -1 alpha-index))
    (assert (/= 0 (aref (simplex-vector-values sd) alpha-index)))
    ;; compute primal step
    (setf (simplex-primal-step sd)
	  (/ (simplex-delta sd)
	     (* (aref (simplex-vector-values sd) alpha-index)
		(simplex-vector-coef sd))))
    (dotimes (k (length (simplex-vector-indices sd)))
      (let ((i (aref (simplex-vector-indices sd) k)))
	;; update basic primal values
	(decf (aref primal-values i) 
	      (* (simplex-primal-step sd)
		 (simplex-vector-coef sd)
		 (aref (simplex-vector-values sd) k)))
	;; update dse weights
	(if (= k alpha-index)
	    (let ((d (* (simplex-vector-coef sd)
			(aref (simplex-vector-values sd) alpha-index))))
	      (assert (= exit-row-weight (aref (basis-dse-weights b) exit-row)))
	      (divf (aref (basis-dse-weights b) exit-row) (* d d)))
	    (let ((ratio (/ (aref (simplex-vector-values sd) k)
			    (aref (simplex-vector-values sd) alpha-index)))
		  (tau-index (find-index (tran-indices dse-tr) i)))
	      (unless (= -1 tau-index)
		(decf (aref (basis-dse-weights b) i)
		      (* 2 ratio 
			 (tran-coef dse-tr) 
			 (aref (tran-values dse-tr) tau-index))))
	      (incf (aref (basis-dse-weights b) i)
		    (* ratio ratio exit-row-weight))))))))



;;;;
(defun simplex-basis-matrix-update (sd)
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (bm (basis-matrix b)))
    (basis-matrix-lu-decomposition bm)
    (lu-check lp bm (basis-header b))))
  


;;;;	     
(defun simplex-basis-update (sd)
  (let* ((b (simplex-basis sd))
	 (header (basis-header b))
	 (lp (simplex-lp sd))
	 (n (length (lp-columns lp)))
	 (q (simplex-pivot-col-ref sd))
	 (exit-row (simplex-pivot-row-index sd))
	 (exit-col-ref (aref header exit-row))
	 (primal-values (basis-primal-values b))
	 (flags (basis-column-flags b))
	 (rcosts (basis-reduced-costs b)))
    ;; compute primal step, update primal basic values, update dse weights
    (simplex-basis-update-primal sd)
    ;; update reduced costs
    (setf (aref rcosts exit-col-ref) (- (simplex-dual-step sd)))
    (dotimes (j n)
      (let ((flag (aref flags j)))
	(when (or (eq flag 'nonbasic-lower-bound)
		  (eq flag 'nonbasic-upper-bound))
	  (decf (aref rcosts j)
		(* (simplex-dual-step sd) 
		   (aref (simplex-pivot-row sd) j))))))
    ;; update exiting variable
    (setf (aref flags exit-col-ref)
	  (basis-exiting-variable-flag (basis-in-phase1 b) 
				       (aref primal-values exit-row)
				       (aref rcosts exit-col-ref)
				       (aref (lp-columns lp) exit-col-ref)))
    ;; update entering variable
    (setf (aref header exit-row) q)
    (setf (aref primal-values exit-row) 
	  (+ (simplex-primal-step sd)
	     (nonbasic-value (basis-in-phase1 b)
			     (aref (lp-columns lp) q)
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
  (vector-push (simplex-pivot-row-index sd) (simplex-vector-indices sd))
  (vector-push 1 (simplex-vector-values sd))
  (simplex-btran sd)
  (check-btran sd)
  (copy-vector-from-tran sd)
  ;; pivot row
  (compute-pivot-row sd)
  (check-pivot-row sd)
  ;; ratio test
  (when (= -1 (choose-entering-basis-index sd))
    (return-from simplex-iteration 'infeasible))
  ;; ftrans
  (simplex-dse-ftran sd)
  (set-entering-column-as-simplex-vector sd)
  (simplex-ftran sd)
  (check-ftran sd)
  (copy-vector-from-tran sd)
  ;; basis change and update
  (simplex-basis-update sd)
  (check-reduced-costs sd)
  (check-dual-feasability sd)
  'dual-feasible)


