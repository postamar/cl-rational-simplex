(defparameter *checks* nil)


(defun make-dense-basis (lp bm bh)
  (let* ((m (basis-matrix-size bm))
	 (db (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (k m db)
      (let* ((col-ref (aref bh k))
	     (col (aref (lp-columns lp) col-ref)))
	(dotimes (l (length (column-row-refs col)))
	  (let* ((row-ref (aref (column-row-refs col) l))
		 (row-ind (aref (lp-active-row-inds lp) row-ref)))
	    (unless (= -1 row-ind)
	      (setf (aref db row-ind k)
		    (* (column-coef col) (aref (column-values col) l))))))))))



;;;; Verifies the LU decomposition
(defun lu-check (lp bm bh)
  (when *checks*
    (let* ((m (basis-matrix-size bm))
	   (ua (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (la (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (ta (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (da (make-dense-basis lp bm bh)))
      ;; fill u
      (dotimes (j (length (basis-matrix-u-columns bm)))
	(let* ((u (aref (basis-matrix-u-columns bm) j)))
	  (dotimes (r (length (hsv-is u)))
	    (let ((i (aref (hsv-is u) r)))
	      (setf (aref ua i j)
		    (* (hsv-coef u) (aref (hsv-vis u) r)))))))
      ;; compute
      (dotimes (k (basis-matrix-n-l-file bm))
	;; reset l
	(let* ((l (aref (basis-matrix-l-file bm) k))
	       (lj (hsv-j l)))
	  (dotimes (i m)
	    (dotimes (j m)
	      (setf (aref la i j) (if (and (= i j) (/= j lj)) 1 0))))
	  (dotimes (kl (length (hsv-is l)))
	    (setf (aref la (aref (hsv-is l) kl) lj)
		  (* (hsv-coef l) (aref (hsv-vis l) kl)))))
	;; matrix multiplication
	(dotimes (i m)
	  (dotimes (j m)
	    (let ((v 0))
	      (dotimes (km m)
		(incf v (* (aref la i km) (aref da km j))))
	      (setf (aref ta i j) v))))
	;; matrix copy
	(dotimes (i m)
	  (dotimes (j m)
	    (setf (aref da i j) (aref ta i j)))))
      ;; check-
      (dotimes (i m)
	(dotimes (j m)
	  (assert (= (aref da i j) (aref ua i j)))))
      ;; check orders
      (dotimes (j m)
	(let ((u (aref (basis-matrix-u-columns bm) j))
	      (u-seq (aref (basis-matrix-u-seqs bm) j)))
	  (dotimes (k (- (hsv-length u) 1))
	  (assert (< (aref (basis-matrix-i->pi bm) (aref (hsv-is u) (aref u-seq k)))
		     (aref (basis-matrix-i->pi bm) (aref (hsv-is u) (aref u-seq (+ k 1))))))))))))
  
  
  

;;;;
(defun check-btran (sd)
  (when *checks*
    (let* ((db (make-dense-basis (simplex-lp sd) (basis-matrix (simplex-basis sd))
				 (basis-header (simplex-basis sd))))
	   (tr (simplex-tran sd))
	   (m (basis-matrix-size (basis-matrix (simplex-basis sd))))
	   (vs (make-nvector m 0 rational))
	   (vr (make-nvector m 0 rational))
	   (vt (make-nvector m 0 rational)))
      (dotimes (k (length (simplex-vector-indices sd)))
	(let ((i (aref (simplex-vector-indices sd) k)))
	  (setf (aref vs i) (* (simplex-vector-coef sd)
			       (aref (simplex-vector-values sd) k)))))
      (dotimes (k (length (tran-indices (simplex-tran sd))))
	(setf (aref vr (aref (tran-indices tr) k))
	      (* (tran-coef tr) (aref (tran-values tr) k))))
      (dotimes (j m)
	(dotimes (i m)
	  (incf (aref vt j)
		(* (aref vr i)
		   (aref db i j)))))
					;   (print (basis-header (simplex-basis sd)))
      #|
      (print vr)
      (print db)
      (print vt)
      (print vs)
      (print '---)
      |#    (dotimes (i m t)
	      (assert (= (aref vt i) (aref vs i)))))))
  

(defun check-ftran (sd)
  (when *checks*
    (let* ((db (make-dense-basis (simplex-lp sd) (basis-matrix (simplex-basis sd))
				 (basis-header (simplex-basis sd))))
	   (tr (simplex-tran sd))
	   (m (basis-matrix-size (basis-matrix (simplex-basis sd))))
	   (vs (make-nvector m 0 rational))
	   (vr (make-nvector m 0 rational))
	   (vt (make-nvector m 0 rational)))
      (dotimes (k (length (simplex-vector-indices sd)))
	(let ((i (aref (simplex-vector-indices sd) k)))
	  (setf (aref vs i) (* (simplex-vector-coef sd)
			       (aref (simplex-vector-values sd) k)))))
      (dotimes (k (length (tran-indices (simplex-tran sd))))
	(setf (aref vr (aref (tran-indices tr) k))
	      (* (tran-coef tr) (aref (tran-values tr) k))))
      (dotimes (i m)
	(dotimes (j m)
	  (incf (aref vt i)
		(* (aref vr j)
		   (aref db i j)))))
					;    (print (list vr db vt vs))
      (dotimes (i m t)
	(assert (= (aref vt i) (aref vs i)))))))
  




(defun print-u-file (bm)
  (let* ((m (basis-matrix-size bm))
	 (ua (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m) 
      (let ((u (aref (basis-matrix-u-columns bm) j)))
	(dotimes (k (length (hsv-is u)))
	  (setf (aref ua 
		      (aref (basis-matrix-i->pi bm) (aref (hsv-is u) k))
		      (aref (basis-matrix-j->pj bm) j))
		(* (hsv-coef u) (aref (hsv-vis u) k))))))
    (print ua)
    t))


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
  
  
(defun check-primal-feasability (sd)
  (when *checks*
    (let* ((b (simplex-basis sd))
	   (bh (basis-header b))
	   (lp (simplex-lp sd)))
      (check-primal-values sd)
      (dotimes (k (length bh))
	(let* ((col (aref (lp-columns lp) (aref bh k)))
	       (xk (aref (basis-primal-values b) k)))
	  (when (column-has-l col)
	    (assert (<= (column-l col) xk)))
	  (when (column-has-u col)
	    (assert (<= xk (column-u col)))))))))
  
  

(defun check-primal-values (sd)
  (when *checks*
    (let* ((b (simplex-basis sd))
	   (bh (basis-header b))
	   (lp (simplex-lp sd))
	   (x (make-nvector (length (lp-columns lp)) 0 rational)))
      ;; fill x
      (dotimes (k (length bh))
	(let* ((col (aref (lp-columns lp) (aref bh k)))
	       (xk (aref (basis-primal-values b) k)))
	  (setf (aref x (column-ref col)) xk)))
      (dotimes (j (length (lp-columns lp)))
	(let ((col (aref (lp-columns lp) j))
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
      (dotimes (row-ref (length (lp-rows lp)))
	(let* ((row (aref (lp-rows lp) row-ref))
	       (total 0))
	  (when (row-is-active row)
	    (assert (= (length (row-col-refs row)) (length (row-col-indices row))))
	    (dotimes (k (length (row-col-refs row)))
	      (let* ((j (aref (row-col-refs row) k))
		     (col (aref (lp-columns lp) j))
		     (aij (rational-in-column col (aref (row-col-indices row) k))))
		(incf total (* aij (aref x j)))))
	    (unless (zerop total)
	      (print total)
	      (print x)
	      (print row)))
	  (assert (zerop total)))))))
  
  
(defun check-reduced-costs (sd)
  (when *checks*
    (let ((orig-rcosts (copy-seq (basis-reduced-costs (simplex-basis sd))))
	  (orig-prow (copy-seq (simplex-pivot-row sd))))
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
	  (print '-)
	  (print (simplex-dual-step sd))
	  (print orig-prow)
	  (print orig-rcosts)
	  (print '--should-be)
	  (print (basis-reduced-costs (simplex-basis sd)))
	  (error "reduced costs"))))))


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
(defun check-phase1-objective (sd)
  (when *checks*
    (let* ((lp (simplex-lp sd))
	   (b (simplex-basis sd))
	   (flags (basis-column-flags b))
	   (n (length (lp-columns lp)))
	   (bh (basis-header b))
	   (z 0))
      (dotimes (j n)
	(let ((col (aref (lp-columns lp) j))
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
	       (col (aref (lp-columns lp) j)))
	  (assert (eq 'basic (aref flags j)))
	  (incf z (* (- (lp-obj-sense lp)) (column-c col) v))))
      (assert (= z (basis-obj-value b))))))
  
  

(defun check-pivot-row (sd)
  (when *checks*
    (let* ((lp (simplex-lp sd))
	   (b (simplex-basis sd))
	   (m (basis-matrix-size (basis-matrix b)))
	   (vc (make-nvector m 0 rational))
	   (flags (basis-column-flags b)))
      (dotimes (j (length flags))
	(let* ((col (aref (lp-columns lp) j))
	       (flag (aref flags j))
	       (x 0))
	  (when (or (eq flag 'nonbasic-lower-bound)
		    (eq flag 'nonbasic-upper-bound))
	  (dotimes (i m)
	    (setf (aref vc i) 0))
	  (dotimes (k (length (column-row-refs col)))
	    (let* ((row-ref (aref (column-row-refs col) k))
		   (row-ind (aref (lp-active-row-inds lp) row-ref)))
	      (unless (= -1 row-ind)
		(setf (aref vc row-ind)
		      (aref (column-values col) k)))))
	  (dotimes (k (length (simplex-vector-indices sd)))
	    (let* ((row-index (aref (simplex-vector-indices sd) k))
		   (row-value (aref (simplex-vector-values sd) k)))
	      (incf x (* (aref vc row-index) row-value))))
	  (mulf x (column-coef col))
	  (mulf x (simplex-vector-coef sd))
	  (assert (= x (aref (simplex-pivot-row sd) j)))))))))
