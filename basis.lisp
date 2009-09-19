;;;;; Basis info


(splay-tree :val-type rational)

(defstruct (basis
	     (:constructor %make-basis))
  (matrix         nil :type basis-matrix)
  (in-phase1      t   :type boolean)
  (status         'undef :type t)
  (header         #() :type (simple-array fixnum 1))
  (dse-weights    #() :type (simple-array rational 1))
  (reduced-costs  #() :type (simple-array rational 1))
  (column-flags   #() :type (simple-array symbol 1))
  (obj-value      0   :type rational)
  (primal-values  #() :type (simple-array rational 1))
  (primal-infeas  (error "basis constructor") :type splay-tree-fixnum-rational))


;;;; Sets nonbasic variables for initial basis in phase 1
(defun phase1-initial-nonbasic-flags (lp flags)
  (dotimes (j (adjvector-column-fill-pointer (lp-columns lp)))
    (let ((col (adjvector-column-ref (lp-columns lp) j)))
      (setf (aref flags j)
	    (cond ((not (column-is-active col))
		   'inactive)
		  ((column-is-slack col)
		   'basic)
		  ((and (column-has-l col) (column-has-u col))
		   'boxed)
		  (t 
		   (if (< (* (- (lp-obj-sense lp)) (column-c col)) 0)
		       'nonbasic-upper-bound
		       'nonbasic-lower-bound)))))))


;;;; Computes phase 1 objective from scratch
(defun phase1-initial-objective-value (lp flags)
  (let ((n (adjvector-column-fill-pointer (lp-columns lp)))
	(z 0))
    (dotimes (j n z)
      (let ((col (adjvector-column-ref (lp-columns lp) j))
	    (flag (aref flags j)))
	(cond ((zerop (column-c col)))
	      ((and (eq flag 'nonbasic-upper-bound)
		    (not (column-has-u col)))
	       (incf z (* (- (lp-obj-sense lp)) (column-c col))))
	      ((and (eq flag 'nonbasic-lower-bound)
		    (not (column-has-l col)))
	       (decf z (* (- (lp-obj-sense lp)) (column-c col)))))))))


;;;;
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



;;;; Generates a dual feasible basis using the slack variables
(defun make-phase1-initial-basis (lp)
  (let* ((m (adjvector-fixnum-fill-pointer (lp-active-row-refs lp)))
	 (n (adjvector-column-fill-pointer (lp-columns lp)))
	 (header      (make-array m :initial-element -1 :element-type 'fixnum))
	 (rcosts      (make-array n :initial-element 0 :element-type 'rational))
	 (flags       (make-array n :initial-element 'unknown :element-type 'symbol))
	 (weights     (make-array m :initial-element 1 :element-type 'rational))
	 (values      (make-array m :initial-element 0 :element-type 'rational))
	 (infeas      (make-splay-tree-fixnum-rational)))
    ;; set flags to best bound
    (phase1-initial-nonbasic-flags lp flags)
    ;; set phase1 reduced costs
    (dotimes (j n)
      (let ((flag (aref flags j)))
	(when (or (eq flag 'nonbasic-lower-bound)
		  (eq flag 'nonbasic-upper-bound))
	  (setf (aref rcosts j) (lp-get-cost lp j)))))
    ;; build basis header from slack variables
    ;; and set basic variable values
    (dotimes (i m)
      (let* ((v 0)
	     (row-ref (adjvector-fixnum-ref (lp-active-row-refs lp) i))
	     (row (adjvector-row-ref (lp-rows lp) row-ref))
	     (slack-col-ref (row-slack-col-ref row))
	     (slack-col (adjvector-column-ref (lp-columns lp) slack-col-ref)))
	(dotimes (k (adjvector-fixnum-fill-pointer (row-col-refs row)))
	  (let* ((col-ref (adjvector-fixnum-ref (row-col-refs row) k))
		 (flag (aref flags col-ref))
		 (col (adjvector-column-ref (lp-columns lp) col-ref))
		 (a (rational-in-column col (adjvector-fixnum-ref (row-col-indices row) k))))
	    (cond ((= col-ref slack-col-ref))
		  ((and (eq flag 'nonbasic-lower-bound)
			(not (column-has-l col)))
		   (decf v a))
		  ((and (eq flag 'nonbasic-upper-bound)
			(not (column-has-u col)))
		   (incf v a)))))
	(setf v (- v))
	(setf (aref values i) v
	      (aref header i) slack-col-ref)
	(cond ((and (column-has-l slack-col) (< v 0))
	       (splay-tree-fixnum-rational-set infeas i (* v v)))
	      ((and (column-has-u slack-col) (< 0 v))
	       (splay-tree-fixnum-rational-set infeas i (* v v))))))
    ;; return basis instance
    (let* ((refac-period (min 200 (+ 10 (floor m 20))))
	   (matrix (make-basis-matrix :lp lp :refactorization-period refac-period)))
      (fill-basis-matrix matrix lp header)
      (basis-matrix-lu-decomposition matrix)
      (%make-basis
       :matrix matrix
       :header header
       :obj-value (phase1-initial-objective-value lp flags)
       :dse-weights weights
       :reduced-costs rcosts 
       :column-flags flags
       :primal-values values
       :primal-infeas infeas))))
    


;;;;
(defun nonbasic-value (in-phase1 col flag)
  (if in-phase1
      (cond ((eq flag 'nonbasic-lower-bound)
	     (if (column-has-l col) 0 -1))
	    ((eq flag 'nonbasic-upper-bound)
	     (if (column-has-u col) 0 1))
	    (t
	     (error "inappropriate flag for phase 1")))
      (cond ((and (eq flag 'nonbasic-lower-bound)
		  (column-has-l col))
	     (column-l col))
	    ((and (eq flag 'nonbasic-upper-bound)
		  (column-has-u col))
	     (column-u col))
	    (t
	     (error "inappropriate flag for phase 2")))))


;;;;
(defun column-u-minus-l (in-phase1 col)
  (if in-phase1 
      (let ((range 2))
	(when (column-has-u col)
	  (decf range))
	(when (column-has-l col)
	  (decf range))
	range)
      (- (column-u col) (column-l col))))


;;;;
(defun column-l-minus-u (in-phase1 col)
  (if in-phase1 
      (let ((range (- 2)))
	(when (column-has-u col)
	  (incf range))
	(when (column-has-l col)
	  (incf range))
	range)
      (- (column-l col) (column-u col))))



;;;;  
(defun basis-update-primal-infeasability (primal-infeas i x col inphase1)
  (multiple-value-bind (sti there)
      (splay-tree-fixnum-rational-find-key primal-infeas i)
    (cond 
      ((and inphase1 
	    (or (and (column-has-l col) (< x 0))
		(and (column-has-u col) (< 0 x))))
       (if there
	   (setf (splay-tree-fixnum-rational-value primal-infeas sti) (* x x))
	   (splay-tree-fixnum-rational-set primal-infeas i (* x x))))
      ((and inphase1 there)
       (splay-tree-fixnum-rational-remove primal-infeas i))
      (inphase1)
      ((and (column-has-l col) (< x (column-l col)))
       (let ((delta (- (column-l col) x)))
	 (mulf delta delta)
	 (if there
	     (setf (splay-tree-fixnum-rational-value primal-infeas sti) delta)
	     (splay-tree-fixnum-rational-set primal-infeas i delta))))
      ((and (column-has-u col) (< (column-u col) x))
       (let ((delta (- x (column-u col))))
	 (mulf delta delta)
	 (if there
	     (setf (splay-tree-fixnum-rational-value primal-infeas sti) delta)
	     (splay-tree-fixnum-rational-set primal-infeas i delta))))
      (there 
       (splay-tree-fixnum-rational-remove primal-infeas i)))))
