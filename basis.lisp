;;;;; Basis info


(defstruct (basis
	     (:constructor %make-basis))
  (matrix         nil :type basis-matrix)
  (in-phase1      t   :type boolean)
  (status         'undef :type t)
  (header         #() :type vector)
  (dse-weights    #() :type vector)
  (reduced-costs  #() :type vector)
  (column-flags   #() :type vector)
  (obj-value      0   :type rational)
  (primal-values  #() :type vector))
    


;;;; Sets nonbasic variables for initial basis in phase 1
(defun phase1-initial-nonbasic-flags (lp flags)
  (dotimes (j (length (lp-columns lp)))
    (let ((col (aref (lp-columns lp) j)))
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
  (let ((n (length (lp-columns lp)))
	(z 0))
    (dotimes (j n z)
      (let ((col (aref (lp-columns lp) j))
	    (flag (aref flags j)))
	(cond ((zerop (column-c col)))
	      ((and (eq flag 'nonbasic-upper-bound)
		    (not (column-has-u col)))
	       (incf z (* (- (lp-obj-sense lp)) (column-c col))))
	      ((and (eq flag 'nonbasic-lower-bound)
		    (not (column-has-l col)))
	       (decf z (* (- (lp-obj-sense lp)) (column-c col)))))))))



;;;; Generates a dual feasible basis using the slack variables
(defun make-phase1-initial-basis (lp)
  (let* ((m (length (lp-active-row-refs lp)))
	 (n (length (lp-columns lp)))
	 (header (make-nvector m 0 fixnum))
	 (rcosts (make-nvector n 0 rational))
	 (flags  (make-nvector n 'unknown))
	 (values (make-nvector m 0 rational)))
    ;; set flags to best bound
    (phase1-initial-nonbasic-flags lp flags)
    ;; set phase1 reduced costs
    (dotimes (j n)
      (let ((flag (aref flags j)))
	(when (or (eq flag 'nonbasic-lower-bound)
		  (eq flag 'nonbasic-upper-bound))
	  (setf (aref rcosts j)
		(* (- (lp-obj-sense lp)) (column-c (aref (lp-columns lp) j)))))))
    ;; build basis header from slack variables
    ;; and set basic variable values
    (dotimes (i m)
      (let* ((v 0)
	     (row-ref (aref (lp-active-row-refs lp) i))
	     (row (aref (lp-rows lp) row-ref))
	     (slack-col-ref (row-slack-col-ref row)))
	(dotimes (k (length (row-col-refs row)))
	  (let* ((col-ref (aref (row-col-refs row) k))
		 (flag (aref flags col-ref))
		 (col (aref (lp-columns lp) col-ref))
		 (a (rational-in-column col (aref (row-col-indices row) k))))
	    (cond ((= col-ref slack-col-ref))
		  ((and (eq flag 'nonbasic-lower-bound)
			(not (column-has-l col)))
		   (decf v a))
		  ((and (eq flag 'nonbasic-upper-bound)
			(not (column-has-u col)))
		   (incf v a)))))
	(setf (aref values i) (- v)
	      (aref header i) slack-col-ref)))
    ;; return basis instance
    (let ((matrix (make-basis-matrix :lp lp)))
      (fill-basis-matrix matrix lp header)
      (basis-matrix-lu-decomposition matrix)
      (%make-basis
       :matrix matrix
       :header header
       :obj-value (phase1-initial-objective-value lp flags)
       :dse-weights (make-nvector m 1 rational)
       :reduced-costs rcosts 
       :column-flags flags
       :primal-values values))))
    


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



