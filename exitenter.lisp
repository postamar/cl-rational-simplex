(in-package :rationalsimplex)

;;;;; Dual simplex functions for choosing exiting and entering variables
;;;;;




;;;; Chooses exiting variable
(defun choose-exiting-basis-index (sd)
  (let* ((lp (simplex-lp sd))
	 (b (simplex-basis sd))
	 (infeas (basis-primal-infeas b))
	 (best-k '())
	 (infty nil)
	 (best-val 0))
    (map-splay-tree-fixnum-rational
     #'(lambda (k delta-square)
	 (let ((weight (aref (basis-dse-weight-vis b) k)))
	   (assert (< 0 (basis-dse-coef b)))
	   (cond ((< weight 0)
		  (error "negative dse weight"))
		 ((< delta-square 0)
		  (error "negative delta-square"))
		 ((zerop delta-square)
		  (error "zero delta in primal infeasabilities"))
		 ((and (zerop weight) infty)
		  (push k best-k))
		 ((zerop weight)
		  (setf best-k (list k)
			infty t))
		 (infty)
		 (t
		  (let ((val (/ delta-square (* (basis-dse-coef b) weight))))
		    (cond
		      ((< best-val val)
		       (setf best-val val
			     best-k (list k)))
		      ((= best-val val)
		       (push k best-k))))))))
       infeas)
    (if (eq '() best-k)
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
  (declare (optimize (speed 1) (safety 0) (debug 0)))
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
	       (declare (fixnum root end))
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
	   (incf (the integer (aref (simplex-pivot-row-values sd)
				    (simplex-pivot-row-length sd)))
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
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (let* ((b (simplex-basis sd))
	 (lp (simplex-lp sd))
	 (best-val 0)
	 (best-k -1)
	 (delta-obj 0)
	 (boxed nil)
	 (slope (abs (simplex-delta sd)))
	 (sign (signum (simplex-delta sd))))
    (declare (rational slope delta-obj))
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
			      (* (the rational (column-u-minus-l (basis-in-phase1 b) colq))
				 (abs (aref (simplex-pivot-row-values sd) best-k))))))
	   (when (< next-slope 0)
	     (return))
	   (unless (zerop (the rational (column-u-minus-l (basis-in-phase1 b) colq)))
	     (setf (aref (simplex-flip-col-refs sd) (simplex-n-flips sd)) q)
	     (incf (simplex-n-flips sd))
	     (incf delta-obj
		   (* (column-c colq)
		      (the rational
			(if (eq (aref (basis-column-flags b) q) 'nonbasic-lower-bound)
			    (column-u-minus-l (basis-in-phase1 b) colq)
			    (column-l-minus-u (basis-in-phase1 b) colq)))))
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
  





