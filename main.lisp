(in-package :rationalsimplex)

;;;;; Front-end of the rational simplex solver:
;;;;; * create new LP instance with lp-make-new
;;;;; * create new LP instance from file with lp-load-from-mps
;;;;; * add new column with lp-add-variable
;;;;; * add new row with lp-add-constraint
;;;;; * solve with lp-solve
;;;;; * profile the solver with rational-simplex-profiling
;;;;; * get solution values with lp-primal-solution,
;;;;;   lp-dual-solution, lp-slack-values, lp-reduced-costs.


;;;; 
(defparameter *lp-name-counter* 0
  "Counter for generating unique LP instance names.")


;;;; 
(defun lp-make-new (&key (name "") (obj-name "OBJ") (obj-sense -1))
  "Creates new, empty LP instance, arguments are optional:
  :NAME is the linear program name,
  :OBJ-NAME is the objective function name,
  :OBJ-SENSE is the optimization sense,
    -1 by default for minimization, else set to 1 for maximization.
Returns LP object instance."
  (declare (string name obj-name)
	   (real obj-sense))
  (make-lp :name (if (string= name "") 
		     (concatenate 'string
				  "UNNAMED-LP-" 
				  (princ-to-string (incf *lp-name-counter*)))
		     name)
	   :obj-name obj-name
	   :obj-sense obj-sense
	   :active-row-refs (make-empty-adjvector-fixnum)
	   :active-col-refs (make-empty-adjvector-fixnum)
	   :active-row-inds (make-empty-adjvector-fixnum)
	   :active-col-inds (make-empty-adjvector-fixnum)
	   :columns (make-adjvector-column (make-column :hsv (make-hsv)) 0)
	   :rows (make-adjvector-row (make-row) 0)
	   :col-ref-by-name (make-hash-table :test 'equal)
	   :row-ref-by-name (make-hash-table :test 'equal)))


;;;; Auxilliary function which parses constraint matrix data (row or col)
;;;; Returns ref-val associative list and success flag
;;;; ref is row/col reference number and val is a rational 
(defun parse-coef-seq (lp coefs err-string ref-hash-table ref-ubd)
  (declare (lp lp) 
	   (sequence coefs) 
	   (integer ref-ubd))
  (macrolet 
      ;; graceful exit macro
      ((errquit (format-str &rest format-args)
	 `(progn (format t (concatenate 'string 
					"In LP ~A,"
					err-string
					,format-str ".~%")
			 (lp-name lp) ,@format-args)
		 (return-from parse-coef-seq (values nil nil)))))
    (let ((coef-alist))
      ;; parse the assoc. list
      (map nil
	   #'(lambda (pair)
	       (unless (typep pair 'cons)
		 (errquit "~A in sequence :COEF is not of type CONS"
			  pair))
	       (destructuring-bind (id . coef-val) pair
		 (cond 
		   ((not (typep coef-val 'real))
		    (errquit
		     "coef value in ~A in sequence :COEF is not of type REAL"
		     pair))
		   ((typep id 'string)
		    (multiple-value-bind (coef-ref present-p)
			(gethash id ref-hash-table)
		      (unless present-p
			(errquit 
			 "~A in sequence :COEF is not a valid name"
			 id))
		      (unless (zerop coef-val)
			(push (cons coef-ref (rationalize coef-val)) coef-alist))))
		   ((typep id 'integer)
		    (unless (<= 0 id ref-ubd)
		      (errquit
		       "ref ~A in sequence :COEF is out of bounds (0,~A)"
		       id
		       ref-ubd))
		    (unless (zerop coef-val)
		      (push (cons id (rationalize coef-val)) coef-alist)))
		   (t
		    (errquit
		     "~A is not a valid row ref number or row name"
		     id)))))
	   coefs)
      ;; return result
      (values coef-alist t))))



;;;; 
(defun lp-add-variable (lp &key (name "") (cost 0) (coefs) (lbound 0) (ubound))
  "Adds column to LP instance, all other arguments are optional:
  :NAME is the column name,
  :COST is the objective cost of the variable,
  :COEF is an associative list of (REF . VALUE) pairs,
    with REF the name or ref-number of the contraints, 
    and VALUE the constraint coefficient of the variable.
  :LBOUND and :UBOUND, if set, are its lower and upper bounds,
    by default the variable is non-negative.
Returns T on success, NIL on failure."
  (declare (lp lp)
	   (string name)
	   (real cost)
	   (sequence coefs))
  (macrolet 
      ;; graceful exit macro
      ((errquit (format-str &rest format-args)
	 `(progn (format t (concatenate 'string 
					"In LP ~A, failed to add variable: " 
					,format-str ".~%")
			 (lp-name lp) ,@format-args)
		 (return-from lp-add-variable nil))))
    (let* ((ref (adjvector-column-fill-pointer (lp-columns lp))) ; column is appended
	   (hsv (make-hsv))
	   (name (if (string= name "") 
		     (concatenate 'string "UNNAMED-VAR-" (princ-to-string ref))
		     name)))
      ;; check bounds
      (unless (or (not lbound) (typep lbound 'real))
	(errquit "lbound is not NIL and not of type REAL"))
      (unless (or (not ubound) (typep ubound 'real))
	(errquit "ubound is not NIL and not of type REAL"))
      (when (and lbound ubound (> lbound ubound))
	(errquit "lbound (~A) > ubound (~A)" lbound ubound))
      ;; check for name conflicts
      (when (gethash name (lp-col-ref-by-name lp))
	(errquit "name ~A already used" name))
      ;; parse coef sequence 
      (multiple-value-bind (hsv-alist parse-p)
	  (parse-coef-seq lp coefs " failed to add variable: " 
			  (lp-row-ref-by-name lp)
			  (1- (adjvector-row-fill-pointer (lp-rows lp))))
	(unless parse-p
	  (return-from lp-add-variable nil))
	;; check coef alist refs and build hsv
	(setf hsv-alist (sort hsv-alist #'< :key #'car))
	(let ((last -1)
	      (cdenom 1))
	  (declare (integer cdenom))
	  (dolist (pair hsv-alist)
	    (destructuring-bind (ref . val) pair
	      (when (= last ref)
		(errquit "duplicate refs in sequence :COEF"))
	      ;; compute common denominator
	      (setf cdenom (/ (denominator val) (gcd cdenom (denominator val)))
		    last ref)))
	  ;; build hyper-sparse vector 
	  (dolist (pair hsv-alist)
	    (destructuring-bind (ref . val) pair
	      (assert (= 1 (denominator (* cdenom val))))
	      (hsv-add ref (* cdenom val) hsv)))
	  (setf (hsv-coef hsv) (/ 1 cdenom))
	  ;; simplify its fractions 
	  (hsv-normalize hsv))
	;; add column into LP data structure instance
	(setf (gethash name (lp-col-ref-by-name lp)) ref)
	(adjvector-fixnum-push-extend ref (lp-active-col-refs lp))
	(adjvector-fixnum-push-extend ref (lp-active-col-inds lp))
	(adjvector-column-push-extend
	 (make-column
	  :name name
	  :ref ref
	  :is-active t
	  :c (rationalize cost)
	  :has-l (when lbound t)
	  :has-u (when ubound t)
	  :l (if lbound lbound 0)
	  :u (if ubound ubound 0)
	  :hsv hsv)
	 (lp-columns lp))
	;; update row information
	(let ((col-index -1))
	  (dolist (pair hsv-alist)
	    (let ((row (adjvector-row-ref (lp-rows lp) (car pair))))
	      (adjvector-fixnum-push-extend ref (row-col-refs row))
	      (adjvector-fixnum-push-extend (incf col-index) (row-col-indices row))))))
      ;; return T on success
      t)))

    

;;;; 
(defun lp-add-constraint (lp &key (name "") (lhs) (rhs 0) (coefs)) 
  "Adds row to LP instance, all other arguments are optional:
  :NAME is the row name,
  :COEF is an associative list of (REF . VALUE) pairs,
    with REF the name or ref-number of the variable, 
    and VALUE the constraint coefficient of the variable,
  :LHS and :RHS are the lower and upper bounds on the row activity,
    i.e.  :LHS <= (x^T . :COEFS) <= :RHS.
Returns T on success, NIL on failure."
  (declare (lp lp)
	   (string name)
	   (sequence coefs))
  (macrolet 
      ;; graceful exit macro
      ((errquit (format-str &rest format-args)
	 `(progn (format t (concatenate 'string 
					"In lp ~A, failed to add constraint: " 
					,format-str ".~%")
			 (lp-name lp) ,@format-args)
		 (return-from lp-add-constraint nil))))
    (let* ((ref (adjvector-row-fill-pointer (lp-rows lp))) ; row is appended
	   (name (if (string= name "") 
		     (concatenate 'string "UNNAMED-CST-" (princ-to-string ref))
		     name))
	   (slack-ref (adjvector-column-fill-pointer (lp-columns lp))) ;slack column too
	   (slack-hsv (make-hsv))
	   (slack-name (concatenate 'string "SLACK-" name))
	   (slack-col (make-column
		       :name slack-name
		       :ref slack-ref
		       :is-active t
		       :is-slack t
		       :has-l (when rhs t)
		       :has-u (when lhs t)
		       :l (if rhs (- rhs) 0)
		       :u (if lhs (- lhs) 0)
		       :hsv (prog1 slack-hsv (hsv-add ref 1 slack-hsv))))
	   (row (make-row 
		 :name name
		 :ref ref
		 :is-active t
		 :slack-col-ref slack-ref
		 :col-refs (make-empty-adjvector-fixnum)
		 :col-indices (make-empty-adjvector-fixnum))))
      ;; check bounds
      (unless (or (not lhs) (typep lhs 'real))
	(errquit "lhs is not NIL and not of type REAL"))
      (unless (or (not rhs) (typep rhs 'real))
	(errquit "rhs is not NIL and not of type REAL"))
      (when (and lhs rhs (> lhs rhs))
	(errquit "lhs (~A)  >  rhs (~A)" lhs rhs))
      ;; check for name conflicts
      (when (gethash name (lp-row-ref-by-name lp))
	(errquit "constraint name ~A already used" name))
      (when (gethash slack-name (lp-col-ref-by-name lp))
	(errquit "variable name ~A already used" slack-name))
      ;; parse coef sequence 
      (multiple-value-bind (coef-alist parse-p)
	  (parse-coef-seq lp coefs " failed to add constraint: " 
			  (lp-col-ref-by-name lp)
			  (1- slack-ref))
	(unless parse-p
	  (return-from lp-add-constraint nil))
	;; check coef alist refs 
	(setf coef-alist (sort coef-alist #'< :key #'car))
	(let ((last -1))
	  (dolist (pair coef-alist)
	    (let ((col-ref (car pair)))
	      (when (= last col-ref)
		(errquit "duplicate refs in sequence :COEF"))
	      (setf last col-ref))))
	;; update columns
	(dolist (pair coef-alist)
	  (destructuring-bind (col-ref . val) pair
	    (let* ((hsv (column-hsv (adjvector-column-ref (lp-columns lp) col-ref)))
		   (col-index (hsv-length hsv)))
	      (dotimes (k col-index)
		(setf (aref (hsv-vis hsv) k) (* (denominator val) 
						(numerator (hsv-coef hsv))
						(aref (hsv-vis hsv) k))))
	      (setf (hsv-coef hsv) (/ 1 (* (denominator (hsv-coef hsv)) 
					   (denominator val))))
	      (hsv-add ref (numerator val) hsv)
	      (hsv-normalize hsv)
	      (adjvector-fixnum-push-extend col-ref (row-col-refs row))
	      (adjvector-fixnum-push-extend col-index (row-col-indices row)))))
	;; add slack column
	(setf (gethash slack-name (lp-col-ref-by-name lp)) slack-ref)
	(adjvector-fixnum-push-extend slack-ref (lp-active-col-refs lp))
	(adjvector-fixnum-push-extend slack-ref (lp-active-col-inds lp))
	(adjvector-column-push-extend slack-col (lp-columns lp))
	;; add row
	(adjvector-fixnum-push-extend slack-ref (row-col-refs row))
	(adjvector-fixnum-push-extend 0 (row-col-indices row))
	(setf (gethash name (lp-row-ref-by-name lp)) ref)
	(adjvector-fixnum-push-extend ref (lp-active-row-refs lp))
	(adjvector-fixnum-push-extend ref (lp-active-row-inds lp))
	(adjvector-row-push-extend row (lp-rows lp))
	;; return T on success
	t))))



;;;;
(defun lp-load-from-mps (full-mps-file-path)
  "Creates LP instance from data in MPS format,
  argument is a string containing a full path to the .mps file.
Returns LP object instance on success, NIL on failure."
  (declare (string full-mps-file-path))
  (let ((mps (load-from-mps full-mps-file-path)))
    (if mps
	(mps->lp mps)
	(progn
	  (format t "Failed loading MPS data from file ~A.~%" full-mps-file-path)
	  ;; return NIL on failure
	  nil))))




;;;;
(defun lp-solve (lp &key 
		 (z-print-freq 1) 
		 (min-z) 
		 (max-z)
		 (max-total-time)
		 (max-phase-time)
		 (max-total-iters)
		 (max-phase-iters))
  "Solves the linear program provided as argument using the dual-simplex method:
  :Z-PRINT-FREQ sets how often the current objective value is displayed,
  :MIN-Z sets a lower limit on the objective value,
  :MAX-Z sets an upper limit,
  :MAX-TOTAL-TIME sets a limit on total solve time,
  :MAX-PHASE-TIME sets a limit on phase1 and phase2 solve time,
  :MAX-TOTAL-ITERS sets a limit on total solve iteration count,
  :MAX-PHASE-ITERS sets a limit on phase1 and phase2 iteration count."
  (declare (lp lp)
	   (integer z-print-freq))
  ;; type checking
  (macrolet 
      ((errmsg (keyword)
	 (let ((string (concatenate 'string 
				    "Solve error: :"
				    (princ-to-string keyword)
				    " must be of type REAL~%")))
	   `(when (and ,keyword (not (typep ,keyword 'real)))
	      (format t ,string)))))
    (errmsg min-z)
    (errmsg max-z)
    (errmsg max-total-time)
    (errmsg max-phase-time)
    (errmsg max-total-iters)
    (errmsg max-phase-iters)
    ;; reset random state
    (setf *random-state* (make-random-state *simplex-random-state*))
    ;; create basis and simplex data structure instances
    (let ((sd nil)
	  (basis (make-phase1-initial-basis lp)))
      ;; check basis was created succesfully
      (unless basis
	(format t "Solve error: could not create initial basis for LP ~A.~%" (lp-name lp))
	(return-from lp-solve nil))
      ;; check simplex-data was created succesfully
      (unless (setf sd (make-simplex lp basis))
	(format t "Solve error: could not initialize the dual-simplex method for LP ~A.~%" 
		(lp-name lp))
	(return-from lp-solve nil))
      ;; solve the lp
      (let ((status (dual-simplex sd
				  :min-z min-z
				  :max-z max-z
				  :z-print-freq z-print-freq
				  :max-total-time max-total-time
				  :max-phase-time max-phase-time
				  :max-total-iters max-total-iters
				  :max-phase-iters max-phase-iters)))
	(format t "~%")
	(values status (simplex-basis sd) (simplex-stats sd))))))



;;;;
(defun lp-primal-solution (lp basis &optional (by-name nil))
  "Computes the primal solution values for a given basis.
Returns an associative list of (column-ref . val) pairs, 
or (name . val) pairs if by-name is not NIL."
  (declare (lp lp)
	   (basis basis))
  (macrolet
      ((exit (format-str &rest format-args)
	 `(progn
	    (format t (concatenate 'string 
				   "Error computing primal solution: basis is not dual-feasible, " 
				   ,format-str ".~%") 
		    ,@format-args)
	    (return-from lp-primal-solution nil))))
    (when (basis-in-phase1 basis)
      (exit "(phase 1 basis)"))
    (let* ((n (length (basis-column-flags basis)))
	   (m (length (basis-header basis)))
	   (x nil))
      (flet ((add-to-list (col val)
	       (unless (column-is-slack col)
		 (push (cons (if by-name (column-name col) (column-ref col)) val) x))))
	;; get basic variable values
	(dotimes (i m)
	  (let ((j (aref (basis-header basis) i)))
	    (unless (eq (aref (basis-column-flags basis) j) 'basic)
	      (exit "variable #~D is in header yet not set to 'BASIC" j))
	    (add-to-list (adjvector-column-ref (lp-columns lp) j) 
			 (aref (basis-primal-values basis) i))))
	;; get non-basic variable values
	(dotimes (j n x)		; return assoc list when done
	  (let ((flag (aref (basis-column-flags basis) j))
		(col (adjvector-column-ref (lp-columns lp) j)))
	    (cond ((eq 'basic flag)) 
		  ((eq 'nonbasic-upper-bound flag)
		   (if (column-has-u col)
		       (add-to-list col (column-u col))
		       (exit "variable #~D has no upper bound, yet is set to 'NONBASIC-UPPER-BOUND" j)))
		  ((eq 'nonbasic-lower-bound flag)
		   (if (column-has-l col)
		       (add-to-list col (column-l col))
		       (exit "variable #~D has no lower bound, yet is set to 'NONBASIC-LOWER-BOUND" j)))
		  (t
		   (exit "variable #~D is set to ~A" j flag)))))))))



;;;;
(defun lp-slack-values (lp basis &optional (by-name nil))
  "Computes the primal slack values in each constraint for a given basis.
Returns an associative list of (row-ref . val) pairs, 
or (name . val) pairs if by-name is not NIL."
  (declare (lp lp)
	   (basis basis))
  (macrolet
      ((exit (format-str &rest format-args)
	 `(progn
	    (format t (concatenate 'string 
				   "Error computing slack values: basis is not dual-feasible, " 
				   ,format-str ".~%") 
		    ,@format-args)
	    (return-from lp-slack-values nil))))
    (when (basis-in-phase1 basis)
      (exit "(phase 1 basis)"))
    (let ((slack nil)
	  (nrows (adjvector-row-fill-pointer (lp-rows lp))))
      (flet ((add-to-list (row val)
	       (push (cons (if by-name (row-name row) (row-ref row)) (- val)) slack)))
	;; go through each row and get slack values
	(dotimes (row-ref nrows slack)	; return assoc list when done
	  (let* ((row (adjvector-row-ref (lp-rows lp) row-ref))
		 (slack-col-ref (row-slack-col-ref row))
		 (slack-col (adjvector-column-ref (lp-columns lp) slack-col-ref))
		 (flag (aref (basis-column-flags basis) slack-col-ref)))
	    (cond ((eq 'basic flag)
		   (let ((val 0) (val-found nil))
		     (dotimes (i (length (basis-header basis)))
		       (when (= (aref (basis-header basis) i) slack-col-ref)
			 (setf val (aref (basis-primal-values basis) i)
			       val-found t)))
		     (unless val-found
		       (exit "slack variable of row #~D is 'BASIC, yet cannot be found in basis header" row-ref))
		     (add-to-list row val)))
		  ((eq 'nonbasic-upper-bound flag)
		   (if (column-has-u slack-col)
		       (add-to-list row (column-u slack-col))
		       (exit "slack variable of row #~D has no upper bound, yet is set to 'NONBASIC-UPPER-BOUND" row-ref)))
		  ((eq 'nonbasic-lower-bound flag)
		   (if (column-has-l slack-col)
		       (add-to-list row (column-l slack-col))
		       (exit "slack variable of row #~D has no lower bound, yet is set to 'NONBASIC-LOWER-BOUND" row-ref)))
		  (t
		   (exit "slack variable of row #~D is set to ~A" row-ref flag)))))))))
      


;;;;
(defun lp-dual-solution (lp basis &optional (by-name nil))
  "Computes the full dual solution vector for a given basis
Returns an associative list of (row-ref . val) pairs, 
or (name . val) pairs if by-name is not NIL."
  (declare (lp lp)
	   (basis basis))
  (macrolet
      ((exit (format-str &rest format-args)
	 `(progn
	    (format t (concatenate 'string 
				   "Error computing dual solution: basis is not dual-feasible, " 
				   ,format-str ".~%") 
		    ,@format-args)
	    (return-from lp-dual-solution nil))))
    (when (basis-in-phase1 basis)
      (exit "(phase 1 basis)"))
    (let* ((rhs (make-hsv))
	   (tr (make-tran (basis-matrix basis)))
	   (bh (basis-header basis))
	   (m (length bh))
	   (cdenom 1)
	   (vals (make-array (adjvector-row-fill-pointer (lp-rows lp)) :element-type 'rational)))
      ;; get common denominator
      (dotimes (i m)
	(let* ((col (adjvector-column-ref (lp-columns lp) (aref bh i)))
	       (c (column-c col))
	       (cd (denominator c)))
	  (unless (zerop c)
	    (mulf cdenom (/ cd (gcd cd cdenom))))))
      ;; create right-hand-side for BTRAN
      (setf (hsv-coef rhs) (/ 1 cdenom))
      (dotimes (i m)
	(let* ((col (adjvector-column-ref (lp-columns lp) (aref bh i)))
	       (c (column-c col)))
	  (unless (zerop c)
	    (hsv-add i (* (numerator c) (/ cdenom (denominator c))) rhs))))
      ;; compute and store multipliers
      (btran tr rhs)
      (check-btran basis lp tr rhs)
      (dotimes (k (hsv-length (tran-hsv tr)))
	(setf (aref vals (adjvector-fixnum-ref (lp-active-row-refs lp) 
					       (aref (hsv-is (tran-hsv tr)) k))) 
	      (* (hsv-coef (tran-hsv tr)) (aref (hsv-vis (tran-hsv tr)) k))))
      ;; return associative list
      (let ((y nil))
	(dotimes (row-ref (length vals) y)
	  (let ((row (adjvector-row-ref (lp-rows lp) row-ref)))
	    (push (cons (if by-name (row-name row) row-ref) (aref vals row-ref)) y)))))))
  
  

;;;;
(defun lp-reduced-costs (lp basis &optional (by-name nil))
  "Computes the reduced costs for a given basis.
Returns an associative list of (column-ref . val) pairs, 
or (name . val) pairs if by-name is not NIL."
  (declare (lp lp)
	   (basis basis))
  (macrolet
      ((exit (format-str &rest format-args)
	 `(progn
	    (format t (concatenate 'string "Error computing reduced costs: basis is not dual-feasible, " 
				   ,format-str ".~%") 
		    ,@format-args)
	    (return-from lp-reduced-costs nil))))
    (when (basis-in-phase1 basis)
      (exit "(phase 1 basis)"))
    (let* ((n (length (basis-column-flags basis)))
	   (x nil))
      (flet ((add-to-list (col val)
	       (unless (column-is-slack col)
		 (push (cons (if by-name (column-name col) (column-ref col)) val) x))))
	(dotimes (j n x) ; return assoc list when done
	  (add-to-list (adjvector-column-ref (lp-columns lp) j) 
		       (aref (basis-reduced-costs basis) j)))))))




;;;; Switches profiling ON for important functions of the solver
(defun rational-simplex-profiling ()
  (profile hsv-normalize)
  (profile btran-solve-eta)
  (profile btran-multiply-eta)
  (profile ftran-solve-eta)
  (profile ftran-multiply-eta)
  (profile simplex-iteration)
  (profile btran)
  (profile ftran)
  (profile choose-exiting-basis-index)
  (profile choose-entering-basis-index)
  (profile compute-pivot-row)
  (profile simplex-basis-update)
  (profile simplex-basis-update-primal)
  (profile simplex-basis-matrix-update)
  (profile lu-update)
  (profile basis-matrix-lu-factorization))
