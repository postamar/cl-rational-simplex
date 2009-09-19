


;;;; Step 1 of the simplex algorithm
;;;; Solve yB = c_B
(defun solve-yB=c_B (y)
  (dotimes (i (lp-data-m *lp*))
    (let ((c_Bi (* (lp-data-obj-sense *lp*)
		   (aref (lp-data-c *lp*) (aref *basis-header* i)))))
      (setf (aref *dense-column* i) c_Bi)))
  (let ((sparse-y (sparse-column-from-dense *dense-column*)))
    ;; go through eta file right to left
    (BTRAN sparse-y)
    ;; (check-yB=c_B y)
    (replace-sparse-column y sparse-y)
    (solver-log 'price "y = ~A~%" y)))




;;;; Step 2 of the simplex algorithm
;;;; Pricing - find j such that y . a_j < c_j and x_j = l_j, j \in N
;;;;                         or y . a_j > c_j and x_j = u_j, j \in N.
;;;; If no j exists, solution is optimal, return -1, else return j
(defun find-entering-column (y)
  (dotimes (j (lp-data-n *lp*) -1)
    (let ((x_j (aref *x* j))
	  (l_j (aref (lp-data-l *lp*) j))
	  (u_j (aref (lp-data-u *lp*) j))
	  (l_j-p (bit (lp-data-l-p *lp*) j))
	  (u_j-p (bit (lp-data-u-p *lp*) j))
	  (rhs (* (lp-data-obj-sense *lp*)
		  (aref *c-numerators* j)
		  (sparse-column-denominator y))))
      ;; compute dot product y . a_j
      (let ((num-dot-product 0)
	    (A_.j (aref *a-sparse-columns* j))
	    (si-y   -1) (si-A_.j -1) 
	    (di-y   -1) (di-A_.j -1) 
	    (y_i     0) (A_ij     0))
	(iterate-through-sparse-columns
	 (y    A_.j
	  si-y si-A_.j
	  di-y di-A_.j
	  y_i  A_ij)
	 (incf num-dot-product
	       (* y_i A_ij)))
	;; compare to cost
	(when (or (and (not (zerop l_j-p))
		       (= x_j l_j)
		       (< num-dot-product rhs))
		  (and (not (zerop u_j-p))
		       (= x_j u_j)
		       (> num-dot-product rhs)))
	  (solver-log 'price 
		      "entering column j = ~A~%num(y).num(A_.j) = ~A, c_j . sense . denom(y) = ~A~%" 
		      j num-dot-product rhs)
	  (return-from find-entering-column j))))))
	  


;;;; Step 3 of the simplex algorithm
;;;; Solve Bd = a_j
(defun solve-Bd=a_j (d j)
  ;; set d = a_j
  (replace-sparse-column d (deep-copy-sparse-column (aref *a-sparse-columns* j)))
  ;; go through eta file left to right
  (FTRAN d)
 ; (check-Bd=a_j d j)
  (solver-log 'pivot "d = ~A~%" d))
		     
	  

;;;; Step 4 of the simplex algorithm
;;;; Find largest step and the index of the leaving column
;;;;     with x_j(step) = x_j + step and x_B(step) = x_B - step d when x_j = l_j
;;;;  or with x_j(step) = x_j - step and x_B(step) = x_B + step d when x_j = u_j.
;;;; Problem is unbounded if l_j <= x_j(step) <= u_j 
;;;;                     and l_B <= x_B(step) <= u_B hold for all positive step.



;;;; Returns acceptable bound on step for column at header-index in the basis, or -1
(defun step-bound (d d-factor sparse-index header-index)
  (let* ((Bi (aref *basis-header* header-index))
	 (x_Bi (aref *x* Bi))
	 (l_Bi (aref (lp-data-l *lp*) Bi))
	 (u_Bi (aref (lp-data-u *lp*) Bi))
	 (l_Bi-p (bit (lp-data-l-p *lp*) Bi)) 
	 (u_Bi-p (bit (lp-data-u-p *lp*) Bi))
	 (num-d_Bi (* d-factor 
		      (aref (sparse-column-numerators d) sparse-index))))
    (cond ((and (not (zerop u_Bi-p))
		(> num-d_Bi 0))
	   (/ (* (- u_Bi x_Bi)
		 (sparse-column-denominator d)) 
	      num-d_Bi))
	  ((and (not (zerop l_Bi-p))
		(< num-d_Bi 0))
	   (/ (* (- l_Bi x_Bi) 
		 (sparse-column-denominator d)) 
	      num-d_Bi))
	  (t
	   -1))))



;;;; Maximizes step and determines the leaving column accordingly.
(defun pivot (d j)
  (let ((step 0)
	(step-p nil)
	(exiting-header-index -1)
	(l_j (aref (lp-data-l *lp*) j))
	(u_j (aref (lp-data-u *lp*) j))
	(l_j-p (bit (lp-data-l-p *lp*) j))
	(u_j-p (bit (lp-data-u-p *lp*) j)))
    ;; check if step is bounded by l_j or u_j
    (when (and (not (zerop l_j-p))
	       (not (zerop u_j-p)))
      (setf step (- u_j l_j)
	    step-p t))
    ;; check if step is bounded by l_B and u_B
    (let ((d-nonzeros (length (sparse-column-indices d)))
	  (d-factor (if (and (not (zerop l_j-p)) 
			     (= (aref *x* j) l_j)) 
			-1 
			1)))
      (dotimes (k d-nonzeros (values exiting-header-index step))
	(let* ((header-index (aref (sparse-column-indices d) k))
	       (bound (step-bound d d-factor k header-index)))
	  (when (and (/= -1 bound)
		     (or (not step-p)
			 (< bound step)))
	    (setf step bound
		  step-p t
		  exiting-header-index header-index)))))))
	    


;;;; Step 5 of the simplex algorithm
;;;; Updates solution component values
;;;; x_j +=/-= step
;;;; x_B +=/-= step * d
(defun update-solution (j step d)
  (let ((l_j-p (bit (lp-data-l-p *lp*) j))
	(l_j (aref (lp-data-l *lp*) j))
	(x_j (aref *x* j))
	(d-nonzeros (length (sparse-column-indices d)))
	(d-indices (sparse-column-indices d))
	(d-numerators (sparse-column-numerators d)))
    (if (and (not (zerop l_j-p)) 
	     (= x_j l_j))
	(incf (aref *x* j) step)
	(decf (aref *x* j) step))
    (dotimes (k d-nonzeros)
      (let ((i (aref d-indices k))
	    (num-d_i (aref d-numerators k)))
	(incf (aref *x* (aref *basis-header* i))
	      (/ (* (if (and (not (zerop l_j-p))
			     (= x_j l_j))
			(- step) 
			step) 
		    num-d_i)
		 (sparse-column-denominator d)))))))



;;;; Step 5 bis of the simplex algorithm
;;;; Updates basis information
(defun update-basis (j exiting-header-index d)
  ;; update the basis header
  (setf (aref *basis-header* exiting-header-index) j)
  (solver-log 'sol "new basis header: ~A~%" *basis-header*)
  ;; update eta file
  (vector-push-extend d *eta-column-vector-file*)
  (vector-push-extend exiting-header-index *eta-column-index-file*))
      
    
    
;;;; Performs the simplex algorithm
(defun simplex ()
  (let ((skip-step-1 nil)
	(y (make-sparse-column))
	(d (make-sparse-column)))
    (solver-log 'gen "~%Beginning simplex iterations~%")
    (solver-log 'sol "basis header: ~A~%" *basis-header*)
    (dotimes (simplex-iteration *iteration-cap*)
      (solver-log 'iter "~%----- Iteration ~A:~%" simplex-iteration)
      ;; check objective limit
      (let ((z (current-value)))
	(solver-log 'obj "z = ~A~%" (float z))
	(solver-log 'sol "x = ~A~%" *x*)
	(when (and *obj-limit*
		   (>= (* (lp-data-obj-sense *lp*) z) 
		       (* (lp-data-obj-sense *lp*) *obj-limit*)))
	  (solver-log 'gen "Objective limit reached.~%")
	  (setf *solver-status* 'feasible)
	  (return-from simplex *solver-status*)))
      ;; check refactorization limit
      (unless (< (fill-pointer *eta-column-index-file*) *refactorization-limit*)
	(solver-log 'gen "Refactoring basis~%")
	(refactorize-basis))
      ;; step 1
      (if skip-step-1
	  (setf skip-step-1 nil)
	  (solve-yB=c_B y))
      ;; step 2
      (let ((j (find-entering-column y)))
	;; check for optimality
	(when (= -1 j)
	  (setf *solver-status* 'optimal)
	  (return-from simplex *solver-status*))
	;; step 3
	(solve-Bd=a_j d j)
	(fill-dense-column-from-sparse *dense-column* d)
	;; step 4
	(multiple-value-bind (exiting-header-index step) (pivot *dense-column* j)
	  ;; check if problem is unbounded
	  (when (= -1 step)
	    (setf *solver-status* 'unbounded)
	    (return-from simplex *solver-status*))
	  (solver-log 'pivot "step = ~A~%" step)
	  (unless (= -1 exiting-header-index)
	    (solver-log 'pivot "exiting column i = ~A~%" 
			(aref *basis-header* exiting-header-index)))
	  ;; step 5
	  (update-solution j step *dense-column*)
	  (if (= -1 exiting-header-index)
	      (setf skip-step-1 t) ; skip step 1 of next simplex iteration
	      (update-basis j exiting-header-index d)))))))

      
	 


