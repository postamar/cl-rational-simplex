;;;; LP preprocessing
(defparameter *presolve-with-bounds* nil)


;;;; Holds computed bounds during presolve
(defstruct presolve-bounds
  m
  n
  row-has-l
  row-has-u
  dual-row-has-l
  dual-row-has-u
  col-has-l
  col-has-u
  dual-col-has-l
  dual-col-has-u
  row-l
  row-u
  dual-row-l
  dual-row-u
  col-l
  col-u
  dual-col-l
  dual-col-u)



;;;; Constructor
(defun init-bounds (m n)
  (make-presolve-bounds
   :m m
   :n n
   :row-has-l (make-array m :initial-element nil :element-type 'boolean)
   :row-has-u (make-array m :initial-element nil :element-type 'boolean)
   :dual-row-has-l (make-array n :initial-element nil :element-type 'boolean)
   :dual-row-has-u (make-array n :initial-element nil :element-type 'boolean)
   :col-has-l (make-array n :initial-element nil :element-type 'boolean)
   :col-has-u (make-array n :initial-element nil :element-type 'boolean)
   :dual-col-has-l (make-array m :initial-element nil :element-type 'boolean)
   :dual-col-has-u (make-array m :initial-element nil :element-type 'boolean)
   :row-l (make-array m :initial-element 0 :element-type 'rational)
   :row-u (make-array m :initial-element 0 :element-type 'rational)
   :dual-row-l (make-array n :initial-element 0 :element-type 'rational)
   :dual-row-u (make-array n :initial-element 0 :element-type 'rational)
   :col-l (make-array n :initial-element 0 :element-type 'rational)
   :col-u (make-array n :initial-element 0 :element-type 'rational)
   :dual-col-l (make-array m :initial-element 0 :element-type 'rational)
   :dual-col-u (make-array m :initial-element 0 :element-type 'rational)))



;;;;
(defun compute-row-bounds (lp bounds)
  (dotimes (i (presolve-bounds-m bounds))
    (let ((row_i (aref (lp-rows lp) i)))
      (unless (row-is-removed row_i)
	(let ((has-rl t)
	      (has-ru t)
	      (rl 0)
	      (ru 0))
	  (dorow (lp i :visit-slack nil :j j :index index :col col_j)
	    (when (not (or has-rl has-ru))
	      (return))
	    (let ((a_ij (aref (column-values col_j) index)))
	      (when (< a_ij 0)
		(if (column-has-u col_j)
		    (incf rl (* a_ij (column-u col_j)))
		    (setf has-rl nil))
		(if (column-has-l col_j)
		    (incf ru (* a_ij (column-l col_j)))
		    (setf has-ru nil)))
	      (when (> a_ij 0)
		(if (column-has-u col_j)
		    (incf ru (* a_ij (column-u col_j)))
		    (setf has-ru nil))
		(if (column-has-l col_j)
		    (incf rl (* a_ij (column-l col_j)))
		    (setf has-rl nil)))))
	  (setf (aref (presolve-bounds-row-has-l bounds) i) has-rl
		(aref (presolve-bounds-row-has-u bounds) i) has-ru
		(aref (presolve-bounds-row-l bounds) i) rl
		(aref (presolve-bounds-row-u bounds) i) ru))))))



;;;;
(defun compute-column-bounds (lp bounds)
  (dotimes (i (presolve-bounds-m bounds))
    (let ((row_i (aref (lp-rows lp) i))
	  (sense (get-row-sense lp i)))
      (unless (row-is-removed row_i)
	(let ((has-li (aref (presolve-bounds-row-has-l bounds) i))
	      (has-ui (aref (presolve-bounds-row-has-u bounds) i))
	      (li (aref (presolve-bounds-row-l bounds) i))
	      (ui (aref (presolve-bounds-row-u bounds) i))
	      (bi (row-b row_i)))
	  (when (or has-li has-ui)
	    (dorow (lp i :visit-slack nil :j j :index index :col col_j)
	      (let ((a_ij (aref (column-values col_j) index))
		    (has-lj (column-has-l col_j))
		    (has-uj (column-has-u col_j))
		    (lj (column-l col_j))
		    (uj (column-u col_j)))
		(symbol-macrolet
		    ((imp-has-l (aref (presolve-bounds-col-has-l bounds) j))
		     (imp-has-u (aref (presolve-bounds-col-has-u bounds) j))
		     (imp-l (aref (presolve-bounds-col-l bounds) j))
		     (imp-u (aref (presolve-bounds-col-u bounds) j)))
		  (unless (= -1 sense)
		    (when (and has-li (< a_ij 0) has-uj)
		      (if imp-has-l
			  (maxf imp-l (+ uj (/ (- bi li) a_ij)))
			  (setf imp-l (+ uj (/ (- bi li) a_ij))
				imp-has-l t)))
		    (when (and has-li (> a_ij 0) has-lj)
		      (if imp-has-u
			  (maxf imp-u (+ lj (/ (- bi li) a_ij)))
			  (setf imp-u (+ lj (/ (- bi li) a_ij))
				imp-has-u t))))
		  (unless (= 1 sense)
		    (when (and has-ui (> a_ij 0) has-uj)
		      (if imp-has-l
			  (maxf imp-l (+ uj (/ (- bi ui) a_ij)))
			  (setf imp-l (+ uj (/ (- bi ui) a_ij))
				imp-has-l t)))
		    (when (and has-ui (< a_ij 0) has-lj)
		      (if imp-has-u
			  (maxf imp-u (+ lj (/ (- bi ui) a_ij)))
			  (setf imp-u (+ lj (/ (- bi ui) a_ij))
				imp-has-u t)))))))))))))
		      
  

;;;; Computes all bounds
(defun compute-bounds (lp bounds)
  (compute-row-bounds lp bounds)
  (compute-column-bounds lp bounds))
    


;;;; Cleanly removes column from LP
(defun remove-column (lp col)
  (setf (column-is-removed col) t)
  (if (< 0 (* (lp-obj-sense lp) (column-c col)))
      (if (column-has-u col)
	  (if (column-has-l col)
	      (setf (column-l col) (column-u col))
	      (setf (lp-is-unbounded lp) t))
	  (setf (lp-is-unbounded lp) t))
      (if (column-has-l col)
	  (if (column-has-u col)
	      (setf (column-u col) (column-l col))
	      (setf (lp-is-unbounded lp) t))
	  (setf (lp-is-unbounded lp) t))))
	


;;;;
(defun fix-column (col val)
  (setf (column-is-removed col) t
	(column-l col) val
	(column-u col) val
	(column-has-l col) t
	(column-has-u col) t))


(defun update-column-lower-bound (lp col val)
  (if (column-has-l col)
      (when (< (column-l col) val)
	(setf (column-l col) val))
      (setf (column-has-l col) t
	    (column-l col) val))
  (when (and (column-has-u col)
	     (< (column-u col) (column-l col)))
    (setf (lp-is-infeasible lp) t)))
    
  

(defun update-column-upper-bound (lp col val)
  (if (column-has-u col)
      (when (< (column-u col) val)
	(setf (column-u col) val))
      (setf (column-has-u col) t
	    (column-u col) val))
  (when (and (column-has-l col)
	     (< (column-u col) (column-l col)))
    (setf (lp-is-infeasible lp) t)))


(defun get-row-sense (lp i)
  (let* ((row_i (aref (lp-rows lp) i))
	 (slack-index (row-slack-index row_i)))
    (if (= -1 slack-index) 
	0
	(aref (column-values (aref (lp-columns lp) slack-index)) 0))))      
		  


;;;;
(defun singleton-row (lp i only-j)
  (let* ((row_i    (aref (lp-rows lp) i))
	 (col_j    (aref (lp-columns lp) only-j))
	 (index    (find-index (column-indices col_j) i))
	 (sense    (get-row-sense lp i))
	 (constant (row-b row_i)))
    (dorow (lp i :visit-only-fixed t :a a_ij :u v)
      (decf constant (* a_ij v)))
    (divf constant (aref (column-values col_j) index))
    (if (zerop sense)
	(fix-column col_j constant)
	(if (< 0 (* sense (aref (column-values col_j) index)))
	    (update-column-upper-bound lp col_j constant)
	    (update-column-lower-bound lp col_j constant)))))
    
		
	
		  
;;;;
(defun presolve-empty-and-singleton-rows (lp)
  (let ((m (length (lp-rows lp)))
	(is-modified nil))
    (dotimes (i m is-modified)
      (let ((row_i     (aref (lp-rows lp) i))
	    (only-j    -1)
	    (row-count 0))
	(unless (or (lp-is-infeasible lp)
		    (row-is-removed row_i))
	  (dorow (lp i :visit-slack nil :j j :index index :col col_j)
	    (setf only-j j)
	    (incf row-count)
	    (when (= 2 row-count)
	      (return)))
	  (when (< row-count 2)
	    (setf is-modified t
		  (row-is-removed row_i) t)
	    (when (= 1 row-count)
	      (singleton-row lp i only-j))))))))



;;;;
(defun presolve-empty-columns (lp)
  (let ((n (length (lp-columns lp)))
	(is-modified nil))
    (dotimes (j n is-modified)
      (let ((col_j  (aref (lp-columns lp) j))
	    (col-count 0))
	(unless (or (lp-is-infeasible lp)
		    (column-is-removed col_j))
	  (cond 
	    ((and (column-has-l col_j) (column-has-u col_j)
		  (= (column-l col_j) (column-u col_j)))
	     (setf (column-is-removed col_j) t))
	    (t
	     (docol (lp col_j)
	       (incf col-count)
	       (return))
	     (when (zerop col-count)
	       (when (or (and (not (column-has-u col_j))
			      (< 0 (* (column-c col_j) (lp-obj-sense lp))))
			 (and (not (column-has-l col_j))
			      (> 0 (* (column-c col_j) (lp-obj-sense lp)))))
		 (setf (lp-is-unbounded lp) t
		       (lp-is-infeasible lp) t))
	       (remove-column lp col_j)
	       (setf is-modified t)))))))))
    



;;;;
(defun presolve-redundant-and-forcing-rows (lp bounds)
  (let ((is-modified nil))
    (dotimes (i (presolve-bounds-m bounds) is-modified)
      (let ((row_i (aref (lp-rows lp) i))
	    (sense (get-row-sense lp i))
	    (has-li (aref (presolve-bounds-row-has-l bounds) i))
	    (has-ui (aref (presolve-bounds-row-has-u bounds) i))
	    (li (aref (presolve-bounds-row-l bounds) i))
	    (ui (aref (presolve-bounds-row-u bounds) i)))
	(unless (or (lp-is-infeasible lp)
		    (row-is-removed row_i))
	  (let ((bi (row-b row_i)))
	    (unless (= -1 sense)
	      (when has-li
		(when (< bi li)
		  (setf (lp-is-infeasible lp) t
			is-modified t))
		(when (= bi li)
		  (dorow (lp i :index index :col col_j)
		    (if (< 0 (aref (column-values col_j) index))
			(setf (column-u col_j) (column-l col_j)
			      (column-has-u col_j) t)
			(setf (column-l col_j) (column-u col_j)
			      (column-has-l col_j) t))
		    (setf (column-is-removed col_j) t))
		  (setf (row-is-removed row_i) t
			is-modified t)))
	      (when (and has-ui (<= ui bi))
		(setf (row-is-removed row_i) t
		      is-modified t)))
	    (unless (= 1 sense)
	      (when has-ui
		(when (< ui bi)
		  (setf (lp-is-infeasible lp) t
			is-modified t))
		(when (= bi ui)
		  (dorow (lp i :index index :col col_j)
		    (if (> 0 (aref (column-values col_j) index))
			(setf (column-u col_j) (column-l col_j)
			      (column-has-u col_j) t)
			(setf (column-l col_j) (column-u col_j)
			      (column-has-l col_j) t))
		    (setf (column-is-removed col_j) t))
		  (setf (row-is-removed row_i) t
			is-modified t)))
	      (when (and has-li (<= bi li))
		(setf (row-is-removed row_i) t
		      is-modified t)))))))))
  


;;;;
(defun presolve-tighten-column-bounds (lp bounds)
  (let ((is-modified nil))
    (dotimes (j (length (lp-columns lp)) is-modified)
      (let ((col_j (aref (lp-columns lp) j))
	    (imp-has-lj (aref (presolve-bounds-col-has-l bounds) j))
	    (imp-has-uj (aref (presolve-bounds-col-has-u bounds) j))
	    (imp-lj (aref (presolve-bounds-col-l bounds) j))
	    (imp-uj (aref (presolve-bounds-col-u bounds) j)))
	(unless (or (lp-is-infeasible lp) 
		    (column-is-removed col_j))
	  (when imp-has-lj
	    (if (column-has-l col_j)
		(when (< (column-l col_j) imp-lj)
		  (setf (column-l col_j) imp-lj
			is-modified t))
		(setf (column-l col_j) imp-lj
		      (column-has-l col_j) t
		      is-modified t)))
	  (when imp-has-uj
	    (if (column-has-u col_j)
		(when (< imp-uj (column-u col_j))
		  (setf (column-u col_j) imp-uj
			is-modified t))
		(setf (column-u col_j) imp-uj
		      (column-has-u col_j) t
		      is-modified t)))
	  (when (and (column-has-l col_j)
		     (column-has-u col_j))
	    (when (= (column-l col_j) (column-u col_j))
	      (setf (column-is-removed col_j) t))
	    (when (> (column-l col_j) (column-u col_j))
	      (setf (lp-is-infeasible lp) t
		    is-modified t))))))))
	
	    


;;;;
(defun presolve-implied-free-columns (lp bounds)
  (let ((is-modified nil))
    (dotimes (j (length (lp-columns lp)) is-modified)
      (let ((col_j (aref (lp-columns lp) j))
	    (imp-has-lj (aref (presolve-bounds-col-has-l bounds) j))
	    (imp-has-uj (aref (presolve-bounds-col-has-u bounds) j))
	    (imp-lj (aref (presolve-bounds-col-l bounds) j))
	    (imp-uj (aref (presolve-bounds-col-u bounds) j)))
	(unless (or (lp-is-infeasible lp) 
		    (column-is-removed col_j))
	  (when (or (column-has-l col_j) (column-has-u col_j))
	    (unless (and (not imp-has-lj) (column-has-l col_j))
	      (unless (and imp-has-lj (column-has-l col_j)
			   (< imp-lj (column-l col_j)))
		(unless (and (not imp-has-uj) (column-has-u col_j))
		  (unless (and imp-has-uj (column-has-u col_j)
			       (> imp-uj (column-u col_j)))
		    ;; TODO: complete this
		    ))))))))))



;;;;
(defun presolve-pass (lp)
  (let* ((is-modified nil))
    ;; check for infeasability
    (when (lp-is-infeasible lp)
      (return-from presolve-pass nil))
    ;; check for empty and singleton rows
    (orf is-modified
	 (presolve-empty-and-singleton-rows lp))
    (when (lp-is-infeasible lp)
      (return-from presolve-pass nil))
    ;; check for empty columns
    (orf is-modified
	 (presolve-empty-columns lp))
    (when (lp-is-infeasible lp)
      (return-from presolve-pass nil))
    ;; compute bounds for each row and each column
    (when *presolve-with-bounds* 
      (let ((bounds (init-bounds (length (lp-rows lp)) (length (lp-columns lp)))))
	(compute-bounds lp bounds)
	;; check for redundant and forcing rows
	(orf is-modified
	     (presolve-redundant-and-forcing-rows lp bounds))
	(when (lp-is-infeasible lp)
	  (return-from presolve-pass nil))
	;; deal with implied free columns
	(orf is-modified
	     (presolve-implied-free-columns lp bounds))
	(when (lp-is-infeasible lp)
	  (return-from presolve-pass nil))
	;; tighten column bounds
	(orf is-modified
	     (presolve-tighten-column-bounds lp bounds))
	(when (lp-is-infeasible lp)
	  (return-from presolve-pass nil))))
    ;; return true if LP was modified
    is-modified))

	 

;;;; Rescales LP between -1 and 1
(defun rescale (lp)
  (let ((n (length (lp-columns lp)))
	(m (length (lp-rows lp))))
    ;; rescale rows
    (dotimes (i m)
      (let ((row_i (aref (lp-rows lp) i)))
	(unless (row-is-removed row_i)
	  (let ((scal  0))
	    (dorow (lp i :index index :col col_j)
	      (absmaxf scal (aref (column-values col_j) index)))
	    (divf (row-b row_i) scal)
	    (divf (row-squash-coef row_i) scal)))))
    (dotimes (j n)
      (let ((col_j (aref (lp-columns lp) j)))
	(dotimes (index (length (column-indices col_j)))
	  (let* ((i (aref (column-indices col_j) index))
		 (row_i (aref (lp-rows lp) i)))
	    (mulf (aref (column-values col_j) index) (abs (row-squash-coef row_i)))))))
    ;; rescale columns
    (dotimes (j n)
      (let ((col_j (aref (lp-columns lp) j)))
	(unless (column-is-removed col_j)
	  (let ((scal  0))
	    (docol (lp col_j :index index)
	      (absmaxf scal (aref (column-values col_j) index)))
	    (divf (column-c col_j) scal)
	    (dotimes (index (length (column-indices col_j)))
	      (divf (aref (column-values col_j) index) scal))
	    (when (column-has-l col_j)
	      (divf (column-l col_j) scal))
	    (when (column-has-u col_j)
	      (divf (column-u col_j) scal))))))))
	


;;;;
(defun print-preprocess-stats (lp)
  (let ((m  (length (lp-rows lp)))
	(n  (length (lp-columns lp)))
	(mr 0)
	(nr 0))
    (dotimes (i m)
      (when (row-is-removed (aref (lp-rows lp) i))
	(incf mr)))
    (dotimes (j n)
      (when (column-is-removed (aref (lp-columns lp) j))
	(incf nr)))
    (format t "Preprocessing over.~%~A out of ~A rows removed.~%~A out of ~A columns removed.~%~%" mr m nr n)))



;;;;
(defun preprocess (lp)
  (loop
     (when (lp-is-infeasible lp)
       (return))
     (unless (presolve-pass lp)
       (when (lp-is-infeasible lp)
	 (return))
       (rescale lp)
       (return)))
  (if (lp-is-infeasible lp)
      (if (lp-is-unbounded lp)
	  (format t "Preprocessing over.~%LP unbounded.~%~%")
	  (format t "Preprocessing over.~%LP infeasible.~%~%"))
      (print-preprocess-stats lp))
  (not (lp-is-infeasible lp)))
     
	 
	 





