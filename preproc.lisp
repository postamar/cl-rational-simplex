;;;;; LP preprocessing



;;;;
(defun fix-column (lp col val)
  (update-column-lower-bound lp col val)
  (update-column-upper-bound lp col val)
  (lp-remove-column lp (column-ref col)))


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


(defun get-row-sense (lp row-ref)
  (let* ((row (aref (lp-rows lp) row-ref))
	 (slack-col-ref (row-slack-col-ref row)))
    (if (= -1 slack-col-ref) 
	0
	(aref (column-values (aref (lp-columns lp) slack-col-ref)) 0))))
		  


;;;;
(defun presolve-singleton-row-fix (lp row-ref only-col-ref)
  (let ((col (aref (lp-columns lp) only-col-ref))
	(row (aref (lp-rows lp) row-ref)))
    (let ((b (row-b row))
	  (a (rational-in-column col (aref (row-col-indices row) 0))))
      (fix-column lp col (/ b a)))))
	
  
(defun presolve-singleton-row-bound (lp row-ref only-col-ref)
  (let* ((col (aref (lp-columns lp) only-col-ref))
	 (row (aref (lp-rows lp) row-ref))
	 (slack-col (aref (lp-columns lp) (row-slack-col-ref row)))
	 (slack-col-coef (aref (column-values slack-col) 0)))
    (let ((b (row-b row))
	  (a (rational-in-column col (aref (row-col-indices row) 0))))
      (if (= -1 slack-col-coef)
	  (update-column-upper-bound lp col (/ b a))
	  (update-column-lower-bound lp col (/ b a))))))


		  
;;;;
(defun presolve-empty-and-singleton-rows (lp)
  (let* ((active-row-refs (lp-active-row-refs lp))
	 (m (length active-row-refs))
	 (is-modified nil))
    (dotimes (i m is-modified)
      (when (lp-is-infeasible lp)
	(return))
      (let* ((row-ref       (aref active-row-refs i))
	     (row           (aref (lp-rows lp) row-ref))
	     (slack-col-ref (row-slack-col-ref row))
	     (row-col-refs  (row-col-refs row))
	     (row-count     (length row-col-refs )))
	(cond ((or (zerop row-count)
		   (and (= 1 row-count)
			(/= -1 slack-col-ref)))
	       ;; empty column
	       (setf is-modified t)
	       (lp-remove-row lp row-ref))
	      ((= 1 row-count)
	       ;; fix column to rhs
	       (presolve-singleton-row-fix lp row-ref (aref row-col-refs 0))
	       (setf is-modified t)
	       (lp-remove-row lp row-ref))
	      ((and (= 2 row-count)
		    (/= -1 slack-col-ref))
	       ;; fix column bound to rhs
	       (if (= slack-col-ref (aref row-col-refs 0))
		   (presolve-singleton-row-bound lp row-ref (aref row-col-refs 1))
		   (presolve-singleton-row-bound lp row-ref (aref row-col-refs 0)))
	       (setf is-modified t)
	       (lp-remove-row lp row-ref)))))))
	      
	       

;;;;
(defun presolve-empty-columns (lp)
  (let* ((active-col-refs (lp-active-col-refs lp))
	 (n (length active-col-refs))
	 (is-modified nil))
    (dotimes (j n is-modified)
      (let* ((col-ref (aref active-col-refs j))
	     (col (aref (lp-columns lp) col-ref))
	     (col-count 0))
	(when (lp-is-infeasible lp)
	  (return))
	(cond 
	  ((and (column-has-l col) (column-has-u col)
		(= (column-l col) (column-u col)))
	   (lp-remove-column lp col-ref))
	  (t
	   (docol (lp col)
	     (incf col-count)
	     (return))
	   (when (zerop col-count)
	     (when (or (and (not (column-has-u col))
			    (< 0 (* (column-c col) (lp-obj-sense lp))))
		       (and (not (column-has-l col))
			    (> 0 (* (column-c col) (lp-obj-sense lp)))))
	       (setf (lp-is-unbounded lp) t
		     (lp-is-infeasible lp) t))
	     (lp-remove-column lp col-ref)
	     (setf is-modified t))))))))
    


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
    ;; return true if LP was modified
    is-modified))

	 

;;;;
(defun print-preprocess-stats (lp)
  (let ((m  (length (lp-rows lp)))
	(n  (length (lp-columns lp)))
	(mr 0)
	(nr 0))
    (dotimes (i m)
      (unless (row-is-active (aref (lp-rows lp) i))
	(incf mr)))
    (dotimes (j n)
      (unless (column-is-active (aref (lp-columns lp) j))
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
       (return)))
  (if (lp-is-infeasible lp)
      (if (lp-is-unbounded lp)
	  (format t "Preprocessing over.~%LP unbounded.~%~%")
	  (format t "Preprocessing over.~%LP infeasible.~%~%"))
      (print-preprocess-stats lp))
  (not (lp-is-infeasible lp)))
     
	 
	 





