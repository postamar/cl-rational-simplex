;;;;; Preassigned pivot procedure



(defun pip-sort-ref-buffer (ssqm)
  (let ((ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(row-count (sparse-square-matrix-row-count ssqm)))
    (setf (sparse-square-matrix-ref-buffer ssqm)
	  (sort ref-buffer #'(lambda (r1 r2)
			       (< (aref row-count r1) (aref row-count r2)))))))



;;;; Initializes pip data structures and isolates pivotable columns
(defun pip-init (ssqm)
  (let ((col-count (sparse-square-matrix-col-count ssqm))
	(row-count (sparse-square-matrix-row-count ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(col-avail (sparse-square-matrix-col-avail ssqm))
	(ref-perm-col (sparse-square-matrix-ref-perm-col ssqm))
	(ref-perm-row (sparse-square-matrix-ref-perm-row ssqm))
	(ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(col-row-refs (sparse-square-matrix-col-row-ref ssqm))
	(row-col-refs (sparse-square-matrix-row-col-ref ssqm))
	(m (sparse-square-matrix-size ssqm))
	(mu (- (sparse-square-matrix-size ssqm) 1)))
    ;; initialize data structures
    (dotimes (k m)
      (setf (aref row-avail k) t
	    (aref col-avail k) t
	    (aref ref-buffer k) k
	    (aref row-count k) (length (aref row-col-refs k))
	    (aref col-count k) (length (aref col-row-refs k)))) 
    ;; scan column counts
    (loop 
       (let ((col-index (dotimes (k m -1)
			  (when (and (aref col-avail k)
				     (= 1 (aref col-count k)))
			    (return k)))))
	 (when (= -1 col-index)
	   (return mu))
	 (switch-columns-in-matrix ssqm mu (aref ref-perm-col col-index))
	 (setf (aref col-count col-index) most-positive-fixnum
	       (aref col-avail col-index) nil)
	 (let* ((col-row-ref (aref col-row-refs col-index))
		(row-index (dotimes (k (length col-row-ref) -1)
			     (when (aref row-avail (aref col-row-ref k))
			       (return (aref col-row-ref k)))))
		(row-col-ref (aref row-col-refs row-index)))
	   (switch-rows-in-matrix ssqm mu (aref ref-perm-row row-index))
	   (setf (aref row-count row-index) most-positive-fixnum
		 (aref row-avail row-index) nil)
	   (dotimes (k (length row-col-ref))
	     (let ((col-ref (aref row-col-ref k)))
	       (when (aref col-avail col-ref)
		 (when (zerop (decf (aref col-count col-ref)))
		   (setf (sparse-square-matrix-is-singular ssqm) 'redundant-column
			 (sparse-square-matrix-singular-ref ssqm) col-ref)
		   (return -1)))))))
       (decf mu)
       (when (<= mu 0)
	 (return mu)))))



;;;;
(defun pip-pivot-row (ssqm nu pivot-row-ref pivot-col-ref )
  (let ((row-count (sparse-square-matrix-row-count ssqm))
	(col-count (sparse-square-matrix-col-count ssqm))
	(col-avail (sparse-square-matrix-col-avail ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(ref-perm-col (sparse-square-matrix-ref-perm-col ssqm))
	(ref-perm-row (sparse-square-matrix-ref-perm-row ssqm))
	(col-row-ref (aref (sparse-square-matrix-col-row-ref ssqm) pivot-col-ref)))
    (switch-rows-in-matrix ssqm nu (aref ref-perm-row pivot-row-ref))
    (setf (aref row-count pivot-row-ref) most-positive-fixnum
	  (aref row-avail pivot-row-ref) nil)
    (switch-columns-in-matrix ssqm nu (aref ref-perm-col pivot-col-ref))
    (setf (aref col-count pivot-col-ref) most-positive-fixnum
	  (aref col-avail pivot-col-ref) nil)
    (dotimes (k (length col-row-ref) (+ 1 nu))
      (let ((row-ref (aref col-row-ref k)))
	(when (aref row-avail row-ref)
	  (when (zerop (decf (aref row-count row-ref)))
	    (setf (sparse-square-matrix-is-singular ssqm) 'redundant-row
		  (sparse-square-matrix-singular-ref ssqm) row-ref)
	    (return -1)))))))
    


;;;;
(defun pip-scan-row-counts (ssqm mu nu l)
  (let ((row-count (sparse-square-matrix-row-count ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(col-avail (sparse-square-matrix-col-avail ssqm))
	(row-col-refs (sparse-square-matrix-row-col-ref ssqm))
	(m (sparse-square-matrix-size ssqm)))
    (loop
       (let* ((row-ref (dotimes (k m -1)
			 (when (and (aref row-avail k) 
				    (= 1 (aref row-count k)))
			   (return k)))))
	 (when (= -1 row-ref)
	   (return nu))
	 (let* ((row-col-ref (aref row-col-refs row-ref))
		(col-ref (dotimes (k (length row-col-ref) -1)
			   (when (aref col-avail (aref row-col-ref k))
			     (return (aref row-col-ref k))))))
	   (when (= -1 (pip-pivot-row ssqm nu row-ref col-ref))
	     (return -1))))
       (incf nu)
       (when (< mu nu)
	 (return nu))
       (unless (zerop l)
	 (return nu)))))



;;;; tally function, t_k(n)
(defun pip-tally (ssqm k col-ref)
  (let ((row-count (sparse-square-matrix-row-count ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(col-row-ref (aref (sparse-square-matrix-col-row-ref ssqm) col-ref))
	(val 0))
    (dotimes (index (length col-row-ref) val)
      (let ((row-ref (aref col-row-ref index)))
	(when (and (aref row-avail row-ref)
		   (<= (aref row-count row-ref) k))
	  (incf val))))))



;;;;
(defun pip-find-spike (ssqm l)
  (let ((ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(col-avail (sparse-square-matrix-col-avail ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(flag-buffer (sparse-square-matrix-flag-buffer ssqm))
	(row-count (sparse-square-matrix-row-count ssqm))
	(col-count (sparse-square-matrix-col-count ssqm))
	(spikes (sparse-square-matrix-spikes ssqm))
	(m (sparse-square-matrix-size ssqm))
	(spike-col-ref -1)
	(k most-positive-fixnum))
    (setf k (aref row-count (aref ref-buffer 0)))
    (dotimes (col-ref m)
      (setf (aref flag-buffer col-ref) (aref col-avail col-ref)))
    ;; steps 6 and 7
    (loop
       (let ((max-tally-val 0)
	     (max-tally-max-col-count 0)
	     (max-tally-col -1))
	 (dotimes (col-ref m)
	   (when (aref flag-buffer col-ref)
	     (let ((tally (pip-tally ssqm k col-ref)))
	       (when (or (< max-tally-val tally)
			 (and (= tally max-tally-val)
			      (< max-tally-max-col-count (aref col-count col-ref))))
		 (setf max-tally-val tally
		       max-tally-max-col-count (aref col-count col-ref)
		       max-tally-col col-ref)))))
	 (when (< 1 max-tally-val)
	   (setf spike-col-ref max-tally-col)
	   (return))
	 (assert (= 1 max-tally-val))
	 (dotimes (col-ref m)
	   (when (aref col-avail col-ref)
	     (setf (aref flag-buffer col-ref) (= 1 (pip-tally ssqm k col-ref)))))
	 (dotimes (ref m)
	   (when (< k (aref row-count (aref ref-buffer ref)))
	     (setf k (aref row-count (aref ref-buffer ref)))
	     (return)))))
    ;; step 8
    (assert (/= -1 spike-col-ref))
    (setf (aref spikes l) spike-col-ref
	  (aref col-avail spike-col-ref) nil
	  (aref col-count spike-col-ref) most-positive-fixnum)
    (let ((col-row-ref (aref (sparse-square-matrix-col-row-ref ssqm) spike-col-ref)))
      (dotimes (ref (length col-row-ref) (progn (pip-sort-ref-buffer ssqm) t))
	(let ((row-ref (aref col-row-ref ref)))
	  (when (aref row-avail row-ref)
	    (when (zerop (decf (aref row-count row-ref)))
	      (setf (sparse-square-matrix-is-singular ssqm) 'redundant-row
		    (sparse-square-matrix-singular-ref ssqm) row-ref)
	      (return nil))))))))
	    
		   
		
;;;; Pivots the spike
;;;; returns new nu
(defun pip-spike-pivot (ssqm mu nu l)
  (let ((ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(col-avail (sparse-square-matrix-col-avail ssqm))
	(row-avail (sparse-square-matrix-row-avail ssqm))
	(flag-buffer (sparse-square-matrix-flag-buffer ssqm))
	(row-count (sparse-square-matrix-row-count ssqm))
	(col-count (sparse-square-matrix-col-count ssqm))
	(ref-perm-col (sparse-square-matrix-ref-perm-col ssqm))
	(ref-perm-row (sparse-square-matrix-ref-perm-row ssqm))
	(col-row-refs (sparse-square-matrix-col-row-ref ssqm))
	(spikes (sparse-square-matrix-spikes ssqm))
	(m (sparse-square-matrix-size ssqm))
	(pivot-col-ref -1)
	(k 1))
    (dotimes (col-ref m)
      (setf (aref flag-buffer col-ref) (aref col-avail col-ref)))
    ;; steps 10 and 11
    (loop
       (let ((max-tally-val 0)
	     (max-tally-max-col-count 0)
	     (max-tally-col -1))
	 (dotimes (col-ref m)
	   (when (aref flag-buffer col-ref)
	     (let ((tally (pip-tally ssqm k col-ref)))
	       (when (or (< max-tally-val tally)
			 (and (= tally max-tally-val)
			      (< max-tally-max-col-count (aref col-count col-ref))))
		 (setf max-tally-val tally
		       max-tally-max-col-count (aref col-count col-ref)
		       max-tally-col col-ref)))))
	 (when (< 1 max-tally-val)
	   (setf pivot-col-ref max-tally-col)
	   (return))
	 (assert (= 1 max-tally-val))
	 (dotimes (col-ref m)
	   (when (aref col-avail col-ref)
	     (setf (aref flag-buffer col-ref) (= 1 (pip-tally ssqm k col-ref)))))
	 (dotimes (ref m)
	   (when (< k (aref row-count (aref ref-buffer ref)))
	     (setf k (aref row-count (aref ref-buffer ref)))
	     (return)))))
    (let ((q (pip-tally ssqm 1 pivot-col-ref)))
      ;; steps 12 and 13
      (let ((pivot-row-ref -1)
	    (col-row-ref (aref col-row-refs pivot-col-ref)))
	(dotimes (ref (length col-row-ref))
	  (let ((row-ref (aref col-row-ref ref)))
	    (when (and (aref row-avail row-ref)
		       (= 1 (aref row-count row-ref)))
	      (setf pivot-row-ref row-ref))))
	(switch-rows-in-matrix ssqm nu (aref ref-perm-row pivot-row-ref))
	(setf (aref row-count pivot-row-ref) most-positive-fixnum
	      (aref row-avail pivot-row-ref) nil)
	(switch-columns-in-matrix ssqm nu (aref ref-perm-col pivot-col-ref))
	(setf (aref col-count pivot-col-ref) most-positive-fixnum
	      (aref col-avail pivot-col-ref) nil)
	(dotimes (ref (length col-row-ref))
	  (let ((row-ref (aref col-row-ref ref)))
	    (when (aref row-avail row-ref)
	      (decf (aref row-count row-ref)))))
	(pip-sort-ref-buffer ssqm)
	(loop 
	   (incf nu)
	   (when (< mu nu)
	     (pip-sort-ref-buffer ssqm)
	     (return (values nu l)))
	   (decf q)
	   (when (zerop q)
	     (pip-sort-ref-buffer ssqm)
	     (return (values nu l)))
	   (decf l)
	   (switch-columns-in-matrix ssqm nu (aref ref-perm-col (aref spikes l)))
	   (let ((spike-row-ref -1))
	     (dotimes (ref m)
	       (let ((row-ref (aref ref-buffer ref)))
		 (when (zerop (aref row-count row-ref))
		   (setf spike-row-ref row-ref)
		   (return))))
	     (assert (/= -1 spike-row-ref))
	     (switch-rows-in-matrix ssqm nu (aref ref-perm-row spike-row-ref))))))))
	     

	       
		      
;;;;
(defun preassigned-pivot (ssqm)
  ;; step 1
  (let ((row-count (sparse-square-matrix-row-count ssqm))
	(ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(m (sparse-square-matrix-size ssqm))
	(mu (pip-init ssqm))
	(nu 0)
	(l 0))
    (pip-sort-ref-buffer ssqm)
    (cond ((= -1 mu)) ;; singularity, do nothing
	  ((zerop mu) ;; find the first spike
	   (pip-find-spike ssqm l)
	   (incf l))
	  (t ;; scan row counts for the first time
	   (let ((new-nu (pip-scan-row-counts ssqm mu nu l)))
	     (unless (= nu new-nu)
	       (pip-sort-ref-buffer ssqm))
	     (setf nu new-nu))))
    ;; deal with bumps
    (unless (or (sparse-square-matrix-is-singular ssqm) 
		(< mu nu))
      (loop
;	 (print-sparse-square-matrix ssqm)
	 ;; scan row counts for path-decisions
	 (cond ((< 1 (aref row-count (aref ref-buffer 0)))
		;; another spike
;		(format t "finding spike~%")
		(unless (pip-find-spike ssqm l)
		  (return))
;		(format t "found column ~A~%" (aref (sparse-square-matrix-ref-perm-col ssqm) (aref (sparse-square-matrix-spikes ssqm) l)))
		(incf l))
	       ((or (zerop l)
		    (and (= 1 (aref row-count (aref ref-buffer 0)))
			 (< 1 (aref row-count (aref ref-buffer 1)) (+ m 1))))
		;; scan row count
;		(format t "scannning row count~%")
		(let ((new-nu (pip-scan-row-counts ssqm mu nu l)))
		  (unless (= nu new-nu)
		    (pip-sort-ref-buffer ssqm))
		  (setf nu new-nu))
		(when (< mu nu)
		  (return)))
	       ((= 1 (aref row-count (aref ref-buffer 0)))
		;; spike pivot
;		(format t "spike pivot~%")
		(multiple-value-bind (new-nu new-l) (pip-spike-pivot ssqm mu nu l)
		  (setf nu new-nu l new-l))
		(when (< mu nu)
		  (return)))
	       (t
		(error "error in preassigned pivot procedure~%")))))
    ;; return t when done with pivoting
    ;; return nil if singular
 ;   (print-sparse-square-matrix ssqm)
    (< mu nu)))


   

		 

