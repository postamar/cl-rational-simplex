;;;;; The Preassigned Pivot Procedure


;;;; Auxilliary functions

(defun pip-sort-refs (b)
  (let ((refs (basis-refs b)))
    (setf (basis-refs b)
	  (sort refs
		#'(lambda (r1 r2)
		    (or (<= (basis-size b) r2)
			(and (< r1 (basis-size b))
			     (< (aref (basis-row-counts b) r1) 
				(aref (basis-row-counts b) r2)))))))))


;;;;; Initialization

(defun pip-init (b)
  (dotimes (k (basis-size b))
    (when (zerop (aref (basis-row-counts b) k))
      (setf (basis-is-singular b)  'redundant-row
	    (basis-singular-ref b) k)
      (return))
    (when (zerop (aref (basis-col-counts b) k))
      (setf (basis-is-singular b)  'redundant-column
	    (basis-singular-ref b) k)
      (return)))
  ;; scan column counts
  (loop 
     (when (basis-is-singular b)
       (return))
     (when (zerop (decf (basis-pip-last b)))
       (return t))
     (let ((j -1))
       ;; find column with count of 1
       (dotimes (k (basis-size b))
	 (when (and (aref (basis-col-avails b) k)
		    (= 1 (aref (basis-col-counts b) k)))
	   (setf j k)
	   (return)))
       (when (= -1 j)
	 (return t))
       ;; permutate column
       (basis-permutate-columns b j (aref (basis-pj->j b) (basis-pip-last b)))
       (setf (aref (basis-col-counts b) j) most-positive-fixnum
	     (aref (basis-col-avails b) j) nil)
       ;; find pivot row
       (let ((i -1)
	     (col-is (lu-eta-matrix-is (aref (basis-l-columns b) j))))
	 (dotimes (k (length col-is))
	   (when (aref (basis-row-avails b) (aref col-is k))
	     (setf i (aref col-is k))
	     (return)))
	 ;; permutate row
	 (basis-permutate-rows b i (aref (basis-pi->i b) (basis-pip-last b)))
	 (setf (aref (basis-row-counts b) i) most-positive-fixnum
	       (aref (basis-row-avails b) i) nil)
	 ;; decrease column counts
	 (let ((row-js (aref (basis-row-js b) i)))
	   (dotimes (k (length row-js))
	     (let ((j (aref row-js k)))
	       (when (aref (basis-col-avails b) j)
		 (when (zerop (decf (aref (basis-col-counts b) j)))
		   (setf (basis-is-singular b)  'redundant-column
			 (basis-singular-ref b) j)
		   (return))))))))))
     


;;;;
(defun pip-scan-row-counts (b)
  (loop 
     (when (basis-is-singular b)
       (return))
     (when (< (basis-pip-last b) (basis-pip-first b))
       (return t))
     (let ((i -1))
       ;; find row with count of 1
       (dotimes (k (basis-size b))
	 (when (and (aref (basis-row-avails b) k)
		    (= 1 (aref (basis-row-counts b) k)))
	   (setf i k)
	   (return)))
       (when (= -1 i)
	 (return t))
       ;; permutate row
       (basis-permutate-rows b i (aref (basis-pi->i b) (basis-pip-first b)))
       (setf (aref (basis-row-counts b) i) most-positive-fixnum
	     (aref (basis-row-avails b) i) nil)
       ;; find pivot column
       (let ((j -1)
	     (row-js (aref (basis-row-js b) i)))
	 (dotimes (k (length row-js))
	   (when (aref (basis-col-avails b) (aref row-js k))
	     (setf j (aref row-js k))
	     (return)))
	 ;; permutate column 
	 (basis-permutate-columns b j (aref (basis-pj->j b) (basis-pip-first b)))
	 (setf (aref (basis-col-counts b) j) most-positive-fixnum
	       (aref (basis-col-avails b) j) nil)
	 ;; decrease row counts
	 (let ((col-is (lu-eta-matrix-is (aref (basis-l-columns b) j))))
	   (dotimes (k (length col-is))
	     (let ((i (aref col-is k)))
	       (when (aref (basis-row-avails b) i)
		 (when (zerop (decf (aref (basis-row-counts b) i)))
		   (setf (basis-is-singular b)  'redundant-row
			 (basis-singular-ref b) i)
		   (return))))))))
     (incf (basis-pip-first b))
     (unless (zerop (length (basis-pip-spikes b)))
       (return t))))



;;;; Tally function: t_k(j)
(defun pip-tally (b k j)
  (let ((row-is (lu-eta-matrix-is (aref (basis-l-columns b) j)))
	(val    0))
    (dotimes (index (length row-is) val)
      (let ((i (aref row-is index)))
	(when (and (aref (basis-row-avails b) i)
		   (<= (aref (basis-row-counts b) i) k))
	  (incf val))))))



(defun pip-find-spike-j (b k)
  (dotimes (j (basis-size b))
    (setf (aref (basis-flags b) j) (aref (basis-col-avails b) j)))
  (loop
     (let ((max-tally-val 0)
	   (max-col-count 0)
	   (max-tally-j  -1))
       (dotimes (j (basis-size b))
	 (when (aref (basis-flags b) j)
	   (let ((tally     (pip-tally b k j))
		 (col-count (aref (basis-col-counts b) j)))
	     (when (or (< max-tally-val tally)
		       (and (= tally max-tally-val)
			    (< max-col-count col-count)))
	       (setf max-tally-val tally
		     max-col-count col-count
		     max-tally-j   j)))))
       (when (< 1 max-tally-val)
	 (return max-tally-j))
       (assert (= 1 max-tally-val))
       (dotimes (j (basis-size b))
	 (when (aref (basis-col-avails b) j)
	   (setf (aref (basis-flags b) j) (= 1 (pip-tally b k j)))))
       (dotimes (index (basis-size b))
	 (let ((row-count (aref (basis-row-counts b) (aref (basis-refs b) index))))
	   (when (< k row-count)
	     (setf k row-count)
	     (return)))))))



;;;; Selects a spike
(defun pip-find-spike (b)
  (let* ((k (aref (basis-row-counts b) (aref (basis-refs b) 0)))
	 (spike-j (pip-find-spike-j b k)))
    ;; step 8
    (assert (/= -1 spike-j))
    (setf (aref (basis-spikes b) spike-j) (basis-pip-first b))
    (vector-push-extend spike-j (basis-pip-spikes b))
    (setf (aref (basis-col-avails b) spike-j) nil
	  (aref (basis-col-counts b) spike-j) most-positive-fixnum)
    (let ((col-is (lu-eta-matrix-is (aref (basis-l-columns b) spike-j))))
      (dotimes (index (length col-is) (progn (pip-sort-refs b) t))
	(let ((i (aref col-is index)))
	  (when (aref (basis-row-avails b) i)
	    (when (zerop (decf (aref (basis-row-counts b) i)))
	      (setf (basis-is-singular b) 'redundant-row
		    (basis-singular-ref b) i)
	      (return nil))))))))



;;;; Pivots the spike
(defun pip-spike-pivot (b)
  (dotimes (j (basis-size b))
    (setf (aref (basis-flags b) j) (aref (basis-col-avails b) j)))
  (let* ((pivot-j (pip-find-spike-j b 1))
	 (q       (pip-tally b 1 pivot-j))
	 (pivot-i -1)
	 (col-is  (lu-eta-matrix-is (aref (basis-l-columns b) pivot-j))))
    ;; steps 12 and 13
    (dotimes (index (length col-is))
      (let ((i (aref col-is index)))
	(when (and (aref (basis-row-avails b) i)
		   (= 1 (aref (basis-row-counts b) i)))
	  (setf pivot-i i))))
    (assert (/= -1 pivot-i))
    (basis-permutate-rows b pivot-i (aref (basis-pi->i b) (basis-pip-first b)))
    (setf (aref (basis-row-counts b) pivot-i) most-positive-fixnum
	  (aref (basis-row-avails b) pivot-i) nil)
    (basis-permutate-columns b pivot-j (aref (basis-pj->j b) (basis-pip-first b)))
    (setf (aref (basis-col-counts b) pivot-j) most-positive-fixnum
	  (aref (basis-col-avails b) pivot-j) nil)
    (dotimes (index (length col-is))
      (let ((i (aref col-is index)))
	(when (aref (basis-row-avails b) i)
	  (decf (aref (basis-row-counts b) i)))))
    (pip-sort-refs b)
    (loop
       (incf (basis-pip-first b))
       (when (< (basis-pip-last b) (basis-pip-first b))
	 (pip-sort-refs b)
	 (return))
       (decf q)
       (when (zerop q)
	 (pip-sort-refs b)
	 (return))
       (let ((spike-j (vector-pop (basis-pip-spikes b)))
	     (spike-i -1))
	 (basis-permutate-columns b spike-j (aref (basis-pj->j b) (basis-pip-first b)))
	 (dotimes (index (basis-size b))
	   (let ((i (aref (basis-refs b) index)))
	     (when (zerop (aref (basis-row-counts b) i))
	       (setf spike-i i)
	       (return))))
	 (assert (/= -1 spike-i))
	 (setf (aref (basis-row-counts b) spike-i) most-positive-fixnum
	       (aref (basis-row-avails b) spike-i) nil)
	 (basis-permutate-rows b spike-i (aref (basis-pi->i b) (basis-pip-first b)))))))
	       


;;;;
(defun preassigned-pivot (b)
  ;; start by scanning row and column counts for the first time
  (when (pip-init b)
    (if (zerop (basis-pip-last b))
	(pip-find-spike b)
	(when (pip-scan-row-counts b)
	  (pip-sort-refs b)))
    ;; deal with bumps or terminate
    (unless (or (basis-is-singular b)
		(< (basis-pip-last b) (basis-pip-first b)))
      (loop
	 ;; scan row counts for path-decisions
	 (let ((smallest-row-count (aref (basis-row-counts b) 
					 (aref (basis-refs b) 0)))
	       (next-smallest-row-count (aref (basis-row-counts b)
					      (aref (basis-refs b) 1))))
	   (cond ((< 1 smallest-row-count)
		  ;; another spike
		  (unless (pip-find-spike b)
		    (return)))
		 ((or (zerop (length (basis-pip-spikes b)))
		      (and (= 1 smallest-row-count)
			   (< 1 next-smallest-row-count (+ (basis-size b) 1))))
		  ;; scan row count
		  (when (pip-scan-row-counts b)
		    (pip-sort-refs b))
		  (when (< (basis-pip-last b) (basis-pip-first b))
		    (return)))
		 ((= 1 smallest-row-count)
		  ;; spike pivot
		  (pip-spike-pivot b)
		  (when (< (basis-pip-last b) (basis-pip-first b))
		    (return)))
		 (t
		  (error "error in preassigned pivot procedure~%")))))
      ;; return t when done with pivoting
      ;; return nil if singular
      (< (basis-pip-last b) (basis-pip-first b)))))

		    
