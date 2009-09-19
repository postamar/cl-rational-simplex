;;;; Sparse square matrix implementation
;;;;



;;;; Data structure for a sparse square matrix of rational numbers
;;;; nothing is ordered
(defstruct (sparse-square-matrix
	     (:constructor %make-sparse-square-matrix))
  (size           0   :type fixnum)
  (values         #() :type vector)
  (ref-buffer     #() :type vector)
  (flag-buffer    #() :type vector)
  (row-avail      #() :type vector)
  (col-avail      #() :type vector)
  (ppivot-factor  0.0 :type float)
  (is-singular    nil :type symbol)
  (singular-ref   0   :type fixnum) 
  (spikes         #() :type vector)
  (col-perm-ref   #() :type vector)
  (ref-perm-col   #() :type vector)
  (row-perm-ref   #() :type vector)
  (ref-perm-row   #() :type vector)
  (col-row-ref    #() :type vector)
  (col-val-ref    #() :type vector)
  (row-val-ref    #() :type vector)
  (row-col-ref    #() :type vector)
  (row-count      #() :type vector)
  (col-count      #() :type vector))



;;;; Constructor
(defun make-sparse-square-matrix (m &key (ppivot-factor 0.0))
  (let ((ssqm (%make-sparse-square-matrix 
	       :size          m
	       :values        (make-vector rational)
	       :ppivot-factor ppivot-factor 
	       :ref-buffer    (make-array m :initial-element 0 :element-type 'fixnum)
	       :flag-buffer   (make-array m :initial-element nil :element-type 'boolean)
	       :row-avail     (make-array m :initial-element nil :element-type 'boolean)
	       :col-avail     (make-array m :initial-element nil :element-type 'boolean)
	       :col-perm-ref  (make-array m :initial-element 0 :element-type 'fixnum)
	       :ref-perm-col  (make-array m :initial-element 0 :element-type 'fixnum)
	       :row-perm-ref  (make-array m :initial-element 0 :element-type 'fixnum)
	       :ref-perm-row  (make-array m :initial-element 0 :element-type 'fixnum)
	       :col-row-ref   (make-array m :initial-element #() :element-type 'vector)
	       :col-val-ref   (make-array m :initial-element #() :element-type 'vector)
	       :row-val-ref   (make-array m :initial-element #() :element-type 'vector)
	       :row-col-ref   (make-array m :initial-element #() :element-type 'vector)
	       :spikes        (make-array m :initial-element 0 :element-type 'fixnum)
	       :row-count     (make-array m :initial-element 0 :element-type 'fixnum)
	       :col-count     (make-array m :initial-element 0 :element-type 'fixnum))))
    (dotimes (k m ssqm)
      (setf (aref (sparse-square-matrix-row-val-ref ssqm) k) (make-vector fixnum)
	    (aref (sparse-square-matrix-row-col-ref ssqm) k) (make-vector fixnum)
	    (aref (sparse-square-matrix-col-row-ref ssqm) k) (make-vector fixnum)
	    (aref (sparse-square-matrix-col-val-ref ssqm) k) (make-vector fixnum)
	    (aref (sparse-square-matrix-ref-perm-row ssqm) k) k
	    (aref (sparse-square-matrix-row-perm-ref ssqm) k) k
	    (aref (sparse-square-matrix-ref-perm-col ssqm) k) k
	    (aref (sparse-square-matrix-col-perm-ref ssqm) k) k))))
    


;;;; Fills existing matrix with zeros
(defun reset-matrix (ssqm)
  (let ((m (sparse-square-matrix-size ssqm)))
    (setf (fill-pointer (sparse-square-matrix-values ssqm)) 0)
    (dotimes (k m)
      (setf (fill-pointer (aref (sparse-square-matrix-col-row-ref ssqm) k)) 0
	    (fill-pointer (aref (sparse-square-matrix-col-val-ref ssqm) k)) 0
	    (fill-pointer (aref (sparse-square-matrix-row-col-ref ssqm) k)) 0
	    (fill-pointer (aref (sparse-square-matrix-row-val-ref ssqm) k)) 0
	    (sparse-square-matrix-is-singular ssqm) nil
	    (aref (sparse-square-matrix-ref-perm-row ssqm) k) k
	    (aref (sparse-square-matrix-row-perm-ref ssqm) k) k
	    (aref (sparse-square-matrix-ref-perm-col ssqm) k) k
	    (aref (sparse-square-matrix-col-perm-ref ssqm) k) k))))



;;;; Fills matrix using sparse columns
(defun fill-matrix (ssqm scol-indices scol-values)
  (let ((m (sparse-square-matrix-size ssqm)))
    (assert (= m (length scol-indices)))
    (assert (= m (length scol-values)))
    (reset-matrix ssqm)
    (dotimes (j m)
      (let* ((ind_j (aref scol-indices j))
	     (val_j (aref scol-values j))
	     (col-size (length ind_j)))
	(assert (= col-size (length val_j)))
	(symbol-macrolet
	    ((col-row-ref (aref (sparse-square-matrix-col-row-ref ssqm) j))
	     (col-val-ref (aref (sparse-square-matrix-col-val-ref ssqm) j)))
	  (dotimes (k col-size)
	    (let ((val-ref (length (sparse-square-matrix-values ssqm)))
		  (val (aref val_j k))
		  (i (aref ind_j k)))
	      (symbol-macrolet
		  ((row-col-ref (aref (sparse-square-matrix-row-col-ref ssqm) i))
		   (row-val-ref (aref (sparse-square-matrix-row-val-ref ssqm) i)))
		(unless (zerop val)
		  (vector-push-extend val (sparse-square-matrix-values ssqm))
		  (vector-push-extend i col-row-ref)
		  (vector-push-extend val-ref col-val-ref)
		  (vector-push-extend j row-col-ref)
		  (vector-push-extend val-ref row-val-ref))))))))))




;;;; Permutates two columns
(defun switch-columns-in-matrix (ssqm j1 j2)
  (let* ((cpr (sparse-square-matrix-col-perm-ref ssqm))
	 (rpc (sparse-square-matrix-ref-perm-col ssqm))
	 (ref1 (aref cpr j1))
	 (ref2 (aref cpr j2)))
    (setf (aref cpr j1) (aref cpr j2))
    (setf (aref cpr j2) ref1)
    (setf (aref rpc ref1) j2
	  (aref rpc ref2) j1)))



;;;; Permutates two rows
(defun switch-rows-in-matrix (ssqm i1 i2)
  (let* ((ipr (sparse-square-matrix-row-perm-ref ssqm))
	 (rpi (sparse-square-matrix-ref-perm-row ssqm))
	 (ref1 (aref ipr i1))
	 (ref2 (aref ipr i2)))
    (setf (aref ipr i1) (aref ipr i2))
    (setf (aref ipr i2) ref1)
    (setf (aref rpi ref1) i2
	  (aref rpi ref2) i1)))


	       

;;;; Permutates column k with column p >= k if |p| > c |k|
;;;; with c being the partial pivoting factor   
;;;; returns -1 is matrix is singular
(defun decomposition-partial-pivoting (ssqm k)
  (let ((row-ref (aref (sparse-square-matrix-row-perm-ref ssqm) k))
	(values (sparse-square-matrix-values ssqm))
	(ref-perm-col (sparse-square-matrix-ref-perm-col ssqm))
	(c (sparse-square-matrix-ppivot-factor ssqm))
	(current-abs-val 0)
	(largest-abs-val 0)
	(p -1))
    (let ((row-col-ref (aref (sparse-square-matrix-row-col-ref ssqm) row-ref)) 
	  (row-val-ref (aref (sparse-square-matrix-row-val-ref ssqm) row-ref)))
      (dotimes (col-ref-index (length row-col-ref))
	(let* ((col-ref (aref row-col-ref col-ref-index))
	       (abs-val (abs (aref values (aref row-val-ref col-ref-index))))
	       (j (aref ref-perm-col col-ref)))
	  (when (= k j)
	    (setf current-abs-val abs-val))
	  (when (and (<= k j)
		     (not (zerop abs-val))
		     (or (= -1 p) 
			 (< largest-abs-val abs-val)))
	    (setf p j
		  largest-abs-val abs-val))))
      (unless (= -1 p)
	(if (or (= k p)
		(= 0.0 c)
		(> (* current-abs-val c) largest-abs-val))
	    (setf p k)
	    (switch-columns-in-matrix ssqm k p)))
      p)))

  
;;;; LU decomposition
;;;; returns t on success, nil if ssqm is singular
(defun decompose-matrix (ssqm p-file)
  (let ((m (sparse-square-matrix-size ssqm))
	(values (sparse-square-matrix-values ssqm))
	(ref-buffer (sparse-square-matrix-ref-buffer ssqm))
	(flag-buffer (sparse-square-matrix-flag-buffer ssqm))
	(col-perm-ref (sparse-square-matrix-col-perm-ref ssqm))
	(row-perm-ref (sparse-square-matrix-row-perm-ref ssqm))
	(ref-perm-row (sparse-square-matrix-ref-perm-row ssqm))
	(row-col-refs (sparse-square-matrix-row-col-ref ssqm))
	(row-val-refs (sparse-square-matrix-row-val-ref ssqm))
	(col-row-refs (sparse-square-matrix-col-row-ref ssqm))
	(col-val-refs (sparse-square-matrix-col-val-ref ssqm)))
    (dotimes (k m t)
      ;; permutate matrix, check for singularity
      (when (= -1 (setf (aref p-file k) 
			(decomposition-partial-pivoting ssqm k)))
	(return nil))
      ;; pivot matrix
      (let* ((col-ref (aref col-perm-ref k))
	     (col-row-ref (aref col-row-refs col-ref))
	     (col-val-ref (aref col-val-refs col-ref))
	     (pivot-val-ref -1))
	(dotimes (aux-i (- m k))
	  (setf (aref ref-buffer (+ k aux-i)) -1))
	;; process column k
	(dotimes (index (length col-row-ref))
	  (when (= k (aref ref-perm-row (aref col-row-ref index)))
	    (setf pivot-val-ref (aref col-val-ref index))
	    (setf (aref ref-buffer k) pivot-val-ref)
	    (return)))
	(assert (/= -1 pivot-val-ref))
	(setf (aref values pivot-val-ref) (/ 1 (aref values pivot-val-ref)))
	(dotimes (index (length col-row-ref))
	  (let ((i (aref ref-perm-row (aref col-row-ref index))))
	    (when (< k i)
	      (let ((val-ref (aref col-val-ref index)))
		(setf (aref ref-buffer i) val-ref)
		(mulf (aref values val-ref) (- (aref values pivot-val-ref)))))))
	;; process columns j > k
	(dotimes (aux-j (- m 1 k))
	  (let* ((val-ref-kj -1)
		 (j (+ k 1 aux-j))
		 (col-ref (aref col-perm-ref j))
		 (col-row-ref (aref col-row-refs col-ref))
		 (col-val-ref (aref col-val-refs col-ref)))
	    (dotimes (aux-i (- m 1 k))
		(setf (aref flag-buffer (+ k 1 aux-i)) nil))
	    (dotimes (index (length col-row-ref))
	      (when (= k (aref ref-perm-row (aref col-row-ref index)))
		(setf val-ref-kj (aref col-val-ref index))
		(return)))
	    (unless (= -1 val-ref-kj)
	      ;; on rows i > k with nonzeros
	      (dotimes (index (length col-row-ref))
		(let ((i (aref ref-perm-row (aref col-row-ref index))))
		  (when (< k i)
		    (setf (aref flag-buffer i) t)
		    (let ((val-ref-ij (aref col-val-ref index))
			  (val-ref-ik (aref ref-buffer i)))
		      (unless (= -1 val-ref-ik)
			(incf (aref values val-ref-ij)
			      (* (aref values val-ref-ik)
				 (aref values val-ref-kj))))))))
	      ;; on rows i > k with zeros
	      (dotimes (aux-i (- m 1 k))
		(let* ((i (+ k 1 aux-i))
		       (val-ref-ik (aref ref-buffer i)))
		  (unless (or (= -1 val-ref-ik)
			      (aref flag-buffer i))
		    ;; add new cell
		    (let ((val (* (aref values val-ref-ik) (aref values val-ref-kj))))
		      (unless (zerop val)
			(let ((row-ref (aref row-perm-ref i))
			      (val-ref-ij (vector-push-extend val values)))
			  (vector-push-extend row-ref    (aref col-row-refs col-ref))
			  (vector-push-extend val-ref-ij (aref col-val-refs col-ref))
			  (vector-push-extend col-ref    (aref row-col-refs row-ref))
			  (vector-push-extend val-ref-ij (aref row-val-refs row-ref))))))))
	      ;; on row k
	      (mulf (aref values val-ref-kj) (aref values pivot-val-ref)))))))))



;;;;      
(defun sparse-square-matrix-to-array (ssqm)
  (let* ((m (sparse-square-matrix-size ssqm))
	 (ref-perm-row (sparse-square-matrix-ref-perm-row ssqm))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m a)
      (let* ((ref-col (aref (sparse-square-matrix-col-perm-ref ssqm) j))
	     (col-row-ref (aref (sparse-square-matrix-col-row-ref ssqm) ref-col))
	     (col-val-ref (aref (sparse-square-matrix-col-val-ref ssqm) ref-col)))
	(dotimes (index (length col-row-ref))
	  (setf (aref a (aref ref-perm-row (aref col-row-ref index)) j)
		(aref (sparse-square-matrix-values ssqm) (aref col-val-ref index))))))))



;;;;
(defun print-2d-array (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t "   .   ")
	  (format t "~6,2F " (float (aref a i j)))))))



;;;;
(defun print-sparse-square-matrix (ssqm)
  (print-2d-array (sparse-square-matrix-to-array ssqm)))



;;;;
(defun check-LU (original-ssqm decomp-ssqm p-file)
  (let ((oa (sparse-square-matrix-to-array original-ssqm))
	(da (sparse-square-matrix-to-array decomp-ssqm))
	(ref-perm-col (sparse-square-matrix-ref-perm-col decomp-ssqm))
	(ref-perm-row (sparse-square-matrix-ref-perm-row decomp-ssqm))
	(m (sparse-square-matrix-size decomp-ssqm)))
    (let ((ca (make-array (list m m) :initial-element 0 :element-type 'rational)))
      (dotimes (i m)
	(dotimes (j m)
	  (setf (aref ca (aref ref-perm-row i) (aref ref-perm-col j)) 
		(aref oa i j))))
      ;; go ahead
      (dotimes (k m)
	;; permutate columns
	(let ((j (aref p-file k)))
	  (when (/= j k)
	    (dotimes (i m)
	      (rotatef (aref ca i j) (aref ca i k)))))
	;; left-multiply by pivot
	;; rows i > k
	(dotimes (aux-i (- m 1 k))
	  (let ((i (+ k 1 aux-i)))
	    (dotimes (j m)
	      (incf (aref ca i j)
		    (* (aref da i k)
		       (aref ca k j))))))
	;; row k
	(dotimes (j m)
	  (mulf (aref ca k j) (aref da k k))))
      ;; check for equality
      (dotimes (i m t)
	(dotimes (j m)
	  (cond ((= i j)
		 (assert (= (aref ca i j) 1)))
		((< i j)
		 (assert (= (aref ca i j) (aref da i j))))))))))
	    


	    
	    
     
    
