;;;; Data structures for everything regarding current basis


(symbol-macrolet
    ((err (error "basis-matrix constructor")))
  (defstruct (basis-matrix
	       (:constructor %make-basis-matrix))
    (size             0   :type fixnum)
    (refactorization-period 1 :type fixnum)
    (is-singular      nil :type symbol)
    (singular-ref     -1  :type fixnum)
    (refs             err :type (simple-array fixnum 1))
    (flags            err :type (simple-array boolean 1))
    (l-file           err :type (simple-array hsv 1))
    (l-pivot-file     err :type (simple-array fixnum 1))
    (lf-file          err :type (simple-array hsv-float 1))
    (n-l-file         0   :type fixnum)
    (n-l-factor-file  0   :type fixnum)
    (u-columns        err :type (simple-array hsv 1))
    (uf-columns       err :type (simple-array hsv-float 1))
    (u-seqs           err :type (simple-array (simple-array fixnum 1) 1))
    (update-row-vals  err :type (simple-array rational 1))
    (fill-ins         err :type (simple-array fixnum 1))
    (col-is           err :type (simple-array (simple-array fixnum 1) 1))
    (col-nnz          err :type (simple-array fixnum 1))
    (row-js           err :type (simple-array (simple-array fixnum 1) 1))
    (row-cis          err :type (simple-array (simple-array fixnum 1) 1)) 
    (row-nnz          err :type (simple-array fixnum 1))
    (col-buckets      err :type (simple-array (simple-array fixnum 1) 1))
    (row-buckets      err :type (simple-array (simple-array fixnum 1) 1))
    (col-bucket-sizes err :type (simple-array fixnum 1))
    (row-bucket-sizes err :type (simple-array fixnum 1))
    (row-col-max      4   :type fixnum)
    (i->pi            err :type (simple-array fixnum 1))
    (pi->i            err :type (simple-array fixnum 1))
    (j->pj            err :type (simple-array fixnum 1))
    (pj->j            err :type (simple-array fixnum 1))))



;;;; Basis matrix constructor
(defun make-basis-matrix (&key 
			  (m -1) 
			  (lp nil) 
			  (max-l-file-size 2000) 
			  (refactorization-period 100)
			  (row-col-max 4))
  (when lp
    (setf m (adjvector-fixnum-fill-pointer (lp-active-row-refs lp))))
  (incf max-l-file-size m)
  (let* ((hsv0 (make-hsv))
	 (hsvf0 (make-hsv-float))
	 (a0 (make-array 1 :initial-element 0 :element-type 'fixnum))
	 (bm (%make-basis-matrix
	    :size             m
	    :refactorization-period refactorization-period
	    :row-col-max      row-col-max
	    :refs             (make-array m :initial-element -1 :element-type 'fixnum)
	    :flags            (make-array m :initial-element nil :element-type 'boolean)
	    :l-file           (make-array max-l-file-size :initial-element hsv0 :element-type 'hsv)
	    :lf-file          (make-array max-l-file-size :initial-element hsvf0 :element-type 'hsv-float)
	    :l-pivot-file     (make-array max-l-file-size :initial-element -1 :element-type 'fixnum)
	    :u-columns        (make-array m :initial-element hsv0 :element-type 'hsv)
	    :uf-columns       (make-array m :initial-element hsvf0 :element-type 'hsv-float)
	    :u-seqs           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :update-row-vals  (make-array m :initial-element 0 :element-type 'rational)
	    :fill-ins         (make-array m :initial-element 0 :element-type 'fixnum)
	    :col-nnz          (make-array m :initial-element 0 :element-type 'fixnum)
	    :row-nnz          (make-array m :initial-element 0 :element-type 'fixnum)
	    :col-is           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-js           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-cis          (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :col-buckets      (make-array (+ m 1) :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-buckets      (make-array (+ m 1) :initial-element a0 :element-type '(simple-array fixnum 1))
	    :col-bucket-sizes (make-array (+ m 1) :initial-element 0 :element-type 'fixnum)
	    :row-bucket-sizes (make-array (+ m 1) :initial-element 0 :element-type 'fixnum)
	    :j->pj            (make-array m :initial-element -1 :element-type 'fixnum)
	    :i->pi            (make-array m :initial-element -1 :element-type 'fixnum)
	    :pj->j            (make-array m :initial-element -1 :element-type 'fixnum)
	    :pi->i            (make-array m :initial-element -1 :element-type 'fixnum))))
    (dotimes (k max-l-file-size)
      (setf (aref (basis-matrix-l-file bm) k)      (make-hsv)
	    (aref (basis-matrix-lf-file bm) k)     (make-hsv-float)))
    (dotimes (k m bm)
      (setf (aref (basis-matrix-u-columns bm) k)   (make-hsv)
	    (aref (basis-matrix-uf-columns bm) k)  (make-hsv-float)
	    (aref (basis-matrix-u-seqs bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-col-is bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-js bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-cis bm) k)     (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-col-buckets bm) (+ k 1)) (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-buckets bm) (+ k 1)) (make-array m :initial-element -1 :element-type 'fixnum)))))



;;;; Resets basis	 
(defun reset-basis-matrix (bm)
  (dotimes (k (basis-matrix-n-l-file bm))
    (reset-hsv-float (aref (basis-matrix-lf-file bm) k))
    (reset-hsv (aref (basis-matrix-l-file bm) k)))
  (setf (basis-matrix-is-singular bm) nil
	(basis-matrix-singular-ref bm) -1
	(basis-matrix-n-l-factor-file bm) 0
	(basis-matrix-n-l-file bm) 0)
  (let ((m (basis-matrix-size bm)))
    (dotimes (k m)
      (setf (aref (basis-matrix-fill-ins bm) k) 0
	    (aref (basis-matrix-refs bm) k) k
	    (aref (basis-matrix-col-nnz bm) k) 0
	    (aref (basis-matrix-row-nnz bm) k) 0
	    (aref (basis-matrix-col-bucket-sizes bm) (+ k 1)) 0
	    (aref (basis-matrix-row-bucket-sizes bm) (+ k 1)) 0
	    (aref (basis-matrix-i->pi bm) k) k
	    (aref (basis-matrix-pi->i bm) k) k
	    (aref (basis-matrix-j->pj bm) k) k
	    (aref (basis-matrix-pj->j bm) k) k)
      (let ((u_k (aref (basis-matrix-u-columns bm) k))
	    (uf_k (aref (basis-matrix-uf-columns bm) k)))
	(reset-hsv u_k)
	(reset-hsv-float uf_k)))))



;;;; Singularity declaration
(defun basis-matrix-row-is-redundant (bm i)
  (setf (basis-matrix-is-singular bm) 'redundant-row
	(basis-matrix-singular-ref bm) i))

(defun basis-matrix-column-is-redundant (bm j)
  (setf (basis-matrix-is-singular bm) 'redundant-column
	(basis-matrix-singular-ref bm) j))

	

;;;; Fill basis according to basis header
(defun fill-basis-matrix (bm lp basis-header)
  (assert (= (basis-matrix-size bm) (length basis-header)))
  (reset-basis-matrix bm)
  (dotimes (j (basis-matrix-size bm))
    (let ((col-ref (aref basis-header j))
	  (u_j     (aref (basis-matrix-u-columns bm) j)))
      (let* ((col          (adjvector-column-ref (lp-columns lp) col-ref))
	     (col-row-refs (hsv-is (column-hsv col)))
	     (col-vals     (hsv-vis (column-hsv col)))
	     (n-nz         (hsv-length (column-hsv col))))
	(setf (hsv-coef u_j) (hsv-coef (column-hsv col)))
	(dotimes (k n-nz)
	  (let ((i (adjvector-fixnum-ref (lp-active-row-inds lp)  (aref col-row-refs k)))
		(ci (hsv-length u_j))
		(val (aref col-vals k)))
	    (unless (= -1 i)
	      (hsv-add i val u_j)
	      (setf (aref (aref (basis-matrix-col-is bm) j) 
			  (aref (basis-matrix-col-nnz bm) j)) i)
	      (setf (aref (aref (basis-matrix-row-js bm) i)
			  (aref (basis-matrix-row-nnz bm) i)) j)
	      (setf (aref (aref (basis-matrix-row-cis bm) i)
			  (aref (basis-matrix-row-nnz bm) i)) ci)
	      (incf (aref (basis-matrix-col-nnz bm) j))
	      (incf (aref (basis-matrix-row-nnz bm) i))))))))
  (dotimes (k (basis-matrix-size bm))
    (let ((nnz-i (aref (basis-matrix-row-nnz bm) k))
	  (nnz-j (aref (basis-matrix-col-nnz bm) k)))
      (cond ((zerop nnz-j)
	     (basis-matrix-column-is-redundant bm k)
	     (return))
	    ((zerop nnz-i)
	     (basis-matrix-row-is-redundant bm k)
	     (return))
	    (t
	     (setf (aref (aref (basis-matrix-col-buckets bm) nnz-j)
			 (aref (basis-matrix-col-bucket-sizes bm) nnz-j))
		   k)
	     (incf (aref (basis-matrix-col-bucket-sizes bm) nnz-j))
	     (setf (aref (aref (basis-matrix-row-buckets bm) nnz-i)
			 (aref (basis-matrix-row-bucket-sizes bm) nnz-i))
		   k) 
	     (incf (aref (basis-matrix-row-bucket-sizes bm) nnz-i))))))
  (not (basis-matrix-is-singular bm)))



;;;;; Output functions

(defun print-2d-array (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t "   .   ")
	  (format t "~6,2F " (float (aref a i j)))))))

(defun print-2d-array-nz (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t " .")
	  (format t " x")))))

