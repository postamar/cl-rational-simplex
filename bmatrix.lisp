;;;; Data structures for everything regarding current basis

(splay-tree :name "pivot-bucket"
	    :key-type fixnum
	    :val-type fixnum)

(splay-tree :name "markowitz-tree"
	    :key-type fixnum 
	    :val-type pivot-bucket)

(defstruct (basis-matrix
	     (:constructor %make-basis-matrix))
  (size          0   :type fixnum)
  (is-singular   nil :type symbol)
  (singular-ref  -1  :type fixnum)
  (refs          #() :type vector)
  (flags         #() :type vector)
  (spikes        #() :type vector)
  (l-columns     #() :type vector)
  (u-columns     #() :type vector)
  (u-file        #() :type vector)
  (col-is        #() :type vector)
  (col-nnz       #() :type vector)
  (row-js        #() :type vector)
  (row-cis       #() :type vector)
  (row-nnz       #() :type vector)
  (pivot-buckets (error "basis constructor")
		 :type markowitz-tree)
  (pivot-i-queue #() :type vector)
  (pivot-j-queue #() :type vector)
  (pivot-ci-queue #() :type vector)
  (i->pi         #() :type vector)
  (pi->i         #() :type vector)
  (j->pj         #() :type vector)
  (pj->j         #() :type vector))



;;;; Basis matrix constructor
(defun make-basis-matrix (&key (m -1) (lp nil))
  (when lp
    (setf m (length (lp-active-row-refs lp))))
  (let ((b (%make-basis-matrix
	    :size          m
	    :refs          (make-nvector m -1 fixnum)
	    :flags         (make-nvector m nil boolean)
	    :spikes        (make-nvector m -1 fixnum)
	    :l-columns     (make-nvector m (make-lu-eta-matrix) lu-eta-matrix)
	    :u-columns     (make-vector lu-eta-matrix)
	    :u-file        (make-vector lu-eta-matrix)
	    :col-nnz       (make-nvector m 0 fixnum)
	    :row-nnz       (make-nvector m 0 fixnum)
	    :col-is        (make-nvector m #() vector)
	    :row-js        (make-nvector m #() vector)
	    :row-cis       (make-nvector m #() vector)
	    :pivot-buckets (make-markowitz-tree)
	    :pivot-i-queue (make-vector fixnum)
	    :pivot-j-queue (make-vector fixnum)
	    :pivot-ci-queue (make-vector fixnum)
	    :j->pj         (make-nvector m -1 fixnum)
	    :i->pi         (make-nvector m -1 fixnum)
	    :pj->j         (make-nvector m -1 fixnum)
	    :pi->i         (make-nvector m -1 fixnum))))
    (dotimes (k m)
      (setf (aref (basis-matrix-l-columns b) k) (make-lu-eta-matrix)
	    (aref (basis-matrix-col-is b) k)    (make-vector fixnum)
	    (aref (basis-matrix-row-js b) k)    (make-vector fixnum)
	    (aref (basis-matrix-row-cis b) k)   (make-vector fixnum)))
    b))



;;;; Resets basis	 
(defun reset-basis-matrix (b)
  (let ((m (basis-matrix-size b)))
    (setf (basis-matrix-is-singular b) nil
	  (basis-matrix-singular-ref b) -1
	  (fill-pointer (basis-matrix-refs b)) m
	  (fill-pointer (basis-matrix-flags b)) m
	  (fill-pointer (basis-matrix-spikes b)) m
	  (fill-pointer (basis-matrix-col-is b)) m
	  (fill-pointer (basis-matrix-col-nnz b)) m
	  (fill-pointer (basis-matrix-row-js b)) m
	  (fill-pointer (basis-matrix-row-cis b)) m
	  (fill-pointer (basis-matrix-row-nnz b)) m
	  (fill-pointer (basis-matrix-i->pi b)) m
	  (fill-pointer (basis-matrix-j->pj b)) m
	  (fill-pointer (basis-matrix-pi->i b)) m
	  (fill-pointer (basis-matrix-pj->j b)) m
	  (fill-pointer (basis-matrix-u-file b)) 0
	  (fill-pointer (basis-matrix-u-columns b)) 0)
    (when (< (length (basis-matrix-refs b)) m)
      (setf (fill-pointer (basis-matrix-refs b)) m))
    (map-markowitz-tree #'(lambda (mc bucket)
			    (declare (ignore mc))
			    (reset-pivot-bucket bucket))
			(basis-matrix-pivot-buckets b))
    (dotimes (k m)
      (setf (fill-pointer (aref (basis-matrix-row-js b) k)) 0
	    (fill-pointer (aref (basis-matrix-row-cis b) k)) 0
	    (fill-pointer (aref (basis-matrix-col-is b) k)) 0
	    (aref (basis-matrix-refs b) k) k
	    (aref (basis-matrix-spikes b) k) m
	    (aref (basis-matrix-col-nnz b) k) 0
	    (aref (basis-matrix-row-nnz b) k) 0
	    (aref (basis-matrix-i->pi b) k) k
	    (aref (basis-matrix-pi->i b) k) k
	    (aref (basis-matrix-j->pj b) k) k
	    (aref (basis-matrix-pj->j b) k) k)
      (let ((l_k (aref (basis-matrix-l-columns b) k)))
	(setf (lu-eta-matrix-j l_k) k)
	(reset-lu-eta-matrix l_k)))))



;;;; Singularity declaration
(defun basis-matrix-row-is-redundant (b i)
  (setf (basis-matrix-is-singular b) 'redundant-row
	(basis-matrix-singular-ref b) i))

(defun basis-matrix-column-is-redundant (b j)
  (setf (basis-matrix-is-singular b) 'redundant-column
	(basis-matrix-singular-ref b) j))

	

;;;;; Output functions

(defun print-2d-array (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t "   .   ")
	  (format t "~6,2F " (float (aref a i j)))))))

(defun print-basis-matrix-l (b)
  (let* ((m (basis-matrix-size b))
	 (j->pj (basis-matrix-j->pj b))
	 (i->pi (basis-matrix-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m)
      (let* ((l_j  (aref (basis-matrix-l-columns b) j))
	     (is   (lu-eta-matrix-is l_j))
	     (vis  (lu-eta-matrix-vis l_j))
	     (n-nz (length is))
	     (f    (lu-eta-matrix-coef l_j)))
	(dotimes (k n-nz)
	  (setf (aref a (aref i->pi (aref is k)) (aref j->pj j)) (* f (aref vis k))))))
    (print-2d-array a)))

(defun print-basis-matrix-u (b)
  (let* ((m (basis-matrix-size b))
	 (j->pj (basis-matrix-j->pj b))
	 (i->pi (basis-matrix-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (k m)
      (setf (aref a k k) 1))
    (dotimes (k (length (basis-matrix-u-columns b)))
      (let* ((u    (aref (basis-matrix-u-columns b) k))
	     (j    (lu-eta-matrix-j u))
	     (is   (lu-eta-matrix-is u))
	     (vis  (lu-eta-matrix-vis u))
	     (n-nz (length is))
	     (f    (lu-eta-matrix-coef u)))
	(dotimes (r n-nz)
	  (setf (aref a (aref i->pi (aref is r)) (aref j->pj j)) (* f (aref vis r))))))
    (print-2d-array a)))
	
(defun print-basis-nz (b)
  (let* ((m (basis-matrix-size b))
	 (j->pj (basis-matrix-j->pj b))
	 (i->pi (basis-matrix-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m)
      (let* ((l_j  (aref (basis-matrix-l-columns b) j))
	     (is   (lu-eta-matrix-is l_j))
	     (vis  (lu-eta-matrix-vis l_j))
	     (n-nz (length is))
	     (f    (lu-eta-matrix-coef l_j)))
	(dotimes (k n-nz)
	  (setf (aref a (aref i->pi (aref is k)) (aref j->pj j)) (* f (aref vis k))))))
    (dotimes (i (array-dimension a 0) (format t "~%"))
      (dotimes (j (array-dimension a 1) (format t "~%"))
	(if (zerop (aref a i j))
	    (format t " .")
	    (format t " x"))))))
      

