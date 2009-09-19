;;;; Data structures for everything regarding current basis



(defstruct (basis
	     (:constructor %make-basis))
  (size          0   :type fixnum)
  (is-singular   nil :type symbol)
  (singular-ref  -1  :type fixnum)
  (ppivot-coef   0.0 :type float)
  (refs          #() :type vector)
  (flags         #() :type vector)
  (spikes        #() :type vector)
  (l-columns     #() :type vector)
  (u-columns     #() :type vector)
  (col-hrefs     #() :type vector)
  (col-is        #() :type vector)
  (col-nnz       #() :type vector)
  (row-hrefs     #() :type vector)
  (row-js        #() :type vector)
  (row-cis       #() :type vector)
  (row-nnz       #() :type vector)
  (heap-js       #() :type vector)
  (heap-is       #() :type vector)
  (heap-cis      #() :type vector)
  (href->hi      #() :type vector)
  (hi->href      #() :type vector)
  (i->pi         #() :type vector)
  (pi->i         #() :type vector)
  (j->pj         #() :type vector)
  (pj->j         #() :type vector)
  (lu-ppivots1   #() :type vector)
  (lu-ppivots2   #() :type vector)
  )



;;;; Basis constructor
(defun make-basis (&key (m -1) (lp nil) (ppivot-coef 0.0))
  (when lp
    (setf m (length (lp-active-row-refs lp))))
  (let ((b (%make-basis
	    :size          m
	    :ppivot-coef   ppivot-coef
	    :refs          (make-nvector m -1 fixnum)
	    :flags         (make-nvector m nil boolean)
	    :spikes        (make-nvector m -1 fixnum)
	    :l-columns     (make-nvector m (make-lu-eta-matrix) lu-eta-matrix)
	    :u-columns     (make-vector lu-eta-matrix)
	    :col-nnz       (make-nvector m 0 fixnum)
	    :row-nnz       (make-nvector m 0 fixnum)
	    :col-hrefs     (make-nvector m #() vector)
	    :col-is        (make-nvector m #() vector)
	    :row-hrefs     (make-nvector m #() vector)
	    :row-js        (make-nvector m #() vector)
	    :row-cis       (make-nvector m #() vector)
	    :heap-js       (make-vector fixnum)
	    :heap-is       (make-vector fixnum)
	    :heap-cis      (make-vector fixnum)
	    :href->hi      (make-vector fixnum)
	    :hi->href      (make-vector fixnum)
	    :j->pj         (make-nvector m -1 fixnum)
	    :i->pi         (make-nvector m -1 fixnum)
	    :pj->j         (make-nvector m -1 fixnum)
	    :pi->i         (make-nvector m -1 fixnum)
	    :lu-ppivots1   (make-nvector m -1 fixnum)
	    :lu-ppivots2   (make-nvector m -1 fixnum))))
    (dotimes (k m)
      (setf (aref (basis-col-hrefs b) k) (make-vector fixnum)
	    (aref (basis-col-is b) k)    (make-vector fixnum)
	    (aref (basis-row-js b) k)    (make-vector fixnum)
	    (aref (basis-row-cis b) k)   (make-vector fixnum)
	    (aref (basis-row-hrefs b) k) (make-vector fixnum)
	    (aref (basis-l-columns b) k) (make-lu-eta-matrix)))
    b))



;;;; Resets basis	 
(defun reset-basis (b)
  (let ((m (basis-size b)))
    (setf (basis-is-singular b) nil
	  (basis-singular-ref b) -1
	  (fill-pointer (basis-refs b)) m
	  (fill-pointer (basis-flags b)) m
	  (fill-pointer (basis-spikes b)) m
	  (fill-pointer (basis-col-hrefs b)) m
	  (fill-pointer (basis-col-is b)) m
	  (fill-pointer (basis-col-nnz b)) m
	  (fill-pointer (basis-row-hrefs b)) m
	  (fill-pointer (basis-row-js b)) m
	  (fill-pointer (basis-row-cis b)) m
	  (fill-pointer (basis-row-nnz b)) m
	  (fill-pointer (basis-i->pi b)) m
	  (fill-pointer (basis-j->pj b)) m
	  (fill-pointer (basis-pi->i b)) m
	  (fill-pointer (basis-pj->j b)) m
	  (fill-pointer (basis-href->hi b)) 0
	  (fill-pointer (basis-hi->href b)) 0
	  (fill-pointer (basis-heap-js b)) 0
	  (fill-pointer (basis-heap-is b)) 0
	  (fill-pointer (basis-heap-cis b)) 0
	  (fill-pointer (basis-u-columns b)) 0)
    (when (< (length (basis-refs b)) m)
      (setf (fill-pointer (basis-refs b)) m))
    (dotimes (k m)
      (setf (fill-pointer (aref (basis-row-js b) k)) 0
	    (fill-pointer (aref (basis-row-cis b) k)) 0
	    (fill-pointer (aref (basis-row-hrefs b) k)) 0
	    (fill-pointer (aref (basis-col-hrefs b) k)) 0
	    (fill-pointer (aref (basis-col-is b) k)) 0
	    (aref (basis-refs b) k) k
	    (aref (basis-spikes b) k) m
	    (aref (basis-col-nnz b) k) 0
	    (aref (basis-row-nnz b) k) 0
	    (aref (basis-i->pi b) k) k
	    (aref (basis-pi->i b) k) k
	    (aref (basis-j->pj b) k) k
	    (aref (basis-pj->j b) k) k)
      (let ((l_k (aref (basis-l-columns b) k)))
	(setf (lu-eta-matrix-j l_k) k)
	(reset-lu-eta-matrix l_k)))))



;;;; Singularity declaration
(defun basis-row-is-redundant (b i)
  (setf (basis-is-singular b) 'redundant-row
	(basis-singular-ref b) i))

(defun basis-column-is-redundant (b j)
  (setf (basis-is-singular b) 'redundant-column
	(basis-singular-ref b) j))

	
	      


;;;;; Output functions

(defun print-2d-array (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t "   .   ")
	  (format t "~6,2F " (float (aref a i j)))))))

(defun print-basis-l (b)
  (let* ((m (basis-size b))
	 (j->pj (basis-j->pj b))
	 (i->pi (basis-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m)
      (let* ((l_j  (aref (basis-l-columns b) j))
	     (is   (lu-eta-matrix-is l_j))
	     (vis  (lu-eta-matrix-vis l_j))
	     (n-nz (length is))
	     (f    (lu-eta-matrix-coef l_j)))
	(dotimes (k n-nz)
	  (setf (aref a (aref i->pi (aref is k)) (aref j->pj j)) (* f (aref vis k))))))
    (print-2d-array a)))

(defun print-basis-u (b)
  (let* ((m (basis-size b))
	 (j->pj (basis-j->pj b))
	 (i->pi (basis-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (k m)
      (setf (aref a k k) 1))
    (dotimes (k (length (basis-u-columns b)))
      (let* ((u    (aref (basis-u-columns b) k))
	     (j    (lu-eta-matrix-j u))
	     (is   (lu-eta-matrix-is u))
	     (vis  (lu-eta-matrix-vis u))
	     (n-nz (length is))
	     (f    (lu-eta-matrix-coef u)))
	(dotimes (r n-nz)
	  (setf (aref a (aref i->pi (aref is r)) (aref j->pj j)) (* f (aref vis r))))))
    (print-2d-array a)))
	
(defun print-basis-nz (b)
  (let* ((m (basis-size b))
	 (j->pj (basis-j->pj b))
	 (i->pi (basis-i->pi b))
	 (a (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m)
      (let* ((l_j  (aref (basis-l-columns b) j))
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
      

