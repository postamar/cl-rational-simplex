;;;; Data structures for everything regarding current basis



(defstruct (basis
	     (:constructor %make-basis))
  ;; Basis state information
  (size         0   :type fixnum)
  (is-singular  nil :type symbol)
  (singular-ref -1  :type fixnum)
  ;; LU and ETA factorization parameters
  (ppivot-coef  0.0 :type float)
  ;; Common data
  (refs         #() :type vector)
  (flags        #() :type vector)
  ;; Basis matrix data
  (row-js       #() :type vector)
  (row-col-inds #() :type vector)
  (spikes       #() :type vector)
  ;; Preassigned pivot data
  (i->pi        #() :type vector)
  (pi->i        #() :type vector)
  (j->pj        #() :type vector)
  (pj->j        #() :type vector)
  (col-counts   #() :type vector)
  (row-counts   #() :type vector)
  (col-avails   #() :type vector)
  (row-avails   #() :type vector)
  (pip-first    -1  :type fixnum)
  (pip-last     -1  :type fixnum)
  (pip-spikes   #() :type vector)
  ;; LU data
  (lu-blocks    #() :type vector)
  (l-columns    #() :type vector)
  (u-columns    #() :type vector)
  (lu-ppivots   #() :type vector)
  ;; TODO
  ;; ETA file
  )



;;;; Basis constructor
(defun make-basis (&key (m -1) (lp nil) (ppivot-coef 0.0))
  (when lp
    (setf m (length (lp-active-row-refs lp))))
  (let* ((basis 
	  (%make-basis
	   :size         m
	   :ppivot-coef  ppivot-coef
	   :j->pj        (make-nvector m -1 fixnum)
	   :i->pi        (make-nvector m -1 fixnum)
	   :pj->j        (make-nvector m -1 fixnum)
	   :pi->i        (make-nvector m -1 fixnum)
	   :row-js       (make-nvector m #() vector)
	   :row-col-inds (make-nvector m #() vector)
	   :refs         (make-nvector m -1 fixnum)
	   :flags        (make-nvector m nil boolean)
	   :col-counts   (make-nvector m -1 fixnum)
	   :row-counts   (make-nvector m -1 fixnum)
	   :col-avails   (make-nvector m t boolean)
	   :row-avails   (make-nvector m t boolean)
	   :pip-spikes   (make-vector fixnum)
	   :spikes       (make-nvector m -1 fixnum)
	   :lu-blocks    (make-vector fixnum)
	   :l-columns    (make-nvector m (make-lu-eta-matrix) lu-eta-matrix)
	   :u-columns    (make-vector lu-eta-matrix)
	   :lu-ppivots   (make-nvector m -1 fixnum))))
    (dotimes (k m)
      (setf (aref (basis-row-js basis) k)       (make-vector fixnum)
	    (aref (basis-row-col-inds basis) k) (make-vector fixnum)
	    (aref (basis-l-columns basis) k)    (make-lu-eta-matrix)))
    basis))



;;;; Resets basis	 
(defun reset-basis (b)
  (setf (basis-is-singular b) nil
	(basis-singular-ref b) -1
	(basis-pip-first b)     0
	(basis-pip-last b) (basis-size b)
	(fill-pointer (basis-lu-blocks b)) 0
	(fill-pointer (basis-pip-spikes b)) 0
	(fill-pointer (basis-u-columns b)) 0)
  (when (< (length (basis-refs b)) (basis-size b))
    (setf (fill-pointer (basis-refs b)) (basis-size b)))
  (dotimes (k (basis-size b))
    (setf (fill-pointer (aref (basis-row-js b) k)) 0
	  (fill-pointer (aref (basis-row-col-inds b) k)) 0
	  (aref (basis-refs b) k) k
	  (aref (basis-spikes b) k) -1
	  (aref (basis-i->pi b) k) k
	  (aref (basis-pi->i b) k) k
	  (aref (basis-j->pj b) k) k
	  (aref (basis-pj->j b) k) k
	  (aref (basis-col-counts b) k) 0
	  (aref (basis-row-counts b) k) 0
	  (aref (basis-col-avails b) k) t
	  (aref (basis-row-avails b) k) t)
    (let ((l_k (aref (basis-l-columns b) k)))
      (setf (lu-eta-matrix-j l_k) k)
      (reset-lu-eta-matrix l_k))))


	     
;;;; Fills basis according to basis header
(defun fill-basis (b lp basis-header)
  (let ((m (basis-size b))
	(actrinds (lp-active-row-inds lp))
	(row-js (basis-row-js b))
	(row-col-inds (basis-row-col-inds b))
	(row-counts (basis-row-counts b))
	(col-counts (basis-col-counts b)))
    (assert (= m (length basis-header)))
    (reset-basis b)
    (dotimes (j m)
      (let ((col-ref (aref basis-header j))
	    (l_j     (aref (basis-l-columns b) j)))
	(let* ((col          (aref (lp-columns lp) col-ref))
	       (col-row-refs (column-row-refs col))
	       (col-vals     (column-values col))
	       (n-nz         (length col-row-refs)))
	  (setf (lu-eta-matrix-coef l_j)    (column-coef col)
		(lu-eta-matrix-col-ref l_j) col-ref
		(aref col-counts j)         n-nz)
	  (dotimes (k n-nz)
	    (let ((row-ref (aref col-row-refs k))
		  (val     (aref col-vals k)))
	      (let ((i (aref actrinds row-ref)))
		(unless (= -1 i)
		  (vector-push-extend i   (lu-eta-matrix-is  l_j))
		  (vector-push-extend val (lu-eta-matrix-vis l_j))
		  (vector-push-extend 0.0 (lu-eta-matrix-vfs l_j))
		  (vector-push-extend j   (aref row-js i))
		  (vector-push-extend k   (aref row-col-inds i)))))))))
    (dotimes (i m)
      (setf (aref row-counts i) (length (aref row-js i))))))
	    



;;;;; Permutations

(defun basis-permutate-columns (b j1 j2)
  (let ((pj1 (aref (basis-j->pj b) j1))
	(pj2 (aref (basis-j->pj b) j2)))
    (rotatef (aref (basis-j->pj b) j1) (aref (basis-j->pj b) j2))
    (setf (aref (basis-pj->j b) pj1) j2
	  (aref (basis-pj->j b) pj2) j1)))


(defun basis-permutate-rows (b i1 i2)
  (let ((pi1 (aref (basis-i->pi b) i1))
	(pi2 (aref (basis-i->pi b) i2)))
    (rotatef (aref (basis-i->pi b) i1) (aref (basis-i->pi b) i2))
    (setf (aref (basis-pi->i b) pi1) i2
	  (aref (basis-pi->i b) pi2) i1)))


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
	    (format t ".")
	    (format t "X"))))))



      

