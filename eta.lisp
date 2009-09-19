;;;; Data structures for everything regarding eta matrices


(defstruct (lu-eta-matrix
	     (:constructor %make-lu-eta-matrix))
  (j       -1  :type fixnum)
  (col-ref -1  :type fixnum)
  (coef    1   :type rational)
  (is      #() :type vector)
  (vfs     #() :type vector)
  (vis     #() :type vector))



;;;; LU ETA constructor
(defun make-lu-eta-matrix ()
  (%make-lu-eta-matrix
   :is     (make-vector fixnum)
   :vfs    (make-vector float)
   :vis    (make-vector integer)))



;;;; Empties ETA matrix
(defun reset-lu-eta-matrix (eta)
  (setf (lu-eta-matrix-col-ref eta)           -1
	(lu-eta-matrix-coef eta)               1
	(fill-pointer (lu-eta-matrix-is eta))  0
	(fill-pointer (lu-eta-matrix-vis eta)) 0
	(fill-pointer (lu-eta-matrix-vfs eta)) 0))
	
  
(defun lu-eta-get-ratio (eta index)
  (let ((v (aref (lu-eta-matrix-vis eta) index))
	(c (lu-eta-matrix-coef eta)))
    (* c v)))


(defun lu-eta-select-pivot (eta index)
  (unless (zerop index)
    (let ((i  (aref (lu-eta-matrix-is eta) index))
	  (vi (aref (lu-eta-matrix-vis eta) index))
	  (vf (aref (lu-eta-matrix-vfs eta) index)))
      (setf (aref (lu-eta-matrix-is eta)  index) (aref (lu-eta-matrix-is eta)  0)
	    (aref (lu-eta-matrix-vis eta) index) (aref (lu-eta-matrix-vis eta) 0)
	    (aref (lu-eta-matrix-vfs eta) index) (aref (lu-eta-matrix-vfs eta) 0))
      (setf (aref (lu-eta-matrix-is eta)  0) i
	    (aref (lu-eta-matrix-vis eta) 0) vi
	    (aref (lu-eta-matrix-vfs eta) 0) vf))))


(defun lu-eta-remove (eta index)
  (let ((last-i  (vector-pop (lu-eta-matrix-is eta)))
	(last-vi (vector-pop (lu-eta-matrix-vis eta)))
	(last-vf (vector-pop (lu-eta-matrix-vfs eta))))
    (unless (= index (length (lu-eta-matrix-is eta)))
      (setf (aref (lu-eta-matrix-is eta)  index) last-i
	    (aref (lu-eta-matrix-vis eta) index) last-vi
	    (aref (lu-eta-matrix-vfs eta) index) last-vf))))
    

(defun lu-eta-normalize (eta)
  (let* ((vis  (lu-eta-matrix-vis eta))
	 (vfs  (lu-eta-matrix-vfs eta))
	 (n-nz (length vis)))
    (loop 
       (when (dotimes (index n-nz t)
	       (when (zerop (aref vis index))
		 (lu-eta-remove eta index)
		 (decf n-nz)
		 (return nil)))
	 (return)))
    (unless (zerop n-nz)
      (let ((vgcd  0))
	(dotimes (index n-nz)
	  (when (= 1 vgcd)
	    (return))
	  (if (zerop index)
	      (setf vgcd (aref vis 0))
	      (setf vgcd (gcd vgcd (aref vis index)))))
	(unless (= 1 vgcd)
	  (dotimes (index n-nz)
	    (divf (aref vis index) vgcd))
	  (mulf (lu-eta-matrix-coef eta) vgcd)))
      (dotimes (index n-nz)
	(setf (aref vfs index) 
	      (float (* (lu-eta-matrix-coef eta) (aref vis index))))))))

