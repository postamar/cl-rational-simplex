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
    (rotatef (aref (lu-eta-matrix-is eta)  index) (aref (lu-eta-matrix-is eta)  0))
    (rotatef (aref (lu-eta-matrix-vis eta) index) (aref (lu-eta-matrix-vis eta) 0))
    (rotatef (aref (lu-eta-matrix-vfs eta) index) (aref (lu-eta-matrix-vfs eta) 0))))


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
    (unless (zerop n-nz)
      (let ((vgcd 0))
	(dotimes (index n-nz)
	  (let ((vi (aref vis index)))
	    (cond ((= 1 vgcd)
		   (return))
		  ((zerop vi))
		  ((zerop vgcd)
		   (setf vgcd vi))
		  (t
		   (setf vgcd (gcd vgcd vi))))))
	(cond ((zerop vgcd)
	       (setf (lu-eta-matrix-coef eta) 1))
	      ((= 1 vgcd))
	      (t 
	       (dotimes (index n-nz)
		 (divf (aref vis index) vgcd))
	       (mulf (lu-eta-matrix-coef eta) vgcd))))
      (dotimes (index n-nz)
	(setf (aref vfs index) 
	      (float (* (lu-eta-matrix-coef eta) (aref vis index))))))))



(defun lu-eta-sort-indices-increasing (eta)
  (flet ((sift-down (root end)
	   (loop
	      (let* ((child (+ 1 (* 2 root)))
		     (child-i (+ child 1))
		     (root-i (+ root 1)))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (aref (lu-eta-matrix-is eta) child-i) 
			      (aref (lu-eta-matrix-is eta) (+ child-i 1))))
		  (incf child)
		  (incf child-i))
		(when (>= (aref (lu-eta-matrix-is eta) root-i) 
			  (aref (lu-eta-matrix-is eta) child-i))
		  (return))
		(rotatef (aref (lu-eta-matrix-is eta) root-i) 
			 (aref (lu-eta-matrix-is eta) child-i))
		(rotatef (aref (lu-eta-matrix-vis eta) root-i) 
			 (aref (lu-eta-matrix-vis eta) child-i))
		(rotatef (aref (lu-eta-matrix-vfs eta) root-i) 
			 (aref (lu-eta-matrix-vfs eta) child-i))
		(setf root child)))))
    ;; max-heapify
    (let ((count (- (length (lu-eta-matrix-is eta)) 1)))
      (loop for heapify-start from (floor (- count 2) 2) downto 0
	 do (sift-down heapify-start (- count 1)))
      ;; sort
      (loop for end from (- count 1) downto 1
	 do (let ((end-i (+ end 1)))
	      (rotatef (aref (lu-eta-matrix-is eta) 1) 
		       (aref (lu-eta-matrix-is eta) end-i))
	      (rotatef (aref (lu-eta-matrix-vis eta) 1) 
		       (aref (lu-eta-matrix-vis eta) end-i))
	      (rotatef (aref (lu-eta-matrix-vfs eta) 1) 
		       (aref (lu-eta-matrix-vfs eta) end-i))
		   (sift-down 0 (- end 1)))))))
