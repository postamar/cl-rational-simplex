;;;; Data structures for everything regarding hyper-sparse vectors


;;;; hyper-sparse vector data structures
(defstruct (hsv
	     (:constructor %make-hsv))
  (j       -1  :type fixnum)
  (coef    1   :type rational)
  (is      #() :type vector)
  (vis     #() :type vector))

(defstruct (hsv-float
	     (:constructor %make-hsv-float))
  (j       -1  :type fixnum)
  (is      #() :type vector)
  (vfs     #() :type vector))


;;;; HSV constructors
(defun make-hsv ()
  (%make-hsv
   :is  (make-vector fixnum)
   :vis (make-vector integer)))


(defun hsv-length (hsv)
  (length (hsv-is hsv)))


(defun make-hsv-float ()
  (%make-hsv-float
   :is  (make-vector fixnum)
   :vfs (make-vector float)))
		  

(defun copy-hsv-into-hsv-float (hsv hsvf)
  (let* ((n (hsv-length hsv))
	 (c (hsv-coef hsv)))
    (reset-hsv-float hsvf)
    (setf (hsv-float-j hsvf) (hsv-j hsv))
    (dotimes (ci n)
      (vector-push-extend (aref (hsv-is hsv) ci) (hsv-float-is hsvf))
      (vector-push-extend (float (* c (aref (hsv-vis hsv) ci))) (hsv-float-vfs hsvf)))))
    

(defun hsv-float-length (hsvf)
  (length (hsv-float-is hsvf)))


;;;; Empties HSV 
(defun reset-hsv (v)
  (setf (hsv-j v) -1
	(hsv-coef v) 1
	(fill-pointer (hsv-is v))  0
	(fill-pointer (hsv-vis v)) 0))

(defun reset-hsv-float (vf)
  (setf (hsv-float-j vf) -1
	(fill-pointer (hsv-float-is vf))  0
	(fill-pointer (hsv-float-vfs vf)) 0))
	

(defun hsv-ratio (v ci)
  (* (hsv-coef v) (aref (hsv-vis v) ci)))


(defun hsv-normalize (v)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (declare (hsv v))
  (let* ((vis  (hsv-vis v))
	 (n-nz (length vis)))
    (declare ((array fixnum 1) vis))
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
	       (setf (hsv-coef v) 1))
	      ((= 1 vgcd))
	      (t 
	       (dotimes (index n-nz)
		 (divf (aref vis index) vgcd))
	       (mulf (hsv-coef v) vgcd)))))))



(defun hsv-sort-indices-increasing (v)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (flet ((sift-down (root end)
	   (declare (fixnum root end))
	   (loop
	      (let* ((child (+ 1 (* 2 root)))
		     (child-i (+ child 1))
		     (root-i (+ root 1)))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (aref (hsv-is v) child-i) 
			      (aref (hsv-is v) (+ child-i 1))))
		  (incf child)
		  (incf child-i))
		(when (>= (aref (hsv-is v) root-i) 
			  (aref (hsv-is v) child-i))
		  (return))
		(rotatef (aref (hsv-is v) root-i) 
			 (aref (hsv-is v) child-i))
		(rotatef (aref (hsv-vis v) root-i) 
			 (aref (hsv-vis v) child-i))
		(setf root child)))))
    ;; max-heapify
    (let ((count (- (hsv-length v) 1)))
      (loop for heapify-start from (floor (- count 2) 2) downto 0
	 do (sift-down heapify-start (- count 1)))
      ;; sort
      (loop for end from (- count 1) downto 1
	 do (let ((end-i (+ end 1)))
	      (rotatef (aref (hsv-is v) 1) 
		       (aref (hsv-is v) end-i))
	      (rotatef (aref (hsv-vis v) 1) 
		       (aref (hsv-vis v) end-i))
		   (sift-down 0 (- end 1)))))))


;;;;
(defun hsv-remove-zeros (v)
  (declare (optimize (debug 0) (safety 0) (speed 1)))
  (let ((nz-ci 0)
	(v-is (hsv-is v))
	(v-vis (hsv-vis v))
	(n (hsv-length v)))
    (dotimes (ci n)
      (unless (zerop (aref v-vis ci))
	(unless (= nz-ci ci)
	  (setf (aref v-is nz-ci) (aref v-is ci))
	  (setf (aref v-vis nz-ci) (aref v-vis ci)))
	(incf nz-ci)))
    (unless (= n nz-ci)
      (setf (fill-pointer v-is) nz-ci
	    (fill-pointer v-vis) nz-ci))))
	  
