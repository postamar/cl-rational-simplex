(in-package :rationalsimplex)

;;;; Data structures for everything regarding hyper-sparse vectors


;;;; hyper-sparse vector data structures
(defstruct (hsv
	     (:constructor %make-hsv)
	     (:print-function print-hsv))
  (coef    1   :type rational)
  (length  0   :type fixnum)
  (prev-alloc-counter 8 :type fixnum)
  (current-alloc-counter 13 :type fixnum)
  (is      (error "hsv constructor") :type (simple-array fixnum 1))
  (vis     (error "hsv constructor") :type (simple-array integer 1)))

(defun print-hsv (v stream depth)
  (declare (ignore depth))
  (format stream "#HSV{ Length = ~D, Coef = ~A," (hsv-length v) (hsv-coef v))
  (dotimes (k (hsv-length v))
    (format stream " (#~D = ~D)" (aref (hsv-is v) k) (aref (hsv-vis v) k)))
  (format stream "}"))

(defstruct (hsv-float
	     (:constructor %make-hsv-float))
  (length  0   :type fixnum)
  (prev-alloc-counter 8 :type fixnum)
  (current-alloc-counter 13 :type fixnum)
  (is      (error "hsv float constructor") :type (simple-array fixnum 1))
  (vfs     (error "hsv float constructor") :type (simple-array double-float 1)))


;;;; HSV constructors
(defun make-hsv ()
  (%make-hsv
   :is  (make-array 13 :initial-element -1 :element-type 'fixnum)
   :vis (make-array 13 :initial-element 0 :element-type 'integer)))



(defun make-hsv-float ()
  (%make-hsv-float
   :is  (make-array 13 :initial-element -1 :element-type 'fixnum)
   :vfs (make-array 13 :initial-element 0.0d0 :element-type 'double-float)))
		 

 
;;;; Empties HSV 

(declaim (inline reset-hsv))
(defun reset-hsv (v)
  (setf (hsv-coef v) 1
	(hsv-length v) 0))

(declaim (inline reset-hsv-float))
(defun reset-hsv-float (vf)
  (setf (hsv-float-length vf) 0))
	

(defun copy-hsv-into-hsv (hsvsrc hsvdest)
  (let* ((n (hsv-length hsvsrc))
	 (c (hsv-coef hsvsrc)))
    (when (< (hsv-current-alloc-counter hsvdest) n)
      (let ((b (hsv-current-alloc-counter hsvdest))
	    (a (hsv-prev-alloc-counter hsvdest)))
	(loop
	   (when (<= n b)
	     (return))
	   (rotatef a b)
	   (incf b a))
	(setf (hsv-is hsvdest) (make-array b :initial-element -1 :element-type 'fixnum)
	      (hsv-vis hsvdest) (make-array b :initial-element 0 :element-type 'integer)
	      (hsv-prev-alloc-counter hsvdest) a
	      (hsv-current-alloc-counter hsvdest) b)))
    (setf (hsv-coef hsvdest) c
	  (hsv-length hsvdest) n)
    (replace (hsv-is hsvdest) (hsv-is hsvsrc) :end2 n)
    (replace (hsv-vis hsvdest) (hsv-vis hsvsrc) :end2 n)
    t))


(defun copy-hsv-into-hsv-float (hsv hsvf)
  (let* ((n (hsv-length hsv))
	 (c (hsv-coef hsv)))
    (when (< (hsv-float-current-alloc-counter hsvf) n)
      (let ((b (hsv-float-current-alloc-counter hsvf))
	    (a (hsv-float-prev-alloc-counter hsvf)))
	(loop
	   (when (<= n b)
	     (return))
	   (rotatef a b)
	   (incf b a))
	(setf (hsv-float-is hsvf) (make-array b :initial-element -1 :element-type 'fixnum)
	      (hsv-float-vfs hsvf) (make-array b :initial-element 0.0d0 :element-type 'double-float)
	      (hsv-float-prev-alloc-counter hsvf) a
	      (hsv-float-current-alloc-counter hsvf) b)))
    (reset-hsv-float hsvf)
    (setf (hsv-float-length hsvf) n)
    (replace (hsv-float-is hsvf) (hsv-is hsv) :end2 n)
    (dotimes (ci n t)
      (setf (aref (hsv-float-vfs hsvf) ci)
	    (coerce (* c (aref (hsv-vis hsv) ci)) 'double-float)))))
    

(declaim (inline hsv-ratio))
(defun hsv-ratio (v ci)
  (* (hsv-coef v) (aref (hsv-vis v) ci)))


(defun hsv-normalize (v)
  (unless (zerop (hsv-length v))
    (let ((vgcd 0))
      (dotimes (index (hsv-length v))
	(let ((vi (aref (hsv-vis v) index)))
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
	     (dotimes (index (hsv-length v))
	       (divf (aref (hsv-vis v) index) vgcd))
	     (mulf (hsv-coef v) vgcd))))))



(defun hsv-sort-indices-increasing (v)
  (flet ((sift-down (root end)
	   (loop
	      (let ((child (+ 1 (* 2 root))))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (aref (hsv-is v) child) 
			      (aref (hsv-is v) (+ child 1))))
		  (incf child))
		(when (>= (aref (hsv-is v) root) 
			  (aref (hsv-is v) child))
		  (return))
		(rotatef (aref (hsv-is v) root) 
			 (aref (hsv-is v) child))
		(rotatef (aref (hsv-vis v) root) 
			 (aref (hsv-vis v) child))
		(setf root child)))))
    ;; max-heapify
    (let ((count (hsv-length v)))
      (loop for heapify-start from (floor (- count 2) 2) downto 0
	 do (sift-down heapify-start (- count 1)))
      ;; sort
      (loop for end from (- count 1) downto 1
	 do (progn 
	      (rotatef (aref (hsv-is v) 0) 
		       (aref (hsv-is v) end))
	      (rotatef (aref (hsv-vis v) 0) 
		       (aref (hsv-vis v) end))
	      (sift-down 0 (- end 1)))))))


;;;;
(defun hsv-remove-zeros (v)
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
      (setf (hsv-length v) nz-ci))))




;;;;
(defun hsv-add (ind val hsv)
  (when (<= (hsv-current-alloc-counter hsv) (hsv-length hsv))
    (rotatef (hsv-current-alloc-counter hsv) (hsv-prev-alloc-counter hsv))
    (incf (hsv-current-alloc-counter hsv) (hsv-prev-alloc-counter hsv))
    (let ((vi (make-array (hsv-current-alloc-counter hsv) :initial-element -1 :element-type 'fixnum))
	  (vv (make-array (hsv-current-alloc-counter hsv) :initial-element 0 :element-type 'integer)))
      (replace vi (hsv-is hsv))
      (replace vv (hsv-vis hsv))
      (setf (hsv-is hsv) vi
	    (hsv-vis hsv) vv)))
  (setf (aref (hsv-is hsv) (hsv-length hsv)) ind
	(aref (hsv-vis hsv) (hsv-length hsv)) val)
  (incf (hsv-length hsv)))


;;;;
(defun hsv-remove (ci hsv)
  (let ((last-ci (- (hsv-length hsv) 1)))
    (when (< ci last-ci)
      (rotatef (aref (hsv-is hsv) ci) (aref (hsv-is hsv) last-ci))
      (rotatef (aref (hsv-vis hsv) ci) (aref (hsv-vis hsv) last-ci)))
    (decf (hsv-length hsv))))
      

(defun hsv-find (ind hsv)
  (find-index-bounded (hsv-is hsv) (hsv-length hsv) ind))


