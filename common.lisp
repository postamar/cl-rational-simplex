;;;; Common functions and macros
;;;;



(defun absmax (a b)
  (max (abs a) (abs b)))

(define-modify-macro absmaxf (value) absmax)
(define-modify-macro maxf (value) max)
(define-modify-macro minf (value) min)
(define-modify-macro mulf (value) *)
(define-modify-macro divf (value) /)
(define-modify-macro ashf (value) ash)

(defmacro orf (place &rest rest)
  `(setf ,place (or ,@rest ,place)))

(defmacro make-vector (&optional (type t))
  `(make-array 0 :adjustable t :fill-pointer t :element-type (quote ,type)))

(defmacro make-nvector (n ielt &optional (type t))
  `(make-array ,n :initial-element ,ielt 
	       :adjustable t :fill-pointer t :element-type (quote ,type)))

(defun linear-modify (vector old-value new-value)
  (let ((len (length vector)))
    (dotimes (k len)
      (when (= old-value (aref vector k))
	(setf (aref vector k) new-value)
	(return t)))))

(defun find-index (vector value)
  (declare ((array fixnum *) vector)
	   (fixnum value))
  (let ((l 0)
	(u (length vector))
	(p 0)
	(d 0))
    (declare (fixnum l u p d))
    (loop
       (setf p (ash (- u l) -1))
       (when (zerop p)
	 (if (eql value (aref vector l))
	     (return l)
	     (return -1)))
       (incf p l)
       (setf d (aref vector p))
       (if (eql value d)
	   (return p)
	   (if (< d value)
	       (setf l p)
	       (setf u p))))))

(defun insert-in-increasing-vector (value vector)
  (declare ((array fixnum *) vector)
	   (fixnum value))
  (let ((len (length vector)))
    (if (zerop len)
	(vector-push-extend value vector)
	(let* ((lasti (- len 1))
	       (lastv (aref vector lasti)))
	  (if (< value lastv)
	      (let ((insi -1))
		(dotimes (i len)
		  (let ((v (aref vector i)))
		    (cond 
		      ((= v value)
		       (return-from insert-in-increasing-vector i))
		      ((< value v)
		       (setf insi i)
		       (return)))))
		(vector-push-extend lastv vector)
		(loop for i from lasti above insi
		   do (setf (aref vector i) (aref vector (- i 1))))
		(setf (aref vector insi) value)
		insi)
	      (if (= value lastv)
		  lasti
		  (vector-push-extend value vector)))))))
    
    
	   


(defun remove-in-increasing-vector (value vector)
  (declare ((array fixnum *) vector)
	   (fixnum value))
  (let ((len (length vector)))
    (unless (zerop len)
      (let ((lasti (- len 1))
	    (remi (find-index vector value)))
	(unless (= -1 remi)
	  (loop for i from remi below lasti
	     do (setf (aref vector i) (aref vector (+ i 1))))
	  (decf (fill-pointer vector))
	  t)))))

  
       
(defun floor-log2 (n)
  (let ((p 0))
    (unless (< n 65536)
      (ashf n -16)
      (incf p 16))
    (unless (< n 256)
      (ashf n -8)
      (incf p 8))
    (unless (< n 16)
      (ashf n -4)
      (incf p 4))
    (unless (< n 4)
      (ashf n -2)
      (incf p 2))
    (unless (< n 2)
      (ashf n -1)
      (incf p))
    (if (zerop n)
	-1
	p)))



;;;;
(defun in-place-sort-keys-increasing (keys values)
  (flet ((sift-down (root end)
	   (loop
	      (let ((child (+ 1 (* 2 root))))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (aref keys child) (aref keys (+ child 1))))
		  (incf child))
		(when (>= (aref keys root) (aref keys child))
		  (return))
		(rotatef (aref keys root) (aref keys child))
		(rotatef (aref values root) (aref values child))
		(setf root child)))))
    ;; max-heapify
    (let ((count (length keys)))
      (loop for heapify-start from (floor (- count 2) 2) downto 0
	 do (sift-down heapify-start (- count 1)))
      ;; sort
      (loop for end from (- count 1) downto 1
	 do (progn (rotatef (aref keys 0) (aref keys end))
		   (rotatef (aref values 0) (aref values end))
		   (sift-down 0 (- end 1)))))))
      
		     
	    
	     
		     
