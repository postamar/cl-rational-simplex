;;;; Common functions and macros
;;;;



(defun absmax (a b)
  (max (abs a) (abs b)))

(define-modify-macro absmaxf (value) absmax)
(define-modify-macro maxf (value) max)
(define-modify-macro minf (value) min)
(define-modify-macro mulf (value) *)
(define-modify-macro divf (value) /)

(defmacro orf (place &rest rest)
  `(setf ,place (or ,@rest ,place)))

(defmacro make-vector (&optional (type t))
  `(make-array 0 :adjustable t :fill-pointer t :element-type (quote ,type)))

(defmacro make-nvector (n ielt &optional (type t))
  `(make-array ,n :initial-element ,ielt 
	       :adjustable t :fill-pointer t :element-type (quote ,type)))

(defun find-index (vector value)
  (declare ((simple-array fixnum *) vector)
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


(defun remove-from-adjustable-vector (adjvector value)
  (let ((l 0)
	(r (length adjvector)))
    (unless (zerop r)
      (decf r)
      (loop
	 (cond ((<= r l)
		(if (= value (aref adjvector l))
		    (setf (fill-pointer adjvector) l)
		    (setf (fill-pointer adjvector) (+ l 1)))
		(return))
	       ((= value (aref adjvector r))
		(decf r))
	       ((/= value (aref adjvector l))
		(incf l))
	       (t
		(let ((temp (aref adjvector l)))
		  (setf (aref adjvector l) (aref adjvector r))
		  (setf (aref adjvector r) temp))
		(decf r)
		(incf l)))))))
