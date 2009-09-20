(in-package :rationalsimplex)

;;;;; Multithreading support definitions
;;;;;
;;;;; Macro definitions which allow the offloading of 
;;;;; expensive computations to seperate threads, implementation-dependent
;;;;; 

;(defconstant +thread-system+ (car (intersection *features* '(:sb-thread))))
(defconstant +thread-system+ nil)

;;;;; Thread operations

(defmacro thread-launch (thread-holder fn)
  (let* ((str (princ-to-string thread-holder))
	 (name (subseq str 1 (1- (length str)))))
    (case +thread-system+
      ((nil)
       `(setf ,thread-holder (multiple-value-list (funcall ,fn))))
      ((:sb-thread) 
       `(setf ,thread-holder (sb-thread:make-thread ,fn :name ,name))))))
  

(defmacro thread-result (thread-holder)
  `(prog1
       ,(case +thread-system+
	      ((nil)
	       `(values-list ,thread-holder))
	      ((:sb-thread) 
	       `(sb-thread:join-thread ,thread-holder)))
     (setf ,thread-holder nil)))


;;;;; Rational simplex thread lock definitions
  
(defparameter *dse-ftran-thread* nil)
(defparameter *basis-matrix-factor-thread* nil)
(defparameter *dse-weight-update-thread* nil)
