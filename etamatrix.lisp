;;;; Eta matrix implementation



;;;; Holds an eta matrix as a sparse column with index 
;;;; values are in common denominator rationals and in floats
(defstruct (eta-matrix
	     (:constructor %make-eta-matrix))
  (col-ref   -1  :type fixnum)
  (row-ref   #() :type vector)
  (val-num   #() :type vector)
  (val-flt   #() :type vector)
  (col-denom 1   :type integer))



;;;; Constructor
(defun make-eta-matrix ()
  (%make-eta-matrix
   :row-ref (make-vector fixnum)
   :val-num (make-vector integer)
   :val-flt (make-vector float)))



;;;;
(defun reset-eta-matrix (eta &optional (column-index -1))
  (setf (eta-matrix-col-ref eta) column-index
	(fill-pointer (eta-matrix-row-ref eta)) 0
	(fill-pointer (eta-matrix-val-num eta)) 0
	(fill-pointer (eta-matrix-val-flt eta)) 0
	(eta-matrix-col-denom eta) 1))
