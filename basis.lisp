;;;;
;;;;



;;;; Data structure for storing an LP basis
(defstruct (basis
	     (:constructor %make-basis))
  (size        0   :type fixnum)
  (header      #() :type vector)
  (col-row-ref #() :type vector)
  (col-val-ref #() :type vector)
  (l-file      #() :type vector)
  (p-file      #() :type vector)
  (u-file      #() :type vector)
  (eta-file    #() :type vector))

  

;;;; Constructor
(defun make-empty-basis (m)
  (let* ((eta (%make-eta-matrix))
	 (col-row-ref (make-array m :initial-element #() :element-type 'vector))
	 (col-val-ref (make-array m :initial-element #() :element-type 'vector))
	 (l-file (make-array m :initial-element eta :element-type 'eta-matrix))
	 (p-file (make-array m :initial-element 0   :element-type 'fixnum))
	 (u-file (make-array m :initial-element eta :element-type 'eta-matrix)))
    (dotimes (k m)
      (setf (aref col-row-ref k) (vector k)
	    (aref col-val-ref k) (vector 1)
	    (aref p-file k) k
	    (aref l-file k) (make-eta-matrix)
	    (aref u-file k) (make-eta-matrix)))
    (dotimes (k m)
      (reset-eta-matrix (aref l-file k) k)
      (vector-push-extend k (eta-matrix-row-ref (aref l-file k)))
      (vector-push-extend 1 (eta-matrix-val-num (aref l-file k)))
      (vector-push-extend 1.0 (eta-matrix-val-flt (aref l-file k)))
      (reset-eta-matrix (aref u-file k) k)
      (vector-push-extend k (eta-matrix-row-ref (aref u-file k)))
      (vector-push-extend 1 (eta-matrix-val-num (aref u-file k)))
      (vector-push-extend 1.0 (eta-matrix-val-flt (aref u-file k))))
    (%make-basis
     :size m
     :header      (make-array m :initial-element -1  :element-type 'fixnum)
     :col-row-ref col-row-ref 
     :col-val-ref col-val-ref 
     :l-file      l-file
     :p-file      p-file
     :u-file      u-file
     :eta-file    (make-vector eta-matrix))))
    
