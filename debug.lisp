(defvar my-mps)
(defvar my-lp)
(defvar my-b)
(defvar my-ob)

(setf my-mps (load-from-mps "/Users/mariusposta/Code/netlib/afiro.mps"))

(defun test-lu (lp b ob)
  (let ((header (make-nvector (basis-size b) 0 fixnum)))
    (loop for j from 0 below (length header) do (setf (aref header j) j))
    (loop
       (format t "basis size = ~A~%" (basis-size b))
       (unless (and (fill-basis b lp header)
		    (basis-lu-decomposition b))
	 (format t "~A ref ~A~%"
		 (basis-is-singular b) (aref header (basis-singular-ref b))))
       (cond ((eq 'redundant-column (basis-is-singular b))
	      (return (values (basis-is-singular b) (basis-singular-ref b))))
	     ((eq 'redundant-row (basis-is-singular b))
	      (decf (basis-size b))
	      (decf (basis-size ob))
	      (lp-remove-row lp (aref (lp-active-row-refs lp) (basis-singular-ref b)))
	      (vector-pop header))
	     (t 
	      (fill-basis ob lp header)
	      (return (lu-check b ob)))))))


(defun run ()
  (setf my-lp (mps->lp my-mps))
  (preprocess my-lp)
  (setf my-b (make-basis :lp my-lp :ppivot-coef 0.0))
  (setf my-ob (make-basis :lp my-lp :ppivot-coef (basis-ppivot-coef my-b)))
  (test-lu my-lp my-b my-ob))
