(defvar my-mps)
(defvar my-lp)
(defvar my-b)


(defun pip (lp b)
  (let ((header (make-nvector (basis-size b) 0 fixnum)))
    (loop for j from 0 below (length header) do (setf (aref header j) j))
    (loop
       (fill-basis b lp header)
       (cond ((preassigned-pivot b)
	      (return t))
	     ((eq 'redundant-column (basis-is-singular b))
	      (return (values (basis-is-singular b) (basis-singular-ref b))))
	     (t
	      (decf (basis-size b))
	      (lp-remove-row lp (aref (lp-active-row-refs lp) (basis-singular-ref b)))
	      (vector-pop header))))))

(setf my-mps (load-from-mps "/Users/mariusposta/Code/netlib/sc50a.mps"))

(defun run ()
  (setf my-lp (mps->lp my-mps))
  (preprocess my-lp)
  (setf my-b (make-basis :lp my-lp :ppivot-coef 0.0))
  (pip my-lp my-b)
  )
#|
  (print 'done-pip)
  (when (lu-decomposition my-b)
    (let ((ob (make-basis :lp my-lp :ppivot-coef 0.0))
	  (header (make-nvector (basis-size my-b) 0 fixnum)))
      (loop for j from 0 below (length header) do (setf (aref header j) j))
      (fill-basis ob my-lp header)
      (lu-check my-b ob))))
|#
