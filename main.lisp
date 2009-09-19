(defvar my-mps)
(defvar my-lp)
(defvar my-b)
(defvar my-sd)

#|
(setf 
 my-mps 
 #S(MPS 
    :lp-name "TEST1"
    :obj-spec ("COST" . -1)
    :row-spec (("ROW1" . <=) ("ROW2" . <=) ("ROW3" . <=))
    :row-rhss (("ROW1" . 7) ("ROW2" . 12) ("ROW3" . 10))
    :col-spec (("X1" ("COST" . 1) ("ROW1" . 3) ("ROW2" . -2) ("ROW3" . -4))
	       ("X2" ("COST" . -3) ("ROW1" . -1) ("ROW2" . 4) ("ROW3" . 3))
	       ("X3" ("COST" . 2) ("ROW1" . 2) ("ROW3". 8)))))
|#




(defun load-mps (name)
  (when 
      (setf my-mps 
	    (load-from-mps 
	     (concatenate 'string
			  "/Users/mariusposta/Code/netlib/"
			  name
			  ".mps")))
    t))

(defun run ()
  (setf *random-state* (make-random-state *simplex-random-state*))
  (setf my-lp (mps->lp my-mps))
  (setf my-b (make-phase1-initial-basis my-lp))
  (setf my-sd (make-simplex my-lp my-b)))
  ;(dual-simplex my-sd))


