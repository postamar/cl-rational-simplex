(defvar my-mps)
(defvar my-lp)
(defvar my-b)
(defvar my-obm)
(defvar my-sd)


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


(setf my-mps (load-from-mps "/Users/mariusposta/Code/netlib/afiro.mps"))

(defun run ()
  (setf my-lp (mps->lp my-mps))
  (preprocess my-lp)
  (setf my-b (make-phase1-initial-basis my-lp))
  (setf my-obm (make-basis-matrix :lp my-lp))
  (fill-basis-matrix my-obm my-lp (basis-header my-b))
  (setf my-sd (make-simplex my-lp my-b))
  t)



