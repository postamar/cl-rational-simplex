#|
(defun check-yB=c_B (y)
  (dotimes (header-index (lp-data-m *lp*))
    (let ((y.a_j 0)
	  (j (aref *basis-header* header-index)))
      (let ((a_j-indices (aref (lp-data-a-indices *lp*) j))
	    (a_j-values (aref (lp-data-a-values *lp*) j)))
	(dotimes (k (length a_j-indices))
	  (incf y.a_j
		(* (aref a_j-values k)
		   (aref y (aref a_j-indices k))))))
      (assert (= y.a_j
		 (* (lp-data-obj-sense *lp*)
		   (aref (lp-data-c *lp*) j)))))))

|#

#|
(defun check-Bd=a_j (d j)
  (let ((a_j-indices (aref (lp-data-a-indices *lp*) j))
	(a_j-values (aref (lp-data-a-values *lp*) j))
	(a_j (make-array (lp-data-m *lp*) 
			 :initial-element 0 
			 :element-type 'rational))
	(B (make-array (list (lp-data-m *lp*) (lp-data-m *lp*)) 
		       :initial-element 0 
		       :element-type 'rational)))
    (dotimes (index (length a_j-indices)) 
      (setf (aref a_j (aref a_j-indices index))
	    (aref a_j-values index)))
    (dotimes (header-index (length *basis-header*))
      (let ((i (aref *basis-header* header-index)))
	(let ((a_i-indices (aref (lp-data-a-indices *lp*) i))
	      (a_i-values (aref (lp-data-a-values *lp*) i)))
	  (dotimes (index (length a_i-indices)) 
	    (setf (aref B (aref a_i-indices index) header-index)
		  (aref a_i-values index))))))
    (dotimes (i (length *basis-header*))
      (let ((b_i.d 0))
	(dotimes (k (length *basis-header*))
	  (incf b_i.d
		(* (aref d k)
		   (aref B i k))))
	(assert (= b_i.d (aref a_j i)))))))
|#

#|
(defun check-LU ()
  (let* ((m (length *basis-header*))
	 (B (make-array (list m m) :initial-element 0 :element-type 'rational)))
    ;; fill B
    (dotimes (header-index m)
      (let* ((j (aref *basis-header* header-index))
	     (a_j-indices (aref (lp-data-a-indices *lp*) j))
	     (a_j-values (aref (lp-data-a-values *lp*) j)))
	(dotimes (index (length a_j-indices))
	  (setf (aref B (aref a_j-indices index) header-index)
		(aref a_j-values index)))))
    ;; perform left-multiplications
    (dotimes (k m)
      (print B)
      ;; left-multiply by P_k, i.e. exchange rows
      (let ((index (aref *eta-basis-permutation-row-indices* k))
	    (temp 0))
	(dotimes (j m)
	  (setf temp (aref B k j))
	  (setf (aref B k j) (aref B index j))
	  (setf (aref B index j) temp)))
      ;; left-multiply by L_k
      (let ((e_k (aref *eta-basis-lower-columns* k)))
	(dotimes (j m)
	  (let ((B_kj 0))
	    (dotimes (i m)
	      (incf B_kj
		    (* (aref e_k i)
		       (aref B i j))))
	    (setf (aref B k j) B_kj)))))
    ;; compare with U
    (print B)
    (print *eta-basis-upper-columns*)
    (dotimes (j m)
      (let ((u_j (aref *eta-basis-upper-columns* j)))
	(dotimes (i m)
	  (assert (= (aref u_j i)
		     (aref B i j))))))))
|#

;;;; simple test instances

(defun make-very-simple-lp-instance-1 ()
  (make-lp-data
   :lp-name "test1"
   :obj-name "obj"
   :obj-sense 1
   :row-names #("row1" "row2")
   :col-names #("x1" "x2" "slack1" "slack2")
   :n 4
   :m 2
   :c #(1000 1200 0 0)
   :b #(160 120)
   :A-indices #(#(0 1) #(0 1) #(0) #(1))
   :A-values #(#(8 4) #(4 6) #(1) #(1))
   :l #(0 0 0 0)
   :u #(34 14 0 0)
   :l-p #*1111
   :u-p #*1100))

	     
(defun make-very-simple-lp-instance-2 ()
  (make-lp-data
   :lp-name "test2"
   :obj-name "obj"
   :obj-sense -1
   :row-names #("row1" "row2")
   :col-names #("x1" "x2" "x3" "x4" "x5")
   :n 5
   :m 2
   :c #(2 1 3 -2 10)
   :b #(5 9)
   :A-indices #(#(0) #(1) #(0 1) #(0 1) #(0 1))
   :A-values #(#(1) #(1) #(1 2) #(-1 2) #(2 1))
   :l #(0 0 0 0 0)
   :u #(7 10 1 5 3)
   :l-p #*11111
   :u-p #*11111))


(defun make-very-simple-lp-instance-3 ()
  (make-lp-data
   :lp-name "test3"
   :obj-name "obj"
   :obj-sense -1
   :row-names #("row1" "row2")
   :col-names #("x1" "x2" "x3" "x4" "slack1" "slack2")
   :n 6
   :m 2
   :c #(-4 -3 -1 -2 0 0)
   :b #(5 4)
   :A-indices #(#(0 1) #(0 1) #(0 1) #(0 1) #(0) #(1))
   :A-values #(#(4 3) #(2 1) #(1 2) #(1 1) #(1) #(1))
   :l #(0 0 0 0 0 0)
   :u #(0 0 0 0 0 0)
   :l-p #*111111
   :u-p #*000000))


 (defun make-very-simple-lp-instance-4 ()
  (make-lp-data
   :lp-name "test4"
   :obj-name "obj"
   :obj-sense -1
   :row-names #("row1" "row2" "row3")
   :col-names #("x1" "x2" "x3")
   :n 3
   :m 3
   :c #(1 0 0)
   :b #(2 1 7)
   :A-indices #(#(0 1 2) #(0 1 2) #(0 1 2))
   :A-values #(#(1 -1 2) #(2 3 -3) #(1 1 4))
   :l #(0 0 0)
   :u #(0 0 0)
   :l-p #*111
   :u-p #*000))   


(defun make-very-simple-lp-instance-5 ()
  (make-lp-data
   :lp-name "test5"
   :obj-name "obj"
   :obj-sense -1
   :row-names #("row1" "row2" "row3" "row4")
   :col-names #("x1" "x2" "x3" "x4" "x5" "x6")
   :n 6
   :m 4
   :c #(-8 -9 -7 -6 -8 -9)
   :b #(4 2 2 1)
   :A-indices #(#(0 2) #(1 3) #(0) #(1 2) #(0 3) #(1))
   :A-values #(#(1 1) #(1 1) #(1) #(1 1) #(1 1) #(1))
   :l #(0 0 0 0 0 0)
   :u #(0 0 0 0 0 0)
   :l-p #*111111
   :u-p #*000000))
