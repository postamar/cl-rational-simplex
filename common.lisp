;;;;; Common data structures for the rational simplex
;;;;;


;;; This structure holds all the data describing a LP in standard form.
;;; An lp-data object is meant to be read-only.
;;; The matrix A is represented as an array of column vectors.
;;; The column vectors are represented as sparse vectors
;;; Columns corresponding to slack variables have the name "".
;;; Bounds l and u are rationals or '-inf or '+inf.

(defstruct lp-data
   lp-name
   obj-name
   obj-sense
   row-names
   col-names
   n
   m
   c
   b
   A-indices
   A-values
   l
   u
   l-p
   u-p)


;;;; Alternative constructor for lp-data
(defun make-new-lp-data (n m)
  (make-lp-data
   :lp-name   ""
   :obj-name  ""
   :obj-sense 0 
   :row-names (make-array m :initial-element "" :element-type 'string)
   :col-names (make-array n :initial-element "" :element-type 'string)
   :n         n
   :m         m
   :c         (make-array n :element-type 'rational)
   :b         (make-array m :element-type 'rational)
   :a-indices (make-array n :initial-element nil)
   :a-values  (make-array n :initial-element nil)
   :l         (make-array n :initial-element 0 :element-type 'rational) 
   :u         (make-array n :initial-element 0 :element-type 'rational)
   :l-p       (make-array n :initial-element 1 :element-type 'bit)
   :u-p       (make-array n :initial-element 0 :element-type 'bit)))


;;;; Makes a new lp-data object from lp with specified rows removed
(defun copy-lp-remove-rows (lp row-list)
  (let ((n (lp-data-n lp))
	(m (- (lp-data-m lp) (length row-list)))
	(row-table (make-array (lp-data-m lp) :initial-element 0 :element-type 'bit))
	(col-table (make-array (lp-data-n lp) :initial-element 0 :element-type 'bit)))
    (dolist (i row-list)
      (setf (bit row-table i) 1))
    (dotimes (j (lp-data-n lp))
      (let ((a_j-indices (aref (lp-data-a-indices lp) j)))
	(when (dotimes (k (length a_j-indices) t)
		(when (zerop (bit row-table (aref a_j-indices k)))
		  (return nil)))
	  (decf n)
	  (setf (bit col-table j) 1))))
    (let ((new-lp (make-new-lp-data n m)))
      (setf (lp-data-lp-name new-lp) 
	    (concatenate 'string (lp-data-lp-name lp) "-PHASE2"))
      (setf (lp-data-obj-name new-lp)
	    (lp-data-obj-name lp))
      (setf (lp-data-obj-sense new-lp)
	    (lp-data-obj-sense lp))
      (let ((new-i 0))
	(dotimes (i (lp-data-m lp))
	  (when (zerop (bit row-table i))
	    (dolist (fun (list #'lp-data-row-names
			       #'lp-data-b))
	      (setf (aref (funcall fun new-lp) new-i)
		    (aref (funcall fun lp) i)))
	    (incf new-i))))
      (let ((new-j 0))
	(dotimes (j (lp-data-n lp))
	  (when (zerop (bit col-table j))
	    (dolist (fun (list #'lp-data-col-names
			       #'lp-data-c
			       #'lp-data-a-indices
			       #'lp-data-a-values
			       #'lp-data-l
			       #'lp-data-u
			       #'lp-data-l-p
			       #'lp-data-u-p))
	      (setf (aref (funcall fun new-lp) new-j)
		    (aref (funcall fun lp) j)))
	    (incf new-j))))
      new-lp)))
			    


;;;; outputs model from lp-data
(defun output-model (lp)

  (format t "linear program name: ~A~%" (lp-data-lp-name lp))
    
  (format t (if (= -1 (lp-data-obj-sense lp)) " Min  " " Max  "))
  (let ((print+ nil))
    (dotimes (j (lp-data-n lp) (format t "~%subject to~%"))
      (let ((c_j (aref (lp-data-c lp) j)))
	(unless (zerop c_j)
	  (if print+ 
	      (format t " + ")
	      (setf print+ t))
	  (format t "(~A ~A)" c_j (aref (lp-data-col-names lp) j))))))

  (dotimes (i (lp-data-m lp))
    (let ((print+ nil))
      (format t "~A:  " (aref (lp-data-row-names lp) i))
      (dotimes (j (lp-data-n lp))
	(let ((a_ij (dotimes (k (length (aref (lp-data-A-indices lp) j)) 0)
		      (when (= (aref (aref (lp-data-A-indices lp) j) k) i)
			(return (aref (aref (lp-data-A-values lp) j) k))))))
	  (unless (zerop a_ij)
	    (if print+
		(format t " + ")
		(setf print+ t))
	    (format t "(~A ~A)" a_ij (aref (lp-data-col-names lp) j))))))
    (format t " = ~A~%" (aref (lp-data-b lp) i)))

  (dotimes (j (lp-data-n lp))
    (let ((name (aref (lp-data-col-names lp) j))
	  (l (aref (lp-data-l lp) j))
	  (u (aref (lp-data-u lp) j))
	  (l-p (bit (lp-data-l-p lp) j))
	  (u-p (bit (lp-data-u-p lp) j)))
      (cond ((and (zerop l-p)
		  (zerop u-p)))
	    ((zerop l-p)
	     (format t "~A <= ~A" name u))
	    ((zerop u-p)
	     (format t "~A <= ~A" l name))
	    (t
	     (format t "~A <= ~A <= ~A" l name u)))
      (format t "~%"))))


;;;; Copies elements in src into dst as much as possible
(defun copy-into-array (src dst)
  (let ((i-max (min (length src)
		    (length dst))))
    (dotimes (i i-max)
      (setf (aref dst i) (aref src i)))))

