;;;; MPS file parser

;;; data structure containing everything in the MPS file
;;; in a non-awkward manner
(defstruct mps 
  lp-name
  obj-spec
  row-spec
  row-rhss
  col-lbounds
  col-ubounds
  col-spec)

   

;;; scans one line from the MPS file into a list
(defun listify (line)
  (let ((prec 0)
	(space nil)
	(list '())
	(line-length (length line)))
    (dotimes (i line-length 
	      (nreverse (if (= prec i) 
			    list 
			    (cons (subseq line prec line-length) list))))
      (if (or (char-equal (char line i) #\Space)
	      (char-equal (char line i) #\Tab))
	  (if space
	      (incf prec)
	      (setf list (if (= prec i) 
			     list 
			     (cons (subseq line prec i) list))
		    space t 
		    prec (+ i 1)))
	  (setf space nil)))))
	      


;;; reads the next non-empty line in the MPS file
(defun get-next-line-list (mps-stream)
  (multiple-value-bind (line eof) (read-line mps-stream nil t)
	  (if eof 
	      nil
	      (let ((line-list (listify line)))
		(if line-list
		    line-list
		    (get-next-line-list mps-stream))))))
      


;;; reads a block in the MPS file (rows, columns, etc.)
(defmacro read-mps-block (mps-stream block end-tag)
  (let ((line-list (gensym)))
    `(loop
	(let ((,line-list (get-next-line-list ,mps-stream)))
	  (unless ,line-list
	    (setf ,block nil)
	    (return))
	  (when (or (string= ,end-tag (car ,line-list))
		    (string= "ENDATA" (car ,line-list)))
	    (return))
	  (push ,line-list ,block)))))



;;; generates an association list for the data in the column block 
(defun make-coef-alist (list)
  (let ((length (length list)))
    (cond ((= length 2)
	   (destructuring-bind (n1 v1) list
	     (list (cons n1 (rationalize (read-from-string v1))))))
	  ((= length 4)
	   (destructuring-bind (n1 v1 n2 v2) list
	     (list (cons n2 (rationalize (read-from-string v2)))
		   (cons n1 (rationalize (read-from-string v1))))))
	  (t
	   (error "parse error in MPS file, column section")))))


;;; turns the data in the mps structure into something actually usable, in standard form
;;; TODO: model checking
;;; - type the bound arrays properly
;;; - very basic preprocessing
(defun mps->raw (mps-data n-named n-slack m)
  (let ((raw (make-lp-data))
	(n (+ n-named n-slack))
	(slack-col n-named)
	(row-factor (make-array m :initial-element 1 :element-type 'rational))
	(row-index-by-name (make-hash-table :test 'equal))
	(col-index-by-name (make-hash-table :test 'equal)))

    (setf (lp-data-lp-name raw) (mps-lp-name mps-data)
	  (lp-data-obj-name raw) (car (mps-obj-spec mps-data))
	  (lp-data-obj-sense raw) (cdr (mps-obj-spec mps-data))
	  (lp-data-n raw) n
	  (lp-data-m raw) m
	  (lp-data-c raw) (make-array n :initial-element 0 :element-type 'rational)
	  (lp-data-b raw) (make-array m :initial-element 0 :element-type 'rational)
	  (lp-data-l raw) (make-array n :initial-element 0 :element-type 'rational)
	  (lp-data-u raw) (make-array n :initial-element 0 :element-type 'rational)
	  (lp-data-l-p raw) (make-array n :initial-element 1 :element-type 'bit)
	  (lp-data-u-p raw) (make-array n :initial-element 0 :element-type 'bit)
	  (lp-data-row-names raw) (make-array m :initial-element "" :element-type 'string)
	  (lp-data-col-names raw) (make-array n :initial-element "" :element-type 'string)
	  (lp-data-A-indices raw) (make-array n)
	  (lp-data-A-values raw) (make-array n))


    ;; fill in row name array
    (let ((i 0))
      (dolist (row (mps-row-spec mps-data))
	(setf (aref (lp-data-row-names raw) i) (car row))
	(setf (gethash (car row) row-index-by-name) i)
	(incf i)))

    ;; fill in column name array
    (let ((i 0))
      (dolist (col (mps-col-spec mps-data))
	(setf (aref (lp-data-col-names raw) i) (car col))
	(setf (gethash (car col) col-index-by-name) i)
	(incf i)))
    (dotimes (j n-slack)
      (setf (aref (lp-data-col-names raw) (+ n-named j))
	    (concatenate 'string "slack" (princ-to-string j))))

    ;; fill b and check for non-negativity
    (dolist (rhs (mps-row-rhss mps-data))
      (let ((i (gethash (car rhs) row-index-by-name))
	    (b_i (cdr rhs)))
	(when (< b_i 0)
	  (setf b_i (- b_i)
		(aref row-factor i) -1))
	(setf (aref (lp-data-b raw) i) b_i)))

    ;; add slack variable coefficients to A
    (dolist (row (mps-row-spec mps-data))
      (unless (eq (cdr row) '=)
	(let ((i (gethash (car row) row-index-by-name)))
	  (when (eq (cdr row) '>=)
	    (setf (aref row-factor i) (- (aref row-factor i))))
	  (setf (aref (lp-data-A-indices raw) slack-col)
		(make-array 1 
			    :initial-element i
			    :element-type 'fixnum))
	  (setf (aref (lp-data-A-values raw) slack-col)
		(make-array 1
			    :initial-element (aref row-factor i)
			    :element-type 'rational))
	  (incf slack-col))))

    ;; fill c and A
    (dolist (col (mps-col-spec mps-data))
      (let ((j (gethash (car col) col-index-by-name))
	    (c 0)
	    (alist '()))
	(dolist (coef (cdr col))
	  (destructuring-bind (row-name . a_ij) coef
	    (if (equal row-name (lp-data-obj-name raw))
		(setf (aref (lp-data-c raw) j) a_ij)
		(push (cons (gethash row-name row-index-by-name) a_ij) alist))))
	(let ((alist-length (length alist)))
	  (setf (aref (lp-data-A-indices raw) j)
		(make-array alist-length :initial-element -1 :element-type 'fixnum)
		(aref (lp-data-A-values raw) j)
		(make-array alist-length :initial-element 0 :element-type 'rational)))
	(setf alist (sort alist #'< :key #'car))
	(dolist (pair alist) 
	  (destructuring-bind (i . a_ij) pair
	    (setf (aref (aref (lp-data-A-indices raw) j) c) i
		  (aref (aref (lp-data-A-values raw) j) c) a_ij)
	    (incf c)))))
	    

    ;; fill l and u
    (dolist (lbound (mps-col-lbounds mps-data))
      (let ((i (gethash (car lbound) col-index-by-name))
	    (l_i (cdr lbound)))
	(setf (bit (lp-data-l-p raw) i) 1
	      (aref (lp-data-l raw) i) l_i)))
	
    (dolist (ubound (mps-col-ubounds mps-data))
      (let ((i (gethash (car ubound) col-index-by-name))
	    (u_i (cdr ubound)))
	(setf (bit (lp-data-u-p raw) i) 1
	      (aref (lp-data-u raw) i) u_i)))

    ;; return data structure
    raw))
	  


;;; this function opens, reads and parses an MPS file
;;; it returns an lp-data structure
(defun load-from-mps (mps-full-file-name)
  (let ((data)
	(header)
	(rows)
	(columns)
	(rhss)
	(bounds)
	(n-named 1)
	(n-slack 0)
	(m 0))

    (with-open-file (mps-stream (make-pathname :directory "/" 
				      :name mps-full-file-name) 
		       :direction :input)

      (read-mps-block mps-stream header "ROWS")
      (read-mps-block mps-stream rows "COLUMNS")
      (read-mps-block mps-stream columns "RHS")
      (read-mps-block mps-stream rhss "BOUNDS")
      (read-mps-block mps-stream bounds "ENDATA"))

    ;; gross error check
    (unless (and header rows columns)
      (error "corrupted MPS file"))

    ;; parse the header
    (setf data (make-mps)
	  header (nreverse header))
    (setf (mps-lp-name data) (cadar header))
    (pop header)
    (let ((obj-sense -1)
	  (obj-name))
      (loop
	 (unless header
	   (return))
	 (when (string= (caar header) "OBJNAME")
	   (setf obj-name (caadr header)))
	 (when (string= (caar header) "OBJSENSE")
	   (when (or (string= (caadr header) "MAX")
		     (string= (caadr header) "MAXIMIZE"))
	     (setf obj-sense 1)))
	 (pop header)
	 (pop header))
      (setf (mps-obj-spec data) (cons obj-name obj-sense)))

    ;; parse the rows
    (dolist (row rows)
      (cond ((string= (car row) "L")
	     (incf n-slack)
	     (incf m)
	     (push (cons (cadr row) '<=) (mps-row-spec data)))
	    ((string= (car row) "G")
	     (incf n-slack)
	     (incf m)
	     (push (cons (cadr row) '>=) (mps-row-spec data)))
	    ((string= (car row) "E")
	     (incf m)
	     (push (cons (cadr row) '=) (mps-row-spec data)))
	    ((string= (car row) "N")
	     (unless (car (mps-obj-spec data))
	       (setf (mps-obj-spec data) (cons (cadr row) (cdr (mps-obj-spec data))))))
	    (t
	     (print row)
	     (error "parse error in MPS file, row section"))))

    ;; parse the columns
    (let ((current-col-name (caar columns))
	  (current-col-coef (make-coef-alist (cdar columns))))
      (dolist (col (cdr columns))
	(cond ((string/= current-col-name (car col))
	       (incf n-named)
	       (push (cons current-col-name current-col-coef)
		     (mps-col-spec data))
	       (setf current-col-name (car col)
		     current-col-coef (make-coef-alist (cdr col))))
	      (t 
	       (nconc current-col-coef (make-coef-alist (cdr col))))))
      (push (cons current-col-name current-col-coef)
	    (mps-col-spec data)))

    ;; parse the right-hand sides
    (dolist (rhs rhss)
      (let ((length (length rhs)))
	(when (or (= length 3) 
		  (= length 5))
	  (decf length)
	  (pop rhs))
	(cond ((= length 2)
	       (destructuring-bind (n1 v1) rhs
		 (push (cons n1 (rationalize (read-from-string v1))) 
		       (mps-row-rhss data))))
	      ((= length 4)
	       (destructuring-bind (n1 v1 n2 v2) rhs
		 (push (cons n2 (rationalize (read-from-string v2))) 
		       (mps-row-rhss data))
		 (push (cons n1 (rationalize (read-from-string v1))) 
		       (mps-row-rhss data))))
	      (t
	       (error "parse error in MPS file, rhs section")))))

    ;; parse the bounds, if they exist
    (dolist (bound bounds)
      (let ((type (car bound))
	    (name (caddr bound))
	    (val (cadddr bound)))
	(cond ((string= type "UP")
	       (push (cons name (rationalize (read-from-string val))) 
		     (mps-col-ubounds data)))
	      ((string= type "LO")
	       (push (cons name (rationalize (read-from-string val))) 
		     (mps-col-lbounds data)))
	      ((string= type "FX")
	       (push (cons name (rationalize (read-from-string val))) 
		     (mps-col-ubounds data))
	       (push (cons name (rationalize (read-from-string val))) 
		     (mps-col-lbounds data)))
	      ((string= type "FR")
	       (push (cons name nil) (mps-col-ubounds data))
	       (push (cons name nil) (mps-col-lbounds data)))
	      (t
	       (error "parse error in MPS file, bounds section")))))

    (mps->raw data n-named n-slack m)))
	     

  
	    



      
	
	    

      
      

