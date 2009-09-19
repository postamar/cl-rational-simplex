;;;; MPS file parser

;;; data structure containing everything in the MPS file
;;; in a non-awkward manner
(defstruct mps 
  lp-name
  obj-spec
  row-spec
  row-rhss
  row-lhss
  col-lbounds
  col-ubounds
  col-spec)

   

;;; trims whitespace left and right of string
(defun whitespace-trim (line)
  (let ((l 0)
	(r (length line)))
    (dotimes (i (length line))
      (if (or (char-equal (char line l) #\Space)
		  (char-equal (char line l) #\Tab))
	  (incf l)
	  (return)))
    (dotimes (i (length line))
      (if (or (char-equal (char line (- r 1)) #\Space)
	      (char-equal (char line (- r 1)) #\Tab))
	  (decf r)
	  (return)))
    (if (<= r l)
	""
	(subseq line l r))))
	
  

;;; scans one line from the MPS file into a list
(defun listify (line)
  (cond ((zerop (length line))
	 '())
	((not (or (char-equal (char line 0) #\Space)
		  (char-equal (char line 0) #\Tab)))
	 (if (and (<= 4 (length line))
		  (string-equal "NAME" (subseq line 0 4)))
	     (list "NAME" (subseq line 14 (min 22 (length line))))
	     (dotimes (i (length line) (list line))
	       (when (or (char-equal (char line i) #\Space)
			 (char-equal (char line i) #\Tab))
		 (return (list (subseq line 0 i)))))))
	(t 
	 (let ((linelist '()))
	   (when (<= 0 (length line))
	     (push (subseq line 1 (min 3 (length line))) linelist)
	     (when (<= 4 (length line))
	       (push (subseq line 4 (min 12 (length line))) linelist)
	       (when (<= 14 (length line))
		 (push (subseq line 14 (min 22 (length line))) linelist)
		 (when (<= 24 (length line))
		   (push (subseq line 24 (min 36 (length line))) linelist)
		   (when (<= 39 (length line))
		     (push (subseq line 39 (min 47 (length line))) linelist)
		     (when (<= 49 (length line))
		       (push (subseq line 49 (min 61 (length line))) linelist)))))))
	   (nreverse
	    (loop for elt in linelist
	       unless (string-equal "" (whitespace-trim elt))
	       collect (whitespace-trim elt)))))))


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
(defmacro read-mps-block (mps-stream block end-tags)
  (let ((line-list (gensym)))
    `(loop
	(let ((,line-list (get-next-line-list ,mps-stream)))
	  (unless ,line-list
	    (setf ,block nil)
	    (return))
	  (when (or ,@(mapcar 
		       #'(lambda (tag) 
			   `(string= ,tag (car ,line-list)))
		       (cons "ENDATA" end-tags)))
	    (return (car ,line-list)))
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



;;; this function opens, reads and parses an MPS file
;;; it returns an mps structure
(defun load-from-mps (mps-full-file-name)
  (let ((data)
	(header)
	(rows)
	(columns)
	(rhss)
	(lhss)
	(bounds)
	(n-named 1)
	(n-slack 0)
	(m 0))

    (with-open-file (mps-stream (make-pathname :directory "/" 
				      :name mps-full-file-name) 
		       :direction :input)

      (read-mps-block mps-stream header ("ROWS"))
      (read-mps-block mps-stream rows ("COLUMNS"))
      (read-mps-block mps-stream columns ("RHS")) 
      (when (string= "RANGES" (read-mps-block mps-stream rhss ("BOUNDS" "RANGES")))
	(read-mps-block mps-stream lhss ("BOUNDS")))
      (read-mps-block mps-stream bounds nil))

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


    ;; parse the ranges, if they exist
    (when lhss
      (dolist (lhs lhss)
	(let ((length (length lhs)))
	  (when (or (= length 3) 
		    (= length 5))
	    (decf length)
	    (pop lhs))
	  (cond ((= length 2)
		 (destructuring-bind (n1 v1) lhs
		   (push (cons n1 (rationalize (read-from-string v1))) 
			 (mps-row-lhss data))))
		((= length 4)
		 (destructuring-bind (n1 v1 n2 v2) lhs
		   (push (cons n2 (rationalize (read-from-string v2))) 
			 (mps-row-lhss data))
		   (push (cons n1 (rationalize (read-from-string v1))) 
			 (mps-row-lhss data))))
		(t
		 (error "parse error in MPS file, range section"))))))


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
	      ((string= type "PL")
	       (push (cons name nil) (mps-col-ubounds data)))
	      ((string= type "MI")
	       (push (cons name nil) (mps-col-lbounds data)))
	      (t
	       (error "parse error in MPS file, bounds section")))))
    (values data n-named n-slack m)))
	     

  
	    



      
	
	    

      
      

