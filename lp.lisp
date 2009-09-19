;;;; Data structures for everything regarding the LP
;;;;

;;;; 
(defstruct lp
  (name                    ""           :type string)
  (is-infeasible           nil          :type boolean)
  (is-unbounded            nil          :type boolean)
  (obj-name                ""           :type string)
  (obj-sense               -1           :type fixnum)
  (columns                 #()          :type vector)
  (rows                    #()          :type vector)
  (column-indices-by-names (error "lp requires hashtable")
			   :type hash-table)
  (row-indices-by-names    (error "lp requires hashtable") 
			   :type hash-table))



;;;;  
(defstruct row
  (name         ""  :type string)
  (is-removed   t   :type boolean)
  (relation     '=  :type symbol)
  (slack-index  -1  :type fixnum)
  (b            0   :type rational)
  (squash-coef  1   :type rational))

  

;;;;
(defstruct column
  (name         ""  :type string)
  (is-removed   t   :type boolean)
  (is-slack     nil :type boolean)
  (c            0   :type rational)
  (has-l        t   :type boolean)
  (has-u        nil :type boolean)
  (l            0   :type rational)
  (u            0   :type rational)
  (squash-coef  1   :type rational)
  (indices      #() :type vector)
  (values       #() :type vector))



;;;; LP Constructor
(defun mps->lp (mps-data)
  (let  ((n       (length (mps-col-spec mps-data)))
	 (n-slack 0)
	 (m       (length (mps-row-spec mps-data)))
	 (cibyn   (make-hash-table :test 'equal))
	 (ribyn   (make-hash-table :test 'equal))
	 (cols    (make-array 0
			      :initial-element (make-column)
			      :adjustable t 
			      :fill-pointer t 
			      :element-type 'column))
	 (rows    (make-array 0
			      :initial-element (make-row)
			      :adjustable t
			      :fill-pointer t
			      :element-type 'row)))
    ;; make column array
    (let ((i 0))
      (dolist (e (mps-col-spec mps-data))
	(vector-push-extend (make-column 
			     :name (car e)
			     :is-removed nil)
			    cols)
	(setf (gethash (car e) cibyn) i)
	(incf i)))
    ;; make row array
    (let ((i 0))
      (dolist (e (mps-row-spec mps-data))
	(vector-push-extend (make-row 
			     :name (car e) 
			     :is-removed nil 
			     :relation (cdr e)
			     :slack-index (if (eq '= (cdr e)) -1 (+ n n-slack)))
			    rows)
	(unless (eq '= (cdr e))
	  (incf n-slack))
	(setf (gethash (car e) ribyn) i)
	(incf i)))
    ;; fill rows and check for non-negativity
    (dolist (e (mps-row-rhss mps-data))
      (let ((row   (aref rows (gethash (car e) ribyn)))
	    (b_i   (cdr e))
	    (rfact 1))
	(when (< b_i 0)
	  (setf b_i   (- b_i)
		rfact -1))
	(setf (row-b row)            b_i
	      (row-squash-coef row)  rfact)))
    ;; add slack columns
    (let ((j 0))
      (dotimes (i m)
	(let ((row_i (aref rows i)))
	  (unless (eq '= (row-relation row_i))
	    (let ((val (* (if (eq '<= (row-relation row_i)) 1 -1)
			  (row-squash-coef row_i)))
		  (name (concatenate 'string 
				     "slack" 
				     (princ-to-string j))))
	      (setf (gethash name cibyn) (length cols))
	      (vector-push-extend 
	       (make-column :name name
			    :is-removed nil
			    :is-slack t
			    :squash-coef (row-squash-coef row_i)
			    :indices (make-array 1 
						 :initial-element i 
						 :element-type 'fixnum)
			    :values (make-array 1
						:initial-element val
						:element-type 'rational))
	       cols)
	      (incf j))))))
    ;; fill columns
    (dolist (e (mps-col-spec mps-data))
      (let* ((j (gethash (car e) cibyn))
	     (col_j (aref cols j))
	     (obj-name (car (mps-obj-spec mps-data)))
	     (alist '()))
	(dolist (coef (cdr e))
	  (destructuring-bind (row-name . a_ij) coef
	    (if (string= row-name obj-name)
		(setf (column-c col_j) a_ij)
		(push (cons (gethash row-name ribyn) a_ij) alist))))
	(setf alist (sort alist #'< :key #'car))
	(let* ((counter 0)
	       (alist-length (length alist))
	       (indices (make-array alist-length 
				    :initial-element -1 
				    :element-type 'fixnum))
	       (values (make-array alist-length 
				   :initial-element 0 
				   :element-type 'rational)))
	  (dolist (pair alist) 
	    (destructuring-bind (i . a_ij) pair
	      (setf (aref indices counter) i
		    (aref values counter) (* (row-squash-coef (aref rows i)) 
					     a_ij))
	      (incf counter)))
	  (setf (column-indices col_j) indices
		(column-values col_j) values))))
    ;; fill column bounds
    (dolist (lbound (mps-col-lbounds mps-data))
      (let ((l_j (cdr lbound))
	    (col_j (aref cols (gethash (car lbound) cibyn))))
	(if l_j
	    (setf (column-has-l col_j) t
		  (column-l col_j) l_j)
	    (setf (column-has-l col_j) nil))))
    (dolist (ubound (mps-col-ubounds mps-data))
      (let ((u_j (cdr ubound))
	    (col_j (aref cols (gethash (car ubound) cibyn))))
	(if u_j
	    (setf (column-has-u col_j) t
		  (column-u col_j) u_j)
	    (setf (column-has-u col_j) nil))))
    ;; return lp instance
    (make-lp :name                    (mps-lp-name mps-data)
	     :obj-name                (car (mps-obj-spec mps-data))
	     :obj-sense               (cdr (mps-obj-spec mps-data))
	     :columns                 cols
	     :rows                    rows
	     :column-indices-by-names cibyn
	     :row-indices-by-names    ribyn)))



;;;; Column iteration macro
(defmacro docol ((lp 
		  col
		  &key 
		  (index (gensym)) 
		  (i (gensym))
		  (row (gensym)))
		 &body body)
  (let ((n-indices (gensym)))
    `(let ((,n-indices (length (column-indices ,col))))
       (dotimes (,index ,n-indices)
	 (let* ((,i (aref (column-indices ,col) ,index))
		(,row (aref (lp-rows ,lp) ,i)))
	   (unless (row-is-removed ,row)
	     ,@body))))))



;;;; Row iteration macros
(defmacro dorow ((lp 
		  i 
		  &key 
		  (visit-only-fixed nil)
		  (visit-also-fixed nil)
		  (visit-slack t)
		  (j (gensym)) 
		  (index (gensym)) 
		  (col (gensym))
		  (a nil)
		  (has-l nil)
		  (has-u nil)
		  (l nil)
		  (u nil))
		 &body body)
  (let ((n (gensym)))
    `(let ((,n (length (lp-columns ,lp))))
       (dotimes (,j ,n)
	 (let ((,col (aref (lp-columns ,lp) ,j)))
	   (symbol-macrolet ,(remove nil (list 
				       (when has-l `(,has-l (column-has-l ,col)))
				       (when has-u `(,has-u (column-has-u ,col)))
				       (when l `(,l (column-l ,col)))
				       (when u `(,u (column-u ,col)))))
	     (unless ,(if visit-only-fixed
			  `(not (column-is-removed ,col))
			  `(or ,(unless visit-also-fixed `(column-is-removed ,col))
			       ,(unless visit-slack `(column-is-slack ,col))))
	       (let ((,index (find-index (column-indices ,col) ,i)))
		 (symbol-macrolet 
		     ,(remove nil 
			      (list (when a `(,a (aref (column-values ,col) ,index)))))
		   (unless (= -1 ,index)
		     ,@body))))))))))
	     

