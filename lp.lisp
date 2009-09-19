;;;; Data structures for everything regarding the LP
;;;;

;;;; 
(defstruct lp
  (name            ""  :type string)
  (is-infeasible   nil :type boolean)
  (is-unbounded    nil :type boolean)
  (obj-name        ""  :type string)
  (obj-sense       -1  :type fixnum)
  (columns         #() :type vector)
  (rows            #() :type vector)
  (active-row-refs #() :type vector)         
  (active-col-refs #() :type vector)
  (active-row-inds #() :type vector)
  (active-col-inds #() :type vector)
  (col-ref-by-name (error "lp requires hashtable")
		   :type hash-table)
  (row-ref-by-name (error "lp requires hashtable") 
		   :type hash-table))



;;;;  
(defstruct row
  (name          ""  :type string)
  (ref           -1  :type fixnum)
  (is-active     nil   :type boolean)
  (relation      '=  :type symbol)
  (slack-col-ref -1  :type fixnum)
  (b             0   :type rational)
  (coef          1   :type rational)
  (col-refs      #() :type vector)
  (col-indices   #() :type vector))

  

;;;;
(defstruct column
  (name         ""  :type string)
  (ref          -1  :type fixnum)
  (is-active    nil :type boolean)
  (is-slack     nil :type boolean)
  (c            0   :type rational)
  (has-l        t   :type boolean)
  (has-u        nil :type boolean)
  (l            0   :type rational)
  (u            0   :type rational)
  (coef         1   :type rational)
  (row-refs     #() :type vector)
  (values       #() :type vector))



;;;; Ratio in column
(defun rational-in-column (col index)
  (* (column-coef col)
     (aref (column-values col) index)))



;;;; Arithmetic functions
(defun factorize-ratio-pair (f1 nl1 f2 nl2)
  (let ((d1   (denominator f1))
	(d2   (denominator f2))
	(n1   (numerator f1))
	(n2   (numerator f2)))
    (let* ((d    (* d1 d2))
	   (n1d2 (* n1 d2))
	   (n2d1 (* n2 d1))
	   (n    (gcd n1d2 n2d1))
	   (f1   (/ n1d2 n))
	   (f2   (/ n2d1 n)))
      (values (/ n d) 
	      (nconc (mapcar #'(lambda (v1) (cons (car v1) (* f1 (cdr v1)))) nl1)
		     (mapcar #'(lambda (v2) (cons (car v2) (* f2 (cdr v2)))) nl2))))))

(defun factorize-ratio-alist (l)
  (let* ((len (length l))
	 (h (ash len -1)))
    (if (= len 1)
	(values (abs (cdar l)) (list (cons (caar l) (signum (cdar l)))))
	(multiple-value-bind (f1 nl1) 
	    (factorize-ratio-alist (subseq l 0 h))
	  (multiple-value-bind (f2 nl2) 
	      (factorize-ratio-alist (subseq l h))
	    (multiple-value-bind (f nl) 
		(factorize-ratio-pair f1 nl1 f2 nl2)
	      (values f nl)))))))



;;;; LP modification functions
(defun lp-remove-column (lp col-ref)
  (let ((ind (aref (lp-active-col-inds lp) col-ref))
	(last-ind (- (length (lp-active-col-refs lp)) 1)))
    (setf (column-is-active (aref (lp-columns lp) col-ref)) nil
	  (aref (lp-active-col-inds lp) col-ref) -1)
    (unless (= ind last-ind)
      (let ((last-col-ref (aref (lp-active-col-refs lp) last-ind)))
	(setf (aref (lp-active-col-refs lp) ind) last-col-ref)
	(setf (aref (lp-active-col-inds lp) last-col-ref) ind)))
    (vector-pop (lp-active-col-refs lp))))

(defun lp-remove-row (lp row-ref)
  (let ((ind (aref (lp-active-row-inds lp) row-ref))
	(last-ind (- (length (lp-active-row-refs lp)) 1)))
    (setf (row-is-active (aref (lp-rows lp) row-ref)) nil
	  (aref (lp-active-row-inds lp) row-ref) -1)
    (unless (= ind last-ind)
      (let ((last-row-ref (aref (lp-active-row-refs lp) last-ind)))
	(setf (aref (lp-active-row-refs lp) ind) last-row-ref)
	(setf (aref (lp-active-row-inds lp) last-row-ref) ind)))
    (vector-pop (lp-active-row-refs lp))))



;;;; LP Constructor
(defun mps->lp (mps-data)
  (let  ((n       (length (mps-col-spec mps-data)))
	 (n-slack 0)
	 (m       (length (mps-row-spec mps-data)))
	 (cibyn   (make-hash-table :test 'equal))
	 (ribyn   (make-hash-table :test 'equal))
	 (arref   (make-vector fixnum))
	 (acref   (make-vector fixnum))
	 (arind   (make-vector fixnum))
	 (acind   (make-vector fixnum))
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
    (let ((j 0))
      (dolist (e (mps-col-spec mps-data))
	(vector-push-extend j acind)
	(vector-push-extend j acref)
	(vector-push-extend (make-column 
			     :name (car e)
			     :ref j
			     :is-active t)
			    cols)
	(setf (gethash (car e) cibyn) j)
	(incf j)))
    ;; make row array
    (let ((i 0))
      (dolist (e (mps-row-spec mps-data))
	(vector-push-extend i arind)
	(vector-push-extend i arref)
	(let ((row (make-row :name (car e) 
			     :ref i
			     :is-active t
			     :relation (cdr e)
			     :col-refs (make-vector fixnum)
			     :col-indices (make-vector fixnum))))
	  (unless (eq '= (cdr e))
	    (let ((slack-col-ref (+ n n-slack)))
	      (setf (row-slack-col-ref row) slack-col-ref)
	      (vector-push-extend slack-col-ref (row-col-refs row))
	      (vector-push-extend 0 (row-col-indices row))
	      (incf n-slack)))
	  (vector-push-extend row rows))
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
	(setf (row-b row)    b_i
	      (row-coef row) rfact)))
    ;; add slack columns
    (let ((j 0))
      (dotimes (i m)
	(let ((row_i (aref rows i)))
	  (unless (eq '= (row-relation row_i))
	    (let ((val (* (if (eq '<= (row-relation row_i)) 1 -1)
			  (row-coef row_i)))
		  (name (concatenate 'string 
				     "slack" 
				     (princ-to-string j))))
	      (let ((col-ref (length cols)))
		(vector-push-extend 
		 (make-column :name name
			      :ref col-ref
			      :is-active t
			      :is-slack t
			      :row-refs (make-array 1 
						    :initial-element i 
						    :element-type 'fixnum)
			      :values (make-array 1
						  :initial-element val
						  :element-type 'integer))
		 cols)
		(setf (gethash name cibyn) col-ref)
		(vector-push-extend col-ref acref)
		(vector-push-extend col-ref acind)))
	    (incf j)))))
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
	(multiple-value-bind (coef val-alist)
	    (factorize-ratio-alist (sort alist #'< :key #'car))
	  (let ((row-refs (make-vector fixnum))
		(values   (make-vector integer)))
	    (dolist (pair val-alist) 
	      (destructuring-bind (i . v_i) pair
		(let ((row (aref rows i)))
		  (vector-push-extend j (row-col-refs row))
		  (vector-push-extend (length row-refs) (row-col-indices row)))
		(vector-push-extend i row-refs)
		(vector-push-extend (* (row-coef (aref rows i)) v_i)
				    values)))
	    (setf (column-coef col_j) coef
		  (column-row-refs col_j) row-refs
		  (column-values col_j) values)))))
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
    (make-lp :name            (mps-lp-name mps-data)
	     :obj-name        (car (mps-obj-spec mps-data))
	     :obj-sense       (cdr (mps-obj-spec mps-data))
	     :columns         cols
	     :rows            rows
	     :active-row-refs arref
	     :active-col-refs acref
	     :active-row-inds arind
	     :active-col-inds acind
	     :col-ref-by-name cibyn
	     :row-ref-by-name ribyn)))



;;;; Column iteration macro
(defmacro docol ((lp 
		  col
		  &key 
		  (index (gensym)) 
		  (row-ref (gensym))
		  (row (gensym)))
		 &body body)
  (let ((n-indices (gensym)))
    `(let ((,n-indices (length (column-row-refs ,col))))
       (dotimes (,index ,n-indices)
	 (let* ((,row-ref (aref (column-row-refs ,col) ,index))
		(,row (aref (lp-rows ,lp) ,row-ref)))
	   (when (row-is-active ,row)
	     ,@body))))))



;;;; Row iteration macros
(defmacro dorow ((lp 
		  row
		  &key 
		  (visit-only-fixed nil)
		  (visit-also-fixed nil)
		  (visit-slack t)
		  (col-ref (gensym)) 
		  (index (gensym)) 
		  (col (gensym))
		  (a nil)
		  (has-l nil)
		  (has-u nil)
		  (l nil)
		  (u nil))
		 &body body)
  (let ((n (gensym))
	(k (gensym)))
    `(let ((,n (length (row-col-refs ,row))))
       (dotimes (,k ,n)
	 (let* ((,col-ref (aref (row-col-refs ,row) ,k))
		(,index (aref (row-col-indices ,row) ,k))
		(,col (aref (lp-columns ,lp) ,col-ref)))
	   (symbol-macrolet ,(remove nil (list 
				       (when has-l `(,has-l (column-has-l ,col)))
				       (when has-u `(,has-u (column-has-u ,col)))
				       (when l `(,l (column-l ,col)))
				       (when u `(,u (column-u ,col)))))
	     (unless ,(if visit-only-fixed
			  `(not (column-is-active ,col))
			  `(or ,(unless visit-also-fixed `(not (column-is-active ,col)))
			       ,(unless visit-slack `(column-is-slack ,col))))
	       ,@(if a
		    `((symbol-macrolet
			 ((,a (aref (column-values ,col) ,index)))
			,@body))
		    body))))))))
