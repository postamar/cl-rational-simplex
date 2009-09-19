(in-package :rationalsimplex)

;;;; Data structures for everything regarding the LP
;;;;


;;;;
(define-adjustable-vector fixnum)
(declaim (inline make-empty-adjvector-fixnum))
(defun make-empty-adjvector-fixnum ()
  (make-adjvector-fixnum 0 0))



;;;;
(defstruct row
  (name          ""  :type string)
  (ref           -1  :type fixnum)
  (is-active     nil   :type boolean)
  (slack-col-ref -1  :type fixnum)
  (coef          1   :type rational)
  (col-refs      (make-empty-adjvector-fixnum) :type adjvector-fixnum)
  (col-indices   (make-empty-adjvector-fixnum) :type adjvector-fixnum))



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
  (hsv          (error "column constructor") :type hsv))



(define-adjustable-vector integer)
(declaim (inline make-empty-adjvector-integer))
(defun make-empty-adjvector-integer ()
  (make-adjvector-integer 0 0))

(define-adjustable-vector column)
(define-adjustable-vector row)

;;;; 
(symbol-macrolet
    ((err (error "lp constructor")))
  (defstruct lp
    (name            ""  :type string)
    (is-infeasible   nil :type boolean)
    (is-unbounded    nil :type boolean)
    (obj-name        ""  :type string)
    (obj-sense       -1  :type fixnum)
    (columns         err :type adjvector-column)
    (rows            err :type adjvector-row)
    (active-row-refs err :type adjvector-fixnum)
    (active-col-refs err :type adjvector-fixnum)
    (active-row-inds err :type adjvector-fixnum)
    (active-col-inds err :type adjvector-fixnum)
    (col-ref-by-name (error "lp requires hashtable")
		     :type hash-table)
    (row-ref-by-name (error "lp requires hashtable") 
		     :type hash-table)))
  



;;;; Ratio in column
(defun rational-in-column (col index)
  (* (hsv-coef (column-hsv col))
     (aref (hsv-vis (column-hsv col)) index)))



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
  (let ((ind (adjvector-fixnum-ref (lp-active-col-inds lp) col-ref))
	(last-ind (- (adjvector-fixnum-fill-pointer (lp-active-col-refs lp)) 1)))
    (setf (column-is-active (adjvector-column-ref (lp-columns lp) col-ref)) nil
	  (adjvector-fixnum-ref (lp-active-col-inds lp) col-ref) -1)
    (unless (= ind last-ind)
      (let ((last-col-ref (adjvector-fixnum-ref (lp-active-col-refs lp) last-ind)))
	(setf (adjvector-fixnum-ref (lp-active-col-refs lp) ind) last-col-ref)
	(setf (adjvector-fixnum-ref (lp-active-col-inds lp) last-col-ref) ind)))
    (adjvector-fixnum-pop (lp-active-col-refs lp))))

(defun lp-remove-row (lp row-ref)
  (let ((ind (adjvector-fixnum-ref (lp-active-row-inds lp) row-ref))
	(last-ind (- (adjvector-fixnum-fill-pointer (lp-active-row-refs lp)) 1)))
    (setf (row-is-active (adjvector-row-ref (lp-rows lp) row-ref)) nil
	  (adjvector-fixnum-ref (lp-active-row-inds lp) row-ref) -1)
    (unless (= ind last-ind)
      (let ((last-row-ref (adjvector-fixnum-ref (lp-active-row-refs lp) last-ind)))
	(setf (adjvector-fixnum-ref (lp-active-row-refs lp) ind) last-row-ref)
	(setf (adjvector-fixnum-ref (lp-active-row-inds lp) last-row-ref) ind)))
    (adjvector-fixnum-pop (lp-active-row-refs lp))))



;;;; LP Constructor
(defun mps->lp (mps-data)
  (let*  ((n       (length (mps-col-spec mps-data)))
	  (m       (length (mps-row-spec mps-data)))
	  (cibyn   (make-hash-table :test 'equal))
	  (ribyn   (make-hash-table :test 'equal))
	  (arref   (make-empty-adjvector-fixnum))
	  (acref   (make-empty-adjvector-fixnum))
	  (arind   (make-empty-adjvector-fixnum))
	  (acind   (make-empty-adjvector-fixnum))
	  (cols    (make-adjvector-column (make-column :hsv (make-hsv)) 0))
	  (rows    (make-adjvector-row (make-row) 0)))
    ;; make column array
    (let ((j 0))
      (dolist (e (mps-col-spec mps-data))
	(adjvector-column-push-extend (make-column 
				       :name (car e)
				       :ref j
				       :is-active t
				       :hsv (make-hsv))
				      cols)
	(setf (gethash (car e) cibyn) j)
	(adjvector-fixnum-push-extend j acref)
	(adjvector-fixnum-push-extend j acind)
	(incf j))
      (dotimes (i m)
	(let ((name (concatenate 'string "slack" (princ-to-string i)))
	      (slack-hsv (make-hsv)))
	  (hsv-add i 1 slack-hsv)
	  (adjvector-column-push-extend 
	   (make-column :name name
			:ref j
			:is-active t
			:is-slack t
			:has-l t
			:has-u t
			:hsv slack-hsv)
	   cols)
	  (setf (gethash name cibyn) j))
	(adjvector-fixnum-push-extend j acref)
	(adjvector-fixnum-push-extend j acind)
	(incf j))
      (assert (= j (+ n m))))
    ;; make row array
    (let ((i 0))
      (dolist (e (mps-row-spec mps-data))
	(adjvector-fixnum-push-extend i arind)
	(adjvector-fixnum-push-extend i arref)
	(let ((row (make-row :name (car e) 
			     :ref i
			     :is-active t
			     :col-refs (make-empty-adjvector-fixnum)
			     :col-indices (make-empty-adjvector-fixnum))))
	  (setf (row-slack-col-ref row) (+ n i))
	  (adjvector-row-push-extend row rows))
	(setf (gethash (car e) ribyn) i)
	(incf i))
      (assert (= i m)))
    ;; fill columns
    (dolist (e (mps-col-spec mps-data))
      (let* ((j (gethash (car e) cibyn))
	     (col_j (adjvector-column-ref cols j))
	     (obj-name (car (mps-obj-spec mps-data)))
	     (alist '()))
	(dolist (coef (cdr e))
	  (destructuring-bind (row-name . a_ij) coef
	    (if (string= row-name obj-name)
		(setf (column-c col_j) a_ij)
		(push (cons (gethash row-name ribyn) a_ij) alist))))
	(multiple-value-bind (coef val-alist)
	    (factorize-ratio-alist (sort alist #'< :key #'car))
	  (dolist (pair val-alist) 
	    (destructuring-bind (i . v_i) pair
	      (let ((row (adjvector-row-ref rows i)))
		(adjvector-fixnum-push-extend j (row-col-refs row))
		(adjvector-fixnum-push-extend (hsv-length (column-hsv col_j)) (row-col-indices row)))
	      (hsv-add i (* (row-coef (adjvector-row-ref rows i)) v_i) (column-hsv col_j))))
	  (setf (hsv-coef (column-hsv col_j)) coef))))
    ;; fill slack refs in rows
    (dotimes (i m)
      (let ((row (adjvector-row-ref rows i)))
	(adjvector-fixnum-push-extend (row-slack-col-ref row) (row-col-refs row))
	(adjvector-fixnum-push-extend 0 (row-col-indices row))))
    ;; fill column bounds
    (dolist (lbound (mps-col-lbounds mps-data))
      (let ((l_j (cdr lbound))
	    (col_j (adjvector-column-ref cols (gethash (car lbound) cibyn))))
	(if l_j
	    (setf (column-has-l col_j) t
		  (column-l col_j) l_j)
	    (setf (column-has-l col_j) nil))))
    (dolist (ubound (mps-col-ubounds mps-data))
      (let ((u_j (cdr ubound))
	    (col_j (adjvector-column-ref cols (gethash (car ubound) cibyn))))
	(if u_j
	    (setf (column-has-u col_j) t
		  (column-u col_j) u_j)
	    (setf (column-has-u col_j) nil))))
    ;; fill slack column bounds
    (dolist (e (mps-row-rhss mps-data))
      (let* ((row   (adjvector-row-ref rows (gethash (car e) ribyn)))
	     (slack (adjvector-column-ref cols (row-slack-col-ref row)))
	     (val   (- (cdr e))))
	(setf (column-l slack) val
	      (column-u slack) val)))
    (dolist (e (mps-row-spec mps-data))
      (let* ((row      (adjvector-row-ref rows (gethash (car e) ribyn)))
	     (slack    (adjvector-column-ref cols (row-slack-col-ref row)))
	     (relation (cdr e)))
	(cond ((eq relation '>=)
	       (setf (column-has-l slack) nil))
	      ((eq relation '<=)
	       (setf (column-has-u slack) nil))
	      (t (assert (eq relation '=))))))
    (when (mps-row-lhss mps-data)
      (dolist (e (mps-row-lhss mps-data))
	(let* ((row   (adjvector-row-ref rows (gethash (car e) ribyn)))
	       (slack (adjvector-column-ref cols (row-slack-col-ref row)))
	       (val   (cdr e)))
	  (cond ((and (column-has-l slack) (column-has-u slack))
		 ;; row type E
		 (if (< val 0)
		     (incf (column-u slack) (abs val))
		     (decf (column-l slack) (abs val))))
		((column-has-u slack)
		 ;; row type G
		 (setf (column-has-l slack) t)
		 (decf (column-l slack) (abs val)))
		((column-has-l slack)
		 ;; row type L
		 (setf (column-has-u slack) t)
		 (incf (column-u slack) (abs val)))))))
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


;;;;
(declaim (inline lp-get-cost))
(defun lp-get-cost (lp j)
  (* (- (lp-obj-sense lp)) (column-c (adjvector-column-ref (lp-columns lp) j))))


;;;; Column iteration macro
(defmacro docol ((lp 
		  col
		  &key 
		  (index (gensym)) 
		  (row-ref (gensym))
		  (row (gensym)))
		 &body body)
  (let ((n-indices (gensym)))
    `(let ((,n-indices (hsv-length (column-hsv ,col))))
       (dotimes (,index ,n-indices)
	 (let* ((,row-ref (aref (hsv-is (column-hsv ,col)) ,index))
		(,row (adjvector-row-ref (lp-rows ,lp) ,row-ref)))
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
    `(let ((,n (adjvector-fixnum-fill-pointer (row-col-refs ,row))))
       (dotimes (,k ,n)
	 (let* ((,col-ref (adjvector-fixnum-ref (row-col-refs ,row) ,k))
		(,index (adjvector-fixnum-ref (row-col-indices ,row) ,k))
		(,col (adjvector-column-ref (lp-columns ,lp) ,col-ref)))
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
			  ((,a (aref (hsv-vis (column-hsv ,col)) ,index)))
			,@body))
		    body))))))))



