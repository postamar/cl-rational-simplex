(in-package :rationalsimplex)

;;;;; The big data structure macro file
;;;;; Implements the definition of the following:
;;;;; 
;;;;; ADJUSTABLE VECTORS
;;;;; similar to (make-array :adjustable t :fill-pointer t)
;;;;; necessary because of SBCL's inadequate type system 
;;;;; (I'm told SBCL maintainers' laziness and the lisp machines are to blame)
;;;;;
;;;;; PERSISTENT DATA STRUCTURES
;;;;; definition macros are available for:
;;;;; * stacks
;;;;; * linked lists
;;;;; * binary trees
;;;;; * threaded binary trees
;;;;; * splay trees


;;;;; In the rational simplex implementation, only the
;;;;; ADJUSTABLE VECTORS are used, liberally, and the
;;;;; STACKS and SPLAY TREES are used, on occasion
;;;;; Note that the hyper-sparse vector implementation
;;;;; is similar to that of the adj. vectors.



;;;; ADJUSTABLE VECTOR
;;;; definition macro, takes desired content type as arg
(defmacro define-adjustable-vector (type)
  (let* ((namestr (concatenate 'string "ADJVECTOR-" (princ-to-string type)))
	 (constr (concatenate 'string "%MAKE-" namestr))
	 (makefull (concatenate 'string "MAKE-" namestr))
	 (push (concatenate 'string namestr "-PUSH"))
	 (pushext (concatenate 'string namestr "-PUSH-EXTEND"))
	 (pop (concatenate 'string namestr "-POP"))
	 (adjvv (concatenate 'string namestr "-DATA"))
	 (adjvfp (concatenate 'string namestr "-FILL-POINTER"))
	 (adjva (concatenate 'string namestr "-PREV-ALLOC-COUNTER"))
	 (adjvb (concatenate 'string namestr "-CURRENT-ALLOC-COUNTER"))
	 (adjvref (concatenate 'string namestr "-REF"))
	 (setadjvref (concatenate 'string "SET-" namestr "-REF"))
	 (ielt (gensym "INIT-ELEMENT-"))
	 (elt (gensym "ELEMENT-"))
	 (adjv (gensym "ADJVECTOR-"))
	 (i (gensym "INDEX-"))
	 (inalloc (gensym "INIT-ALLOC-SIZE-"))
	 (a (gensym))
	 (b (gensym))
	 (a+b (gensym))
	 (v (gensym)))
    `(progn
       (defstruct (,(intern namestr)
		    (:constructor ,(intern constr)))
	 (fill-pointer 0 :type fixnum)
	 (prev-alloc-counter 0 :type fixnum)
	 (current-alloc-counter 0 :type fixnum)
	 (data (error "adjvector construction") 
	       :type (simple-array ,type 1)))
       (declaim (inline ,(intern makefull)))
       (defun ,(intern makefull) (,ielt ,inalloc)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,type ,ielt)
		  (fixnum ,inalloc))
	 (let ((,a 1)
	       (,b 1))
	   (declare (fixnum ,a ,b))
	   (loop 
	      (when (<= ,inalloc ,b)
		(return))
	      (rotatef ,a ,b)
	      (incf ,b ,a))
	   (,(intern constr)
	     :prev-alloc-counter ,a
	     :current-alloc-counter ,b
	     :data (make-array ,b :initial-element ,ielt :element-type ',type))))
       (declaim (inline ,(intern push)))
       (defun ,(intern push) (,elt ,adjv)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,type ,elt)
		  (,(intern namestr) ,adjv))
	 (setf (aref (,(intern adjvv) ,adjv) (,(intern adjvfp) ,adjv)) ,elt)
	 (incf (,(intern adjvfp) ,adjv)))
       (declaim (inline ,(intern pop)))
       (defun ,(intern pop) (,adjv)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,(intern namestr) ,adjv))
	 (aref (,(intern adjvv) ,adjv) (decf (,(intern adjvfp) ,adjv))))
       (declaim (inline ,(intern pushext)))
       (defun ,(intern pushext) (,elt ,adjv)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,type ,elt)
		  (,(intern namestr) ,adjv))
	 (if (<= (,(intern adjvb) ,adjv) (,(intern adjvfp) ,adjv))
	     (if (zerop (,(intern adjvb) ,adjv))
		 (setf (,(intern adjva) ,adjv) 1
		       (,(intern adjvb) ,adjv) 1
		       (,(intern adjvv) ,adjv) (make-array 1 :initial-element ,elt :element-type ',type)
		       (,(intern adjvfp) ,adjv) 1)
		 (progn 
		   (let* ((,a (,(intern adjva) ,adjv))
			  (,b (,(intern adjvb) ,adjv))
			  (,a+b (+ ,a ,b))
			  (,v (make-array ,a+b :initial-element ,elt :element-type ',type)))
		     (declare (fixnum ,a ,b ,a+b)
			      ((simple-array ,type 1) ,v))
		     (replace ,v (,(intern adjvv) ,adjv))
		     (setf (,(intern adjva) ,adjv) ,b
			   (,(intern adjvb) ,adjv) ,a+b
			   (,(intern adjvv) ,adjv) ,v)
		     (incf (,(intern adjvfp) ,adjv)))))
	     (,(intern push) ,elt ,adjv)))
       (declaim (inline ,(intern adjvref)))
       (defun ,(intern adjvref) (,adjv ,i)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,(intern namestr) ,adjv)
		  (fixnum ,i))
	 (assert (< -1 ,i (,(intern adjvfp) ,adjv)))
	 (aref (,(intern adjvv) ,adjv) ,i))
       (declaim (inline ,(intern setadjvref)))
       (defun ,(intern setadjvref) (,adjv ,i ,elt)
	 (declare (optimize (debug 0) (safety 0) (speed 3)))
	 (declare (,(intern namestr) ,adjv)
		  (fixnum ,i)
		  (,type ,elt))
	 (assert (< -1 ,i (,(intern adjvfp) ,adjv)))
	 (setf (aref (,(intern adjvv) ,adjv) ,i) ,elt))
       (defsetf ,(intern adjvref) ,(intern setadjvref))
       t)))



;;;; PERSISTENT DATA STRUCTURES
;;;; general definition macro
;;;; takes desired name, key and value types and initializers,
;;;; and number of pointers per element:
;;;; for example, a binary tree needs 2, a single linked list 1, and a stack 0
(defmacro defdatastruct (name keytype valtype npointer init-key-elt init-val-elt)
  (let ((constr (gensym))
	(namestr (princ-to-string name))
	(errmsg (concatenate 'string (princ-to-string name) " constructor"))
	(ds (gensym "DS"))
	(key (gensym "KEY"))
	(val (gensym "VAL"))
	(ptr (gensym "PTR"))
	(keyarr (gensym "KEYARR"))
	(valarr (gensym "VALARR"))
	(ptrarr (gensym "PTRARR"))
	(ptrargs (loop for k from 1 upto npointer collect (gensym "PTR")))
	(i (gensym "I")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app)))
	   (slotname (k)
	     (intern (concatenate 'string namestr "-POINTERS-" (princ-to-string k))))
	   (accname (k)
	     (intern (concatenate 'string namestr "-POINTER-" (princ-to-string k))))
	   (setname (k)
	     (intern (concatenate 'string "SET-" namestr "-POINTER-" (princ-to-string k)))))
      `(progn 
	 (defstruct (,name
		      (:constructor ,constr))
	   (keys     (error ,errmsg) :type (simple-array ,keytype 1))
	   (values   (error ,errmsg) :type (simple-array ,valtype 1))
	   ,@(loop for k from 1 upto npointer
		collect (list (intern (concatenate 'string 
						   "POINTERS-"
						   (princ-to-string k)))
			      `(error ,errmsg)
			      ':type
			      '(simple-array fixnum 1)))
	   (garbage  (error ,errmsg) :type (simple-array fixnum 1))
	   (data-length 0 :type fixnum)
	   (garbage-length 0 :type fixnum)
	   (data-prev-alloc-counter 1 :type fixnum)
	   (data-current-alloc-counter 1 :type fixnum)
	   (garbage-prev-alloc-counter 1 :type fixnum)
	   (garbage-current-alloc-counter 1 :type fixnum)
	   (header -1 :type fixnum))
	 (defun ,(fname "MAKE-" "") ()
	   (,constr 
	    :keys     (make-array 1 :initial-element ,init-key-elt :element-type ',keytype)
	    :values   (make-array 1 :initial-element ,init-val-elt :element-type ',valtype)
	    ,@(mapcan #'copy-list
		      (loop for k from 1 upto npointer
			 collect (list (intern (concatenate 'string 
							    "POINTERS-"
							    (princ-to-string k))
					       :keyword)
				       '(make-array 1 :initial-element -1 :element-type 'fixnum))))
	    :garbage  (make-array 1 :initial-element -1 :element-type 'fixnum)))
	 (declaim (inline ,(fname "" "-VALUE")))
	 (defun ,(fname "" "-VALUE") (,ds ,i)
	   (aref ( ,(fname "" "-VALUES") ,ds) ,i))
	 (declaim (inline ,(fname "" "-KEY")))
	 (defun ,(fname "" "-KEY") (,ds ,i)
	   (aref ( ,(fname "" "-KEYS") ,ds) ,i))
	 ,@(loop for k from 1 upto npointer
	      collect `(declaim (inline ,(accname k))))
	 ,@(loop for k from 1 upto npointer
	      collect `(defun ,(accname k)
			   (,ds ,i)
			 (aref (,(slotname k) ,ds) ,i)))
	 (declaim (inline ,(fname "SET-" "-VALUE")))
	 (defun ,(fname "SET-" "-VALUE") (,ds ,i ,val)
	   (setf (aref ( ,(fname "" "-VALUES") ,ds) ,i) ,val))
	 (declaim (inline ,(fname "SET-" "-KEY")))
	 (defun ,(fname "SET-" "-KEY") (,ds ,i ,key)
	   (setf (aref ( ,(fname "" "-KEYS") ,ds) ,i) ,key))
	 ,@(loop for k from 1 upto npointer
	      collect `(declaim (inline ,(setname k))))
	 ,@(loop for k from 1 upto npointer
	      collect `(defun ,(setname k) (,ds ,i ,ptr)
			 (setf (aref (,(slotname k) ,ds) ,i) ,ptr)))
	 (defun ,(fname "ADD-IN-" "") (,ds ,key ,val ,@ptrargs)
	   (if (zerop (,(fname "" "-GARBAGE-LENGTH") ,ds))
	       (let ((,i (,(fname "" "-DATA-LENGTH") ,ds)))
		 (if (<= (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds) ,i)
		     (progn 
		       (rotatef (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds)
				(,(fname "" "-DATA-PREV-ALLOC-COUNTER") ,ds))
		       (incf (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds)
			     (,(fname "" "-DATA-PREV-ALLOC-COUNTER") ,ds))
		       (let ((,keyarr (make-array (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds)
						  :initial-element ,key :element-type ',keytype))
			     (,valarr (make-array (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds)
						  :initial-element ,val :element-type ',valtype)))
			 (replace ,keyarr (,(fname "" "-KEYS") ,ds) :end2 ,i)
			 (replace ,valarr (,(fname "" "-VALUES") ,ds) :end2 ,i)
			 (setf (,(fname "" "-KEYS") ,ds) ,keyarr
			       (,(fname "" "-VALUES") ,ds) ,valarr))
		       ,@(loop for k from 0 below npointer 
			    collect 
			      `(let ((,ptrarr (make-array (,(fname "" "-DATA-CURRENT-ALLOC-COUNTER") ,ds)
							  :initial-element ,(nth k ptrargs) :element-type 'fixnum)))
				 (replace ,ptrarr (,(slotname (+ k 1)) ,ds) :end2 ,i)
				 (setf (,(slotname (+ k 1)) ,ds) ,ptrarr))))
		     (setf ,@(mapcan #'copy-list
				     (loop for k from 0 below npointer 
					collect (list `(aref (,(slotname (+ k 1)) ,ds) ,i)
						      (nth k ptrargs))))
			   (aref (,(fname "" "-KEYS") ,ds) ,i) ,key
			   (aref (,(fname "" "-VALUES") ,ds) ,i) ,val))
		 (incf (,(fname "" "-DATA-LENGTH") ,ds))
		 ,i)
	       (let ((,i (aref (,(fname "" "-GARBAGE") ,ds)
			       (decf (,(fname "" "-GARBAGE-LENGTH") ,ds)))))
		 (setf ,@(mapcan #'copy-list
				 (loop for k from 0 below npointer 
				    collect (list `(aref (,(slotname (+ k 1)) ,ds) ,i)
						  (nth k ptrargs))))
		       (aref (,(fname "" "-KEYS") ,ds) ,i) ,key
		       (aref (,(fname "" "-VALUES") ,ds) ,i) ,val)
		 ,i)))
	 (declaim (inline ,(fname "REMOVE-IN-" "")))
	 (defun ,(fname "REMOVE-IN-" "") (,ds ,i)
	   (cond ((<= (,(fname "" "-DATA-LENGTH") ,ds) ,i)
		  -1)
		 ((<= (,(fname "" "-GARBAGE-CURRENT-ALLOC-COUNTER") ,ds)
		   (,(fname "" "-GARBAGE-LENGTH") ,ds))
		  (rotatef (,(fname "" "-GARBAGE-CURRENT-ALLOC-COUNTER") ,ds)
			   (,(fname "" "-GARBAGE-PREV-ALLOC-COUNTER") ,ds))
		  (incf (,(fname "" "-GARBAGE-CURRENT-ALLOC-COUNTER") ,ds)
			(,(fname "" "-GARBAGE-PREV-ALLOC-COUNTER") ,ds))
		  (let ((,ptrarr (make-array (,(fname "" "-GARBAGE-CURRENT-ALLOC-COUNTER") ,ds)
					     :initial-element ,i :element-type 'fixnum)))
		    (replace ,ptrarr (,(fname "" "-GARBAGE") ,ds) 
			     :end2 (,(fname "" "-GARBAGE-LENGTH") ,ds))
		    (setf (,(fname "" "-GARBAGE") ,ds) ,ptrarr))
		  (incf (,(fname "" "-GARBAGE-LENGTH") ,ds))
		  ,i)
		 (t 
		  (setf (aref (,(fname "" "-GARBAGE") ,ds) 
			      (,(fname "" "-GARBAGE-LENGTH") ,ds)) ,i)
		  (incf (,(fname "" "-GARBAGE-LENGTH") ,ds))
		  ,i)))
	 (declaim (inline ,(fname "" "-COUNT")))
	 (defun ,(fname "" "-COUNT") (,ds)
	   (- (,(fname "" "-DATA-LENGTH") ,ds)
	      (,(fname "" "-GARBAGE-LENGTH") ,ds)))
	 (declaim (inline ,(fname "RESET-" "")))
	 (defun ,(fname "RESET-" "") (,ds)
	   (setf (,(fname "" "-GARBAGE-LENGTH") ,ds) 0
		 (,(fname "" "-DATA-LENGTH") ,ds) 0
		 (,(fname "" "-HEADER") ,ds) -1))
	 (defsetf ,(fname "" "-VALUE") ,(fname "SET-" "-VALUE"))
	 (defsetf ,(fname "" "-KEY") ,(fname "SET-" "-KEY"))
	 ,@(loop for k from 1 upto npointer
	      collect `(defsetf ,(accname k) ,(setname k)))
	 t))))


;;;; Helper macro for defining a DFS tree-mapping function
;;;; Is called only by data structure definition macros for bin trees
(defmacro defmaptree (dsname fn ds initial-i)
  (let* ((keystr (concatenate 'string dsname "-KEY"))
	 (valstr (concatenate 'string dsname "-VALUE"))
	 (lstr (concatenate 'string dsname "-POINTER-1"))
	 (rstr (concatenate 'string dsname "-POINTER-2"))
	 (dfs (gensym "FNDFS"))
	 (tree (gensym "TREE"))
	 (i (gensym "I"))
	 (li (gensym "I"))
	 (ri (gensym "I"))
	 (fnstr (princ-to-string fn))
	 (fncall (cond ((< (length fnstr) (length "#'"))
			`(,fn))
		       ((string= "#'" (subseq fnstr 0 2))
			`(funcall ,fn))
		       (t
			`(,fn)))))
    `(labels ((,dfs (,tree ,i)
		(let ((,li (,(intern lstr) ,tree ,i))
		      (,ri (,(intern rstr) ,tree ,i)))
		  (unless (= -1 ,li)
		    (,dfs ,tree ,li))
		  (,@fncall (,(intern keystr) ,tree ,i) (,(intern valstr) ,tree ,i))
		  (unless (= -1 ,ri)
		    (,dfs ,tree ,ri)))))
       (unless (= -1 ,initial-i)
	 (,dfs ,ds ,initial-i)))))
	 
	   
;;;; Helper macro for defining a general mapping function
;;;; Is called only by data structure definition macros for 
;;;; something other than bin trees
(defmacro defmap (dsname fn-firsti fn-nexti fn ds moreds)
  (let* ((dslist (append (list ds) moreds))
	 (nds (length dslist))
	 (keystr (concatenate 'string dsname "-KEY"))
	 (valstr (concatenate 'string dsname "-VALUE"))
	 (iterlist (loop for k from 1 upto nds collect (gensym "ITER")))
	 (fnstr (princ-to-string fn))
	 (ffnstr (princ-to-string fn-firsti))
	 (nfnstr (princ-to-string fn-nexti))
	 (fncall (cond ((< (length fnstr) (length "#'"))
			`(,fn))
		       ((string= "#'" (subseq fnstr 0 2))
			`(funcall ,fn))
		       (t
			`(,fn))))
	 (ffncall (cond ((< (length ffnstr) (length "#'"))
			 `(,fn-firsti))
			((string= "#'" (subseq ffnstr 0 2))
			 `(funcall ,fn-firsti))
			(t
			 `(,fn-firsti))))
	 (nfncall (cond ((< (length nfnstr) (length "#'"))
			 `(,fn-nexti))
			((string= "#'" (subseq nfnstr 0 2))
			 `(funcall ,fn-nexti))
			(t
			 `(,fn-nexti)))))
    `(let (,@(mapcar #'(lambda (iter ds)
			 `(,iter (,@ffncall ,ds)))
		     iterlist
		     dslist))
       (loop 
	  ,@(mapcar #'(lambda (iter)
			`(when (= -1 ,iter)
			   (return)))
		    iterlist)
	  (,@fncall ,@(mapcan #'copy-list
			       (mapcar #'(lambda (iter ds)
					   `((,(intern keystr) ,ds ,iter)
					     (,(intern valstr) ,ds ,iter)))
				       iterlist
				       dslist)))
	  (setf ,@(mapcan #'copy-list
			  (mapcar #'(lambda (iter ds)
				      `(,iter (,@nfncall ,ds ,iter)))
				  iterlist
				  dslist)))))))
		
	  
    

;;;; A few common data structures are henceforth implemented
;;;; The syntax is unkind but I'm told there's no way around it
;;;; The macroexpansions produce self-explanatory code


(defmacro stack (&key (name "") (key-type 'fixnum) (val-type 'fixnum) (init-key-el 0) (init-val-el 0))
  (let ((namestr (if (string= "" name)
		     (concatenate 'string
				  "STACK-"
				  (princ-to-string key-type)
				  "-"
				  (princ-to-string val-type))
		     (string-upcase name)))
	(s (gensym "STACK"))
	(mores (gensym "MORESTACKS"))
	(i  (gensym "I"))
	(fn (gensym "FN"))
	(key (gensym "KEY"))
	(val (gensym "VAL")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app))))
    `(progn
       (defdatastruct ,(intern namestr) ,key-type ,val-type 0 ,init-key-el ,init-val-el)
       (defun ,(fname "" "-PUSH") (,s ,key ,val)
	 (setf (,(fname "" "-HEADER") ,s) 
	       (,(fname "ADD-IN-" "") ,s ,key ,val)))
       (defun ,(fname "" "-PEEK") (,s)
	 (let ((,i (,(fname "" "-HEADER") ,s)))
	   (if (= -1 ,i)
	       (values nil nil)
	       (values (,(fname "" "-KEY") ,s ,i) (,(fname "" "-VALUE") ,s ,i)))))
       (defun ,(fname "" "-POP") (,s)
	 (let ((,i (,(fname "" "-HEADER") ,s)))
	   (if (= -1 ,i)
	       (values nil nil)
	       (progn 
		 (,(fname "REMOVE-IN-" "") ,s ,i)
		 (decf (,(fname "" "-HEADER") ,s))
		 (values (,(fname "" "-KEY") ,s ,i) (,(fname "" "-VALUE") ,s ,i))))))
       (defmacro ,(fname "MAP-" "") (,fn ,s &rest ,mores)
	 `(defmap ,,namestr ,',(fname "" "-HEADER") #'(lambda (ds i) (decf i)) ,,fn ,,s ,,mores))
       t))))



(defmacro single-linked-list (&key (name "") (key-type 'fixnum) (val-type 'fixnum) (key-equal 'eql) (val-equal 'eql) (init-key-el 0) (init-val-el 0))
  (let ((namestr (if (string= "" name)
		     (concatenate 'string
				  "SINGLE-LINKED-LIST-"
				  (princ-to-string key-type)
				  "-"
				  (princ-to-string val-type))
		     (string-upcase name)))
	(ll (gensym "SLL"))
	(morell (gensym "MORESLL"))
	(i  (gensym "I"))
	(nexti (gensym "I"))
	(key (gensym "KEY"))
	(val (gensym "VAL"))
	(fn (gensym "FN"))
	(ptr (gensym "PTR"))
	(nextptr (gensym "PTR")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app))))
    `(progn
       (defdatastruct ,(intern namestr) ,key-type ,val-type 1 ,init-key-el ,init-val-el)
       (defun ,(fname "" "-PUSH") (,ll ,key ,val)
	 (let ((,ptr (,(fname "" "-HEADER") ,ll)))
	   (setf (,(fname "" "-HEADER") ,ll) 
		 (,(fname "ADD-IN-" "") ,ll ,key ,val ,ptr))))
       (defun ,(fname "" "-POP") (,ll)
	 (let ((,i (,(fname "" "-HEADER") ,ll)))
	   (unless (= -1 ,i)
	     (,(fname "REMOVE-IN-" "") ,ll ,i)
	     (setf (,(fname "" "-HEADER") ,ll) (,(fname "" "-POINTER-1") ,ll ,i)))
	   ,i))
       (defun ,(fname "" "-FIND-KEY") (,ll ,key &optional (,i -1))
	 (when (= -1 ,i)
	   (setf ,i (,(fname "" "-HEADER") ,ll)))
	 (if (= -1 ,i)
	     (values -1 nil)
	     (loop
		(let ((,ptr (,(fname "" "-POINTER-1") ,ll ,i)))
		  (cond ((,key-equal ,key (,(fname "" "-KEY") ,ll ,i))
			 (return (values ,i t)))
			((= -1 ,ptr)
			 (return (values ,i nil)))
			(t
			 (setf ,i ,ptr)))))))
       (defun ,(fname "" "-FIND-VALUE") (,ll ,val &optional (,i -1))
	 (when (= -1 ,i)
	   (setf ,i (,(fname "" "-HEADER") ,ll)))
	 (if (= -1 ,i)
	     (values -1 nil)
	     (loop
		(let ((,ptr (,(fname "" "-POINTER-1") ,ll ,i)))
		  (cond ((,val-equal ,val (,(fname "" "-VALUE") ,ll ,i))
			 (return (values ,i t)))
			((= -1 ,ptr)
			 (return (values ,i nil)))
			(t
			 (setf ,i ,ptr)))))))
       (defun ,(fname "" "-INSERT-AFTER") (,ll ,i ,key ,val)
	 (if (= -1 ,i)
	     (,(fname "" "-PUSH") ,ll ,key ,val)
	     (let* ((,ptr (,(fname "" "-POINTER-1") ,ll ,i))
		    (,nextptr (,(fname "ADD-IN-" "") ,ll ,key ,val ,ptr)))
	       (setf (,(fname "" "-POINTER-1") ,ll ,i) ,nextptr))))
       (defun ,(fname "" "-REMOVE-AFTER") (,ll ,i)
	 (if (= -1 ,i)
	     (,(fname "" "-POP") ,ll)
	     (let ((,nexti (,(fname "" "-POINTER-1") ,ll ,i)))
	       (if (= -1 ,nexti)
		   -1
		   (let ((,nextptr (,(fname "" "-POINTER-1") ,ll ,nexti)))
		     (setf (,(fname "" "-POINTER-1") ,ll ,i) ,nextptr)
		     (,(fname "REMOVE-IN-" "") ,ll ,nexti))))))
       (defmacro ,(fname "MAP-" "") (,fn ,ll &rest ,morell)
	 `(defmap ,,namestr ,',(fname "" "-HEADER") ,',(fname "" "-POINTER-1") ,,fn ,,ll ,,morell))
       t))))



(defmacro binary-tree (&key (name "") (key-type 'fixnum) (val-type 'fixnum) (key-equal 'eql) (key-increasing '<) (init-key-el 0) (init-val-el 0))
  (let ((namestr (if (string= "" name)
		     (concatenate 'string
				  "BINARY-TREE-"
				  (princ-to-string key-type)
				  "-"
				  (princ-to-string val-type))
		     (string-upcase name)))
	(bt (gensym "BT"))
	(i (gensym "I"))
	(ip (gensym "IP"))
	(newi (gensym "NEWI"))
	(key (gensym "KEY"))
	(keyi (gensym "KEY"))
	(val (gensym "VAL"))
	(there (gensym "BOOL"))
	(fn (gensym "FN"))
	(li (gensym "PTR"))
	(ri (gensym "PTR")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app))))
    `(progn
       (defdatastruct ,(intern namestr) ,key-type ,val-type 2 ,init-key-el ,init-val-el)
       (defun ,(fname "" "-FIND-KEY") (,bt ,key &optional (,i -1))
	 (when (= -1 ,i)
	   (setf ,i (,(fname "" "-HEADER") ,bt)))
	 (if (= -1 ,i)
	     (values -1 nil)
	     (loop
		(when (,key-equal ,key (,(fname "" "-KEY") ,bt ,i))
		  (return (values ,i t)))
		(let ((,li (,(fname "" "-POINTER-1") ,bt ,i))
		      (,ri (,(fname "" "-POINTER-2") ,bt ,i))
		      (,keyi (,(fname "" "-KEY") ,bt ,i)))
		  (cond ((and (/= -1 ,li)
			      (,key-increasing ,key ,keyi))
			 (setf ,i ,li))
			((and (/= -1 ,ri)
			      (,key-increasing ,keyi ,key))
			 (setf ,i ,ri))
			(t
			 (return (values ,i nil))))))))
       (defun ,(fname "" "-INSERT") (,bt ,key ,val &optional (,i -1))
	 (when (= -1 ,i)
	   (setf ,i (,(fname "" "-HEADER") ,bt)))
	 (if (= -1 ,i)
	     (values 
	      (setf (,(fname "" "-HEADER") ,bt) 
		    (,(fname "ADD-IN-" "") ,bt ,key ,val -1 -1))
	      t)
	     (multiple-value-bind (,ip ,there)
		 (,(fname "" "-FIND-KEY") ,bt ,key ,i)
	       (if ,there
		   (values ,ip nil)
		   (let ((,newi (,(fname "ADD-IN-" "") ,bt ,key ,val -1 -1)))
		     (if (,key-increasing ,key (,(fname "" "-KEY") ,bt ,ip))
			 (setf (,(fname "" "-POINTER-1") ,bt ,ip) ,newi)
			 (setf (,(fname "" "-POINTER-2") ,bt ,ip) ,newi))
		     (values ,newi t))))))
       (defun ,(fname "" "-FIND-SMALLEST") (,bt)
	 (let ((,i (,(fname "" "-HEADER") ,bt)))
	   (if (= -1 ,i)
	       -1
	       (loop
		  (let ((,li (,(fname "" "-POINTER-1") ,bt ,i)))
		    (when (= -1 ,li)
		      (return ,i))
		    (setf ,i ,li))))))
       (defun ,(fname "" "-FIND-LARGEST") (,bt)
	 (let ((,i (,(fname "" "-HEADER") ,bt)))
	   (if (= -1 ,i)
	       -1
	       (loop
		  (let ((,ri (,(fname "" "-POINTER-2") ,bt ,i)))
		    (when (= -1 ,ri)
		      (return ,i))
		    (setf ,i ,ri))))))
       (defmacro ,(fname "MAP-" "") (,fn ,bt &optional (,i -1))
	 (if (= -1 ,i)
	     `(defmaptree ,,namestr ,,fn ,,bt (,',(fname "" "-HEADER") ,,bt))
	     `(defmaptree ,,namestr ,,fn ,,bt ,,i)))
       t))))



(defmacro threaded-binary-tree (&key (name "") (key-type 'fixnum) (val-type 'fixnum) (key-equal 'eql) (key-increasing '<) (init-key-el 0) (init-val-el 0))
  (let ((namestr (if (string= "" name)
		     (concatenate 'string
				  "THREADED-BINARY-TREE-"
				  (princ-to-string key-type)
				  "-"
				  (princ-to-string val-type))
		     (string-upcase name)))
	(bt (gensym "BT"))
	(morebt (gensym "MOREBT"))
	(fn (gensym "FN"))
	(i (gensym "I"))
	(ip (gensym "IP"))
	(newi (gensym "NEWI"))
	(key (gensym "KEY"))
	(val (gensym "VAL"))
	(there (gensym "BOOL"))
	(previ (gensym "PTR"))
	(nexti (gensym "PTR")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app))))
    `(progn
       (binary-tree :name ,namestr :key-type ,key-type :val-type ,val-type :key-equal ,key-equal :key-increasing ,key-increasing)
       (defdatastruct ,(intern namestr) ,key-type ,val-type 4 ,init-key-el ,init-val-el)
       (defun ,(fname "" "-INSERT") (,bt ,key ,val &optional (,i -1))
	 (when (= -1 ,i)
	   (setf ,i (,(fname "" "-HEADER") ,bt)))
	 (if (= -1 ,i)
	     (values 
	      (setf (,(fname "" "-HEADER") ,bt) 
		    (,(fname "ADD-IN-" "") ,bt ,key ,val -1 -1 -1 -1))
	      t)
	     (multiple-value-bind (,ip ,there)
		 (,(fname "" "-FIND-KEY") ,bt ,key ,i)
	       (cond (,there
		      (values ,ip nil))
		     ((,key-increasing ,key (,(fname "" "-KEY") ,bt ,ip))
		      (let* ((,previ (,(fname "" "-POINTER-3") ,bt ,ip))
			     (,newi (,(fname "ADD-IN-" "") ,bt ,key ,val
				      -1 -1 ,previ ,ip)))
			(setf (,(fname "" "-POINTER-1") ,bt ,ip) ,newi
			      (,(fname "" "-POINTER-3") ,bt ,ip) ,newi)
			(unless (= -1 ,previ)
			  (setf (,(fname "" "-POINTER-4") ,bt ,previ) ,newi))
			(values ,newi t)))
		     (t
		      (let* ((,nexti (,(fname "" "-POINTER-4") ,bt ,ip))
			     (,newi (,(fname "ADD-IN-" "") ,bt ,key ,val
				      -1 -1 ,ip ,nexti)))
			 (setf (,(fname "" "-POINTER-2") ,bt ,ip) ,newi
			       (,(fname "" "-POINTER-4") ,bt ,ip) ,newi)
			 (unless (= -1 ,nexti)
			   (setf (,(fname "" "-POINTER-3") ,bt ,nexti) ,newi))
			 (values ,newi t)))))))
       (defmacro ,(fname "MAP-FORWARD-" "") (,fn ,bt &rest ,morebt)
	 `(defmap ,,namestr ,',(fname "" "-FIND-SMALLEST") ,',(fname "" "-POINTER-4") ,,fn ,,bt ,,morebt))
       (defmacro ,(fname "MAP-BACKWARD-" "") (,fn ,bt &rest ,morebt)
	 `(defmap ,,namestr ,',(fname "" "-FIND-LARGEST") ,',(fname "" "-POINTER-3") ,,fn ,,bt ,,morebt))
       t))))



(defmacro splay-tree (&key (name "") (key-type 'fixnum) (val-type 'fixnum) (key-equal 'eql) (key-increasing '<) (init-key-el 0) (init-val-el 0))
  (let ((namestr (if (string= "" name)
		     (concatenate 'string
				  "SPLAY-TREE-"
				  (princ-to-string key-type)
				  "-"
				  (princ-to-string val-type))
		     (string-upcase name)))
	(st (gensym "ST"))
	(ip (gensym "IP"))
	(newi (gensym "NEWI"))
	(there (gensym "BOOL"))
	(key (gensym "KEY"))
	(val (gensym "VAL"))
	(keyi (gensym "KEY"))
	(keyl (gensym "KEY"))
	(keyr (gensym "KEY"))
	(li (gensym "PTR"))
	(ri (gensym "PTR"))
	(lli (gensym "PTR"))
	(lri (gensym "PTR"))
	(lrli (gensym "PTR"))
	(rli (gensym "PTR"))
	(rri (gensym "PTR")))
    (flet ((fname (pre app)
	     (intern (concatenate 'string pre namestr app))))
    `(progn
       (binary-tree :name ,namestr 
		    :key-type ,key-type 
		    :val-type ,val-type 
		    :key-equal ,key-equal 
		    :key-increasing ,key-increasing 
		    :init-key-el ,init-key-el 
		    :init-val-el ,init-val-el)
       (defun ,(fname "" "-SPLAY") (,st ,key)
	 (if (= -1 (,(fname "" "-HEADER") ,st))
	     (values -1 nil)
	     (loop
		(let* ((,ip (,(fname "" "-HEADER") ,st))
		       (,li (,(fname "" "-POINTER-1") ,st ,ip))
		       (,ri (,(fname "" "-POINTER-2") ,st ,ip))
		       (,keyi (,(fname "" "-KEY") ,st ,ip)))
		  (cond ((,key-equal ,keyi ,key)
			 (return (values ,ip t)))
			((and (/= -1 ,li) (,key-increasing ,key ,keyi))
			 ;; splay left
			 (let ((,lli (,(fname "" "-POINTER-1") ,st ,li))
			       (,lri (,(fname "" "-POINTER-2") ,st ,li))
			       (,keyl (,(fname "" "-KEY") ,st ,li)))
			   (cond ((,key-equal ,key ,keyl)
				  ;; zig
				  (setf (,(fname "" "-POINTER-2") ,st ,li) ,ip
					(,(fname "" "-POINTER-1") ,st ,ip) ,lri
					(,(fname "" "-HEADER") ,st) ,li))
				 ((and (/= -1 ,lli) (,key-increasing ,key ,keyl))
				  ;; zig zig
				  (setf (,(fname "" "-POINTER-1") ,st ,li) (,(fname "" "-POINTER-2") ,st ,lli)
					(,(fname "" "-POINTER-2") ,st ,li) ,ip
					(,(fname "" "-POINTER-1") ,st ,ip) ,lri
					(,(fname "" "-POINTER-2") ,st ,lli) ,li
					(,(fname "" "-HEADER") ,st) ,lli))
				 ((and (/= -1 ,lri) (,key-increasing ,keyl ,key))
				  ;; zig zag
				  (setf (,(fname "" "-POINTER-1") ,st ,ip) (,(fname "" "-POINTER-2") ,st ,lri)
					(,(fname "" "-POINTER-2") ,st ,li) (,(fname "" "-POINTER-1") ,st ,lri)
					(,(fname "" "-POINTER-1") ,st ,lri) ,li
					(,(fname "" "-POINTER-2") ,st ,lri) ,ip
					(,(fname "" "-HEADER") ,st) ,lri))
				 (t
				  (return (values ,ip nil))))))
			((and (/= -1 ,ri) (,key-increasing ,keyi ,key))
			 ;; splay right
			 (let ((,rli (,(fname "" "-POINTER-1") ,st ,ri))
			       (,rri (,(fname "" "-POINTER-2") ,st ,ri))
			       (,keyr (,(fname "" "-KEY") ,st ,ri)))
			   (cond ((,key-equal ,key ,keyr)
				  ;; zig
				  (setf (,(fname "" "-POINTER-1") ,st ,ri) ,ip
					(,(fname "" "-POINTER-2") ,st ,ip) ,rli
					(,(fname "" "-HEADER") ,st) ,ri))
				 ((and (/= -1 ,rli) (,key-increasing ,key ,keyr))
				  ;; zig zag
				  (setf (,(fname "" "-POINTER-2") ,st ,ip) (,(fname "" "-POINTER-1") ,st ,rli)
					(,(fname "" "-POINTER-1") ,st ,ri) (,(fname "" "-POINTER-2") ,st ,rli)
					(,(fname "" "-POINTER-1") ,st ,rli) ,ip
					(,(fname "" "-POINTER-2") ,st ,rli) ,ri
					(,(fname "" "-HEADER") ,st) ,rli))
				 ((and (/= -1 ,rri) (,key-increasing ,keyr ,key))
				  ;; zig zig
				  (setf (,(fname "" "-POINTER-2") ,st ,ri) (,(fname "" "-POINTER-1") ,st ,rri)
					(,(fname "" "-POINTER-1") ,st ,ri) ,ip
					(,(fname "" "-POINTER-2") ,st ,ip) ,rli
					(,(fname "" "-POINTER-1") ,st ,rri) ,ri
					(,(fname "" "-HEADER") ,st) ,rri))
				 (t
				  (return (values ,ip nil))))))
			(t
			 (return (values ,ip nil))))))))
       (defun ,(fname "" "-REMOVE") (,st ,key)
	 (multiple-value-bind (,ip ,there)
	     (,(fname "" "-SPLAY") ,st ,key)
	   (if (not ,there)
	       (values -1 nil)
	       (if (= -1 (,(fname "" "-POINTER-1") ,st ,ip))
		   (let ((,ri (,(fname "" "-POINTER-2") ,st ,ip)))
		     (setf (,(fname "" "-HEADER") ,st) ,ri)
		     (values (,(fname "REMOVE-IN-" "") ,st ,ip) t))
		   (loop
		      (let* ((,li (,(fname "" "-POINTER-1") ,st ,ip))
			     (,lri (,(fname "" "-POINTER-2") ,st ,li)))
			(if (= -1 ,lri)
			    (let ((,ri (,(fname "" "-POINTER-2") ,st ,ip)))
			      (setf (,(fname "" "-HEADER") ,st) ,li
				    (,(fname "" "-POINTER-2") ,st ,li) ,ri)
			      (return (values (,(fname "REMOVE-IN-" "") ,st ,ip) t)))
			    (let ((,lrli (,(fname "" "-POINTER-1") ,st ,lri)))
			    ;; left rotation on left child of root
			    (setf (,(fname "" "-POINTER-1") ,st ,lri) ,li
				  (,(fname "" "-POINTER-1") ,st ,ip) ,lri
				  (,(fname "" "-POINTER-2") ,st ,li) ,lrli)))))))))
       (defun ,(fname "" "-SET") (,st ,key ,val)
	 (if (= -1 (,(fname "" "-HEADER") ,st))
	     (values 
	      (setf (,(fname "" "-HEADER") ,st) 
		    (,(fname "ADD-IN-" "") ,st ,key ,val -1 -1))
	      t)
	     (multiple-value-bind (,ip ,there)
		 (,(fname "" "-SPLAY") ,st ,key)
	       (if ,there
		   (progn
		     (setf (,(fname "" "-VALUE") ,st ,ip) ,val)
		     (values ,ip nil))
		   (let ((,newi (,(fname "ADD-IN-" "") ,st ,key ,val -1 -1))
			 (,li (,(fname "" "-POINTER-1") ,st ,ip))
			 (,ri (,(fname "" "-POINTER-2") ,st ,ip))
			 (,keyi (,(fname "" "-KEY") ,st ,ip)))
		     (cond ((,key-increasing ,key ,keyi)
			    (if (= -1 ,li)
				(setf (,(fname "" "-POINTER-1") ,st ,ip) ,newi)
				(let ((,lli (,(fname "" "-POINTER-1") ,st ,li))
				      (,lri (,(fname "" "-POINTER-2") ,st ,li))
				      (,keyl (,(fname "" "-KEY") ,st ,li)))
				  (cond ((and (= -1 ,lli) (,key-increasing ,key ,keyl))
					 (setf (,(fname "" "-POINTER-1") ,st ,li) ,newi))
					((and (= -1 ,lri) (,key-increasing ,keyl ,key))
					 (setf (,(fname "" "-POINTER-2") ,st ,li) ,newi))
					(t 
					 (error "splay tree error during set"))))))
			   ((,key-increasing ,keyi ,key)
			    (if (= -1 ,ri)
				(setf (,(fname "" "-POINTER-2") ,st ,ip) ,newi)
				(let ((,rli (,(fname "" "-POINTER-1") ,st ,ri))
				      (,rri (,(fname "" "-POINTER-2") ,st ,ri))
				      (,keyr (,(fname "" "-KEY") ,st ,ri)))
				  (cond ((and (= -1 ,rli) (,key-increasing ,key ,keyr))
					 (setf (,(fname "" "-POINTER-1") ,st ,ri) ,newi))
					((and (= -1 ,rri) (,key-increasing ,keyr ,key))
					 (setf (,(fname "" "-POINTER-2") ,st ,ri) ,newi))
					(t 
					 (error "splay tree error during set"))))))
			   (t
			    (error "splay tree error during set")))
		     (values ,newi t))))))
       t))))

 
