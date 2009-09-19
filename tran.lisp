;;;;; BTRAN, FTRAN, etc...


(splay-tree :name hyper-sparse-vector-tree :val-type integer)

(stack :name new-non-zero-stack)


(defstruct (tran
	     (:constructor %make-tran))
  (hsv       (error "tran constructor") :type hsv)
  (hsv-float (error "tran constructor") :type hsv-float)
  (non-zeros (error "tran constructor") :type hyper-sparse-vector-tree)
  (new-nzs   (error "tran constructor") :type new-non-zero-stack)
  (bm        (error "tran constructor") :type basis-matrix)
  (u-file    (error "tran constructor") :type simple-bit-vector)
  (l-file    (error "tran constructor") :type simple-bit-vector))
  
  
(defun make-tran (bm)
  (%make-tran 
   :hsv       (make-hsv)
   :hsv-float (make-hsv-float)
   :non-zeros (make-hyper-sparse-vector-tree)
   :new-nzs   (make-new-non-zero-stack)
   :bm        bm
   :u-file    (make-array (basis-matrix-size bm) 
			  :initial-element 0 :element-type 'bit)
   :l-file    (make-array (length (basis-matrix-l-file bm))
			  :initial-element 0 :element-type 'bit)))




;;;;; Hyper-sparse vector functions

(defun is-hsvt-component-non-zero (hsvt index)
  (multiple-value-bind (hsvt-index there)
      (hyper-sparse-vector-tree-find-key hsvt index)
    (declare (ignore hsvt-index))
    there))
	     

(defmacro do-hsv-l (pivot-fun other-fun scale-fun v eta-l pivot-i)
  (let ((k1 (gensym))
	(k2 (gensym))
	(i1 (gensym))
	(i2 (gensym))
	(n1 (gensym))
	(n2 (gensym)))
    `(let ((,k2 0)
	   (,k1 0)
	   (,n1 (hsv-length ,v))
	   (,n2 (hsv-length ,eta-l))
	   (,i1 -1)
	   (,i2 -1))
       (declare (fixnum ,k2 ,k1 ,n1 ,n2 ,i1 ,i2 ,pivot-i))
       (unless (= 0 ,n1)
	 (setf ,i1 (aref (hsv-is ,v) 0))
	 (when (or (= 0 ,n2)
		   (progn 
		     (setf ,i2 (aref (hsv-is ,eta-l) 0))
		     (loop
			(cond
			  ((= ,i1 ,pivot-i)
			   (funcall ,pivot-fun ,k1)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref (hsv-is ,v) ,k1))))
			  ((= ,i1 ,i2)
			   (funcall ,other-fun ,k1 ,k2)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref (hsv-is ,v) ,k1)))
			   (if (= ,n2 (incf ,k2))
			       (return t)
			       (setf ,i2 (aref (hsv-is ,eta-l) ,k2))))
			  ((< ,i1 ,i2)
			   (funcall ,scale-fun ,k1)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref (hsv-is ,v) ,k1))))
			  ((> ,i1 ,i2)
			   (if (= ,n2 (incf ,k2))
			       (return t)
			       (setf ,i2 (aref (hsv-is ,eta-l) ,k2))))))))
	   (loop 
	      (if (= (aref (hsv-is ,v) ,k1) ,pivot-i)
		     (funcall ,pivot-fun ,k1)
		     (funcall ,scale-fun ,k1))
	      (when (= ,n1 (incf ,k1)) 
		(return))))))))


(defmacro do-hsv-u (pivot-fun other-fun scale-fun v eta-u u-seq i->pi)
  (let ((k1 (gensym))
	(k2 (gensym))
	(i1 (gensym))
	(i2 (gensym))
	(n1 (gensym))
	(n2 (gensym)))
    `(let ((,k2 0)
	   (,k1 0)
	   (,n1 (hsv-length ,v))
	   (,n2 (hsv-length ,eta-u))
	   (,i1 -1)
	   (,i2 (aref ,i->pi (aref (hsv-is ,eta-u) (aref ,u-seq 0)))))
       (declare (fixnum ,k2 ,k1 ,n1 ,n2 ,i1 ,i2))
       (unless (= 0 ,n1)
	 (setf ,i1 (aref (hsv-is ,v) 0))
	 (when (loop 
		  (cond
		    ((zerop (aref (hsv-vis ,eta-u) (aref ,u-seq ,k2)))
		     (if (= ,n2 (incf ,k2))
			 (return t)
			 (setf ,i2 (aref ,i->pi (aref (hsv-is ,eta-u) (aref ,u-seq ,k2))))))
		    ((and (= ,i1 ,i2) (= ,k2 (- ,n2 1)))
		     (funcall ,pivot-fun ,k1 (aref ,u-seq ,k2))
		     (if (= ,n1 (incf ,k1))
			 (return nil)
			 (setf ,i1 (aref (hsv-is ,v) ,k1)))
		     (return t))
		    ((= ,i1 ,i2)
		     (funcall ,other-fun ,k1 (aref ,u-seq ,k2))
		     (if (= ,n1 (incf ,k1))
			 (return nil)
			 (setf ,i1 (aref (hsv-is ,v) ,k1)))
		     (if (= ,n2 (incf ,k2))
			 (return t)
			 (setf ,i2 (aref ,i->pi (aref (hsv-is ,eta-u) (aref ,u-seq ,k2))))))
		    ((< ,i1 ,i2)
		     (funcall ,scale-fun ,k1)
		     (if (= ,n1 (incf ,k1))
			 (return nil)
			 (setf ,i1 (aref (hsv-is ,v) ,k1))))
		    ((> ,i1 ,i2)
		     (if (= ,n2 (incf ,k2))
			 (return t)
			 (setf ,i2 (aref ,i->pi (aref (hsv-is ,eta-u) (aref ,u-seq ,k2))))))))
	   (loop 
	      (funcall ,scale-fun ,k1)
	      (when (= ,n1 (incf ,k1)) 
		(return))))))))
	   

;;;;; TRAN functions

;;;;
(defun tran-prepare-non-zero (tr)
    ;; reset non-zero vector and files 
    (bit-xor (tran-u-file tr) (tran-u-file tr) t)
    (bit-xor (tran-l-file tr) (tran-l-file tr) t)
    (reset-hyper-sparse-vector-tree (tran-non-zeros tr))
    (reset-new-non-zero-stack (tran-new-nzs tr)))


;;;;
(defun tran-prepare-fill-in (tr)
  ;; reset result vector
  (setf (hsv-length (tran-hsv tr)) 0)
  ;; add to result vector
  (map-hyper-sparse-vector-tree
   #'(lambda (v-ind v-val)
       (hsv-add v-ind v-val (tran-hsv tr)))
   (tran-non-zeros tr)))


;;;;
(defun tran-permutation (tr perm)
  (dotimes (k (hsv-length (tran-hsv tr)))
    (let ((perm-ind (aref perm (aref (hsv-is (tran-hsv tr)) k))))
      (setf (aref (hsv-is (tran-hsv tr)) k) perm-ind))))


;;;;; BTRAN functions 

;;;;
(defun is-btran-u-residual-non-zero (tr eta-u u-seq i->pi)
  (let ((eta-k 0)
	(max-eta-k (hsv-length eta-u)))
    (unless (<= max-eta-k 1)
      (map-hyper-sparse-vector-tree
       #'(lambda (v-ind v-val)
	   (declare (ignore v-val))
	   (loop
	      (cond
		((>= eta-k (- max-eta-k 1))
		 (return-from is-btran-u-residual-non-zero nil))
		((zerop (aref (hsv-vis eta-u) (aref u-seq eta-k)))
		 (incf eta-k))
		(t
		 (let ((eta-ind (aref i->pi (aref (hsv-is eta-u) (aref u-seq eta-k)))))
		   (cond 
		     ((= eta-ind v-ind)
		      (return-from is-btran-u-residual-non-zero t))
		     ((> eta-ind v-ind)
		      (return))
		     ((< eta-ind v-ind)
		      (incf eta-k))))))))
       (tran-non-zeros tr)))
    nil)) 
  


(defun is-btran-l-residual-non-zero (tr eta-l)
  (let ((eta-k 0)
	(max-eta-k (hsv-length eta-l)))
    (unless (zerop max-eta-k)
      (map-hyper-sparse-vector-tree
       #'(lambda (v-ind v-val)
	   (declare (ignore v-val))
	   (loop
	      (let ((eta-ind (aref (hsv-is eta-l) eta-k)))
		(cond ((= eta-ind v-ind)
		       (return-from is-btran-l-residual-non-zero t))
		      ((> eta-ind v-ind)
		       (return))
		      ((< eta-ind v-ind)
		       (incf eta-k)
		       (when (>= eta-k max-eta-k)
			 (return-from is-btran-l-residual-non-zero nil)))))))
       (tran-non-zeros tr)))
    nil))


;;;;
(defun btran-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (m (basis-matrix-size bm))
	 (i->pi (basis-matrix-i->pi bm))
	 (pj->j (basis-matrix-pj->j bm))
	 (hsv-nz (tran-non-zeros tr)))
    (dotimes (k m)
      (let* ((j (aref pj->j k))
	     (u (aref (basis-matrix-u-columns bm) j))
	     (u-seq (aref (basis-matrix-u-seqs bm) j)))
	(if (is-btran-u-residual-non-zero tr u u-seq i->pi)
	    (progn 
	      (setf (sbit (tran-u-file tr) j) 1)
	      (unless (is-hsvt-component-non-zero hsv-nz k)
		(hyper-sparse-vector-tree-set hsv-nz k 0)))
	    (when (is-hsvt-component-non-zero hsv-nz k)
	      (setf (sbit (tran-u-file tr) j) 1)))))))



;;;; 
(defun btran-l-f-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    (loop for k from (- (basis-matrix-n-l-factor-file bm) 1) downto 0
       do (let* ((l (aref (basis-matrix-l-file bm) k))
		 (pivot-i (aref (basis-matrix-l-pivot-file bm) k)))
	    (when (is-btran-l-residual-non-zero tr l)
	      (setf (sbit (tran-l-file tr) k) 1)
	      (unless (is-hsvt-component-non-zero hsv-nz pivot-i)
		(hyper-sparse-vector-tree-set hsv-nz pivot-i 0)))))))


;;;;
(defun btran-l-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    (loop for k from (- (basis-matrix-n-l-file bm) 1) 
       downto (basis-matrix-n-l-factor-file bm)
       do (let* ((l (aref (basis-matrix-l-file bm) k))
		 (pivot-i (aref (basis-matrix-l-pivot-file bm) k)))
	    (when (is-hsvt-component-non-zero hsv-nz pivot-i)
	      (find-ftran-l-non-zeros tr l)
	      (setf (bit (tran-l-file tr) k) 1))))))



;;;;
(defun btran-solve-eta (tr eta u-seq)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (declare ((simple-array fixnum 1) u-seq))
  (let ((residue 0)
	(new-d-fact (numerator (hsv-coef eta)))
	(pivot-ci -1))
    (declare (integer residue new-d-fact))
    ;; compute residue
    (do-hsv-u
	#'(lambda (result-pivot-k eta-pivot-k)
	    (setf pivot-ci result-pivot-k)
	    (mulf new-d-fact (aref (hsv-vis eta) eta-pivot-k)))
	#'(lambda (result-k eta-k)
	    (assert (/= eta-k (aref u-seq (- (hsv-length eta) 1))))
	    (incf residue (* (aref (hsv-vis eta) eta-k)
			     (aref (hsv-vis (tran-hsv tr)) result-k))))
	#'(lambda (result-k)
	    (declare (ignore result-k))
	    t)
	(tran-hsv tr)
	eta
	u-seq
	(basis-matrix-i->pi (tran-bm tr)))
    ;; get gcd
    (let* ((pivot-vi (- (* (denominator (hsv-coef eta))
			   (aref (hsv-vis (tran-hsv tr)) pivot-ci))
			(* (numerator (hsv-coef eta)) 
			   residue)))
	   (factor (gcd pivot-vi new-d-fact)))
      ;; update other values coef
      (divf new-d-fact factor)
      ;; update coef
      (divf (hsv-coef (tran-hsv tr)) new-d-fact)
      ;; update result values
      (dotimes (ci (hsv-length (tran-hsv tr)))
	(if (= ci pivot-ci)
	    (setf (aref (hsv-vis (tran-hsv tr)) ci) 
		  (/ pivot-vi factor))
	    (mulf (aref (hsv-vis (tran-hsv tr)) ci) 
		  new-d-fact))))))
	     

;;;;
(defun btran-multiply-eta (tr eta pivot-i)
  (declare (optimize (speed 1) (safety 0) (debug 0)))
  (declare (fixnum pivot-i))
  (let ((residue 0)
	(pivot-ci -1))
    (declare (integer residue))
    ;; compute residue
    (do-hsv-l
	#'(lambda (pivot-result-k)
	    (setf pivot-ci pivot-result-k))
	#'(lambda (result-k eta-k)
	    (incf residue (* (aref (hsv-vis eta) eta-k)
			     (aref (hsv-vis (tran-hsv tr)) result-k))))
	#'identity
	(tran-hsv tr)
	eta
	pivot-i)
    (let* ((new-d-fact (denominator (hsv-coef eta)))
	   (factor (gcd new-d-fact residue)))
      (declare (integer new-d-fact))
      ;; update other values coef
      (divf new-d-fact factor)
      ;; update coef
      (divf (hsv-coef (tran-hsv tr)) new-d-fact)
      ;; update residue
      (divf residue factor)
      ;; update result values
      (dotimes (ci (hsv-length (tran-hsv tr)))
	(if (= ci pivot-ci)
	    (setf (aref (hsv-vis (tran-hsv tr)) ci) 
		  (+ (* new-d-fact
			(aref (hsv-vis (tran-hsv tr)) pivot-ci))
		     (* (numerator (hsv-coef eta)) 
			residue)))
	    (mulf (aref (hsv-vis (tran-hsv tr)) ci) 
		  new-d-fact))))))



;;;;
(defun btran-u-fill-in (tr)
  (let* ((bm (tran-bm tr))
	 (pj->j (basis-matrix-pj->j bm))
	 (m (basis-matrix-size bm)))
    (dotimes (k m)
      (let ((j (aref pj->j k)))
	(unless (zerop (bit (tran-u-file tr) j))
	  (btran-solve-eta tr 
			   (aref (basis-matrix-u-columns bm) j)
			   (aref (basis-matrix-u-seqs bm) j)))))))


;;;;
(defun btran-l-f-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (loop for k from (- (basis-matrix-n-l-factor-file bm) 1) downto 0
       do (unless (zerop (bit (tran-l-file tr) k))
	    (btran-multiply-eta tr 
				(aref (basis-matrix-l-file bm) k)
				(aref (basis-matrix-l-pivot-file bm) k))))))


;;;;
(defun btran-l-u-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (loop for k from (- (basis-matrix-n-l-file bm) 1) 
       downto (basis-matrix-n-l-factor-file bm)
       do (unless (zerop (bit (tran-l-file tr) k))
	   (ftran-multiply-eta tr 
			       (aref (basis-matrix-l-file bm) k) 
			       (aref (basis-matrix-l-pivot-file bm) k))))))


;;;;
(defun btran-u (tr rhs)
  (let* ((init-n-nz (hsv-length rhs))
	 (hsv-nz (tran-non-zeros tr))
	 (perm-col (basis-matrix-j->pj (tran-bm tr)))
	 (inv-perm-row (basis-matrix-pi->i (tran-bm tr))))
    (tran-prepare-non-zero tr)
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref perm-col (aref (hsv-is rhs) k))
				    (aref (hsv-vis rhs) k)))
    (setf (hsv-coef (tran-hsv tr)) (hsv-coef rhs))
    (btran-u-non-zero tr)
    (tran-prepare-fill-in tr)
    (btran-u-fill-in tr)
    (tran-permutation tr inv-perm-row) 
    (hsv-sort-indices-increasing (tran-hsv tr))))


;;;;
(defun btran-l (tr rhs)
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (hsv-length rhs))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref (hsv-is rhs) k)
				    (aref (hsv-vis rhs) k)))
    (setf (hsv-coef (tran-hsv tr)) (hsv-coef rhs))
    (btran-l-u-non-zero tr)
    (btran-l-f-non-zero tr)
    (tran-prepare-fill-in tr)
    (btran-l-u-fill-in tr)
    (btran-l-f-fill-in tr)))



(defun check-btran-u (tr rhs)
  (when *checks*
    (let* ((bm (tran-bm tr))
	   (m (basis-matrix-size bm))
	   (du (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (pu (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (vs (make-array m :initial-element 0 :element-type 'rational))
	   (vr (make-array m :initial-element 0 :element-type 'rational))
	   (vt (make-array m :initial-element 0 :element-type 'rational))
	   (i->pi (basis-matrix-i->pi bm))
	   (j->pj (basis-matrix-j->pj bm)))
      ;; make dense u
      (dotimes (k m)
	(let ((u (aref (basis-matrix-u-columns bm) k)))
	  (dotimes (l (hsv-length u))
	    (setf (aref pu (aref i->pi k) (aref j->pj (aref (hsv-is u) l)))
		  (* (hsv-coef u) (aref (hsv-vis u) l)))
	    (setf (aref du (aref (hsv-is u) l) k)
		(* (hsv-coef u) (aref (hsv-vis u) l))))))
      ;; make dense vectors
      (dotimes (k (hsv-length rhs))
	(let ((i (aref (hsv-is rhs) k)))
	  (setf (aref vs i) (* (hsv-coef rhs) (aref (hsv-vis rhs) k)))))
      (dotimes (k (hsv-length (tran-hsv tr)))
	(setf (aref vr (aref (hsv-is (tran-hsv tr)) k))
	      (* (hsv-coef (tran-hsv tr)) (aref (hsv-vis (tran-hsv tr)) k))))
      (dotimes (j m)
	(dotimes (i m)
	  (incf (aref vt j)
		(* (aref vr i)
		   (aref du i j)))))
      (unless (dotimes (i m t)
		(unless (= (aref vt i) (aref vs i))
		  (return nil)))
	(print vr)
	(print pu)
	(print vt)
	(print vs)
	(print '---))
      (dotimes (i m t)
	(assert (= (aref vt i) (aref vs i)))))))

  
;;;;
(defun btran (tr rhs)
  (btran-u tr rhs)
  (hsv-remove-zeros (tran-hsv tr))
  (hsv-normalize (tran-hsv tr))
  (check-btran-u tr rhs)
  (btran-l tr (tran-hsv tr))
  (hsv-remove-zeros (tran-hsv tr))
  (hsv-normalize (tran-hsv tr)))
    
    
  

	      


;;;;; FTRAN functions	  

;;;;
(defun find-ftran-u-non-zeros (tr eta-u u-seq)
  (let ((hsv-nz (tran-non-zeros tr))
	(i->pi (basis-matrix-i->pi (tran-bm tr)))
	(eta-k 0)
	(max-eta-k (hsv-length eta-u)))
    (unless (<= max-eta-k 1)
      (block find-ftran-u-non-zeros-block
	(map-hyper-sparse-vector-tree
	 #'(lambda (v-ind v-val)
	     (declare (ignore v-val))
	     (loop
		(cond 
		  ((>= eta-k max-eta-k)
		   (return-from find-ftran-u-non-zeros-block))
		  ((zerop (aref (hsv-vis eta-u) (aref u-seq eta-k)))
		   (incf eta-k))
		  (t
		   (let ((eta-ind (aref i->pi (aref (hsv-is eta-u) (aref u-seq eta-k)))))
		     (cond 
		       ((= eta-ind v-ind)
			(incf eta-k)
			(return))
		       ((> eta-ind v-ind)
			(return))
		       ((< eta-ind v-ind)
			(new-non-zero-stack-push (tran-new-nzs tr) eta-ind 0)
			(incf eta-k))))))))
	 (tran-non-zeros tr)))
      (loop for eta-last-k from eta-k below max-eta-k 
	 do (let ((eta-ci (aref u-seq eta-last-k)))
	      (unless (zerop (aref (hsv-vis eta-u) eta-ci))
		(new-non-zero-stack-push (tran-new-nzs tr) 
					 (aref i->pi (aref (hsv-is eta-u) eta-ci)) 
					 0))))
      ;; add new non-zeros
      (loop 
	 (when (= -1 (new-non-zero-stack-header (tran-new-nzs tr)))
	   (return))
	 (multiple-value-bind (ind val)
	     (new-non-zero-stack-pop (tran-new-nzs tr))
	   (hyper-sparse-vector-tree-set hsv-nz ind val))))))



;;;;
(defun find-ftran-l-non-zeros (tr eta-l)
  (let ((hsv-nz (tran-non-zeros tr))
	(eta-k 0)
	(max-eta-k (hsv-length eta-l)))
    (unless (zerop max-eta-k)
      (block find-ftran-l-non-zeros-block
	(map-hyper-sparse-vector-tree
	 #'(lambda (v-ind v-val)
	     (declare (ignore v-val))
	     (loop
		(when (>= eta-k max-eta-k)
		  (return-from find-ftran-l-non-zeros-block))
		(let ((eta-ind (aref (hsv-is eta-l) eta-k)))
		  (cond ((= eta-ind v-ind)
			 (incf eta-k)
			 (return))
			((> eta-ind v-ind)
			 (return))
			((< eta-ind v-ind)
			 (new-non-zero-stack-push (tran-new-nzs tr) eta-ind 0)
			 (incf eta-k))))))
	 (tran-non-zeros tr)))
      (loop for eta-last-k from eta-k below max-eta-k
	 do (let ((eta-ind (aref (hsv-is eta-l) eta-last-k)))
	      (new-non-zero-stack-push (tran-new-nzs tr) eta-ind 0)))
      ;; add new non-zeros
      (loop 
	 (when (= -1 (new-non-zero-stack-header (tran-new-nzs tr)))
	   (return))
	 (multiple-value-bind (ind val)
	     (new-non-zero-stack-pop (tran-new-nzs tr))
	   (hyper-sparse-vector-tree-set hsv-nz ind val))))))



;;;;
(defun ftran-l-f-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    ;; go through PL-file
    (dotimes (k (basis-matrix-n-l-factor-file bm))
      (let* ((l (aref (basis-matrix-l-file bm) k))
	     (pivot-i (aref (basis-matrix-l-pivot-file bm) k)))
	(when (is-hsvt-component-non-zero hsv-nz pivot-i)
	  (find-ftran-l-non-zeros tr l)
	  (setf (bit (tran-l-file tr) k) 1))))))



;;;;
(defun ftran-l-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    (loop for k from (basis-matrix-n-l-factor-file bm)
       below (basis-matrix-n-l-file bm)
       do (let* ((l (aref (basis-matrix-l-file bm) k))
		 (pivot-i (aref (basis-matrix-l-pivot-file bm) k)))
	    (when (is-btran-l-residual-non-zero tr l)
	      (setf (sbit (tran-l-file tr) k) 1)
	      (unless (is-hsvt-component-non-zero hsv-nz pivot-i)
		  (hyper-sparse-vector-tree-set hsv-nz pivot-i 0)))))))


;;;;
(defun ftran-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (pj->j (basis-matrix-pj->j bm))
	 (m (basis-matrix-size bm))
	 (hsv-nz (tran-non-zeros tr)))
    (loop for k from (- m 1) downto 0
       do (let* ((j (aref pj->j k))
		 (u (aref (basis-matrix-u-columns bm) j))
		 (u-seq (aref (basis-matrix-u-seqs bm) j)))
	    (when (is-hsvt-component-non-zero hsv-nz k)
	      (find-ftran-u-non-zeros tr u u-seq)
	      (setf (bit (tran-u-file tr) j) 1))))))




;;;;
(defun ftran-multiply-eta (tr eta eta-pivot-i)
  (declare (optimize (speed 1) (debug 0) (safety 0)))
  (declare (fixnum eta-pivot-i))
  (let* ((pivot-k (hsv-find eta-pivot-i (tran-hsv tr))) 
	 (new-d-fact (denominator (hsv-coef eta)))
	 (factor (gcd new-d-fact (aref (hsv-vis (tran-hsv tr)) pivot-k)))
	 (pivot-v (* (the integer (numerator (hsv-coef eta)))
		     (the integer (/ (aref (hsv-vis (tran-hsv tr)) pivot-k)
				     factor)))))
    (declare (integer pivot-v new-d-fact))
    ;; update values factor
    (divf new-d-fact factor)
    ;; update coef
    (divf (hsv-coef (tran-hsv tr)) new-d-fact)
    ;; update result values
    (do-hsv-l
	#'(lambda (result-pivot-k)
	    (assert (= result-pivot-k pivot-k))
	    (mulf (aref (hsv-vis (tran-hsv tr)) result-pivot-k)
		  new-d-fact))
      #'(lambda (result-k eta-k)
	  (setf (aref (hsv-vis (tran-hsv tr)) result-k)
		(+ (* new-d-fact
		      (aref (hsv-vis (tran-hsv tr)) result-k))
		   (* pivot-v
		      (aref (hsv-vis eta) eta-k)))))
      #'(lambda (result-k)
	  (mulf (aref (hsv-vis (tran-hsv tr)) result-k)
		new-d-fact))
      (tran-hsv tr)
      eta
      eta-pivot-i)))
      
			    



;;;;
(defun ftran-solve-eta (tr eta-u u-seq)
  (declare (optimize (speed 1) (debug 0) (safety 0)))
  (declare ((simple-array fixnum 1) u-seq))
  (let* ((i->pi (basis-matrix-i->pi (tran-bm tr)))
	 (eta-pivot-k (aref u-seq (- (hsv-length eta-u) 1)))
	 (pivot-i (aref i->pi (aref (hsv-is eta-u) eta-pivot-k)))
	 (pivot-k (hsv-find pivot-i (tran-hsv tr)))
	 (eta-pivot-v (aref (hsv-vis eta-u) eta-pivot-k))
	 (pivot-v (aref (hsv-vis (tran-hsv tr)) pivot-k))
	 (factor (gcd pivot-v eta-pivot-v)))
    (declare (integer eta-pivot-v))
    ;; update pivot values
    (divf eta-pivot-v factor)
    (divf pivot-v factor)
    (let* ((new-n-fact (numerator (hsv-coef eta-u)))
	   (factor2 (gcd pivot-v new-n-fact)))
      (declare (integer new-n-fact))
      ;; update other values factor
      (divf new-n-fact factor2)
      (let ((new-n-fact2 (* new-n-fact eta-pivot-v)))
      ;; update coef
	(divf (hsv-coef (tran-hsv tr)) new-n-fact2)
	;; update result values
	(do-hsv-u
	    #'(lambda (result-pivot-k eta-pivot-k)
		(declare (ignore eta-pivot-k))
		(setf (aref (hsv-vis (tran-hsv tr)) result-pivot-k)
		      (* (the integer (/ pivot-v factor2))
			 (denominator (hsv-coef eta-u)))))
	  #'(lambda (result-k eta-k)
	      (setf (aref (hsv-vis (tran-hsv tr)) result-k)
		    (* new-n-fact
		       (- (* eta-pivot-v
			     (aref (hsv-vis (tran-hsv tr)) result-k))
			  (* pivot-v
			     (aref (hsv-vis eta-u) eta-k))))))
	  #'(lambda (result-k)
	      (mulf (aref (hsv-vis (tran-hsv tr)) result-k)
		    new-n-fact2))
	(tran-hsv tr)
	eta-u
	u-seq
	i->pi)))))



;;;;
(defun ftran-l-u-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (loop for k from (basis-matrix-n-l-factor-file bm)
       below (basis-matrix-n-l-file bm)
       do (unless (zerop (bit (tran-l-file tr) k))
	    (btran-multiply-eta tr 
				(aref (basis-matrix-l-file bm) k)
				(aref (basis-matrix-l-pivot-file bm) k))))))

;;;;
(defun ftran-l-f-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (dotimes (k (basis-matrix-n-l-factor-file bm))
      do (unless (zerop (bit (tran-l-file tr) k))
	   (ftran-multiply-eta tr 
			       (aref (basis-matrix-l-file bm) k)
			       (aref (basis-matrix-l-pivot-file bm) k))))))


;;;;
(defun ftran-u-fill-in (tr)
  (let* ((bm (tran-bm tr))
	 (pj->j (basis-matrix-pj->j bm))
	 (m (basis-matrix-size bm)))
    (loop for k from (- m 1) downto 0
       do (let ((j (aref pj->j k)))
	    (unless (zerop (bit (tran-u-file tr) j))
	      (ftran-solve-eta tr 
			       (aref (basis-matrix-u-columns bm) j)
			       (aref (basis-matrix-u-seqs bm) j)))))))



;;;;
(defun ftran-u (tr rhs)
  (let* ((init-n-nz (hsv-length rhs))
	 (hsv-nz (tran-non-zeros tr))
	 (perm-row (basis-matrix-i->pi (tran-bm tr)))
	 (inv-perm-col (basis-matrix-pj->j (tran-bm tr))))
    (tran-prepare-non-zero tr)
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref perm-row (aref (hsv-is rhs) k))
				    (aref (hsv-vis rhs) k)))
    (setf (hsv-coef (tran-hsv tr)) (hsv-coef rhs))
    (ftran-u-non-zero tr)
    (tran-prepare-fill-in tr)
    (ftran-u-fill-in tr)
    (tran-permutation tr inv-perm-col)
    (hsv-sort-indices-increasing (tran-hsv tr))))
  


;;;;
(defun ftran-l (tr rhs)
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (hsv-length rhs))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref (hsv-is rhs) k)
				    (aref (hsv-vis rhs) k)))
    (setf (hsv-coef (tran-hsv tr)) (hsv-coef rhs)))
  (ftran-l-f-non-zero tr)
  (ftran-l-u-non-zero tr)
  (tran-prepare-fill-in tr)
  (ftran-l-f-fill-in tr)
  (ftran-l-u-fill-in tr))




(defun check-ftran-u (tr rhs)
  (when *checks*
    (let* ((bm (tran-bm tr))
	   (m (basis-matrix-size bm))
	   (du (make-array (list m m) :initial-element 0 :element-type 'rational))
	  ; (pu (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (vs (make-array m :initial-element 0 :element-type 'rational))
	   (vr (make-array m :initial-element 0 :element-type 'rational))
	   (vt (make-array m :initial-element 0 :element-type 'rational)))
      ;; make dense u
      (dotimes (k m)
	(let ((u (aref (basis-matrix-u-columns bm) k)))
	  (dotimes (l (hsv-length u))
	    (setf (aref du (aref (hsv-is u) l) k)
		(* (hsv-coef u) (aref (hsv-vis u) l))))))
      ;; make dense vectors
      (dotimes (k (hsv-length rhs))
	(let ((i (aref (hsv-is rhs) k)))
	  (setf (aref vs i) (* (hsv-coef rhs) (aref (hsv-vis rhs) k)))))
      (dotimes (k (hsv-length (tran-hsv tr)))
	(setf (aref vr (aref (hsv-is (tran-hsv tr)) k))
	      (* (hsv-coef (tran-hsv tr)) (aref (hsv-vis (tran-hsv tr)) k))))
      (dotimes (i m)
	(dotimes (k m)
	  (incf (aref vt i)
		(* (aref vr k)
		   (aref du i k)))))
      (unless (dotimes (i m t)
	      (unless (= (aref vt i) (aref vs i))
		(return)))
	(print vr)
	(print du)
	(print vt)
	(print vs)
	(print '---))
      (dotimes (i m t)
	(assert (= (aref vt i) (aref vs i)))))))


;;;;
(defun ftran (tr v)
  (ftran-l tr v)
  (hsv-remove-zeros (tran-hsv tr))
  (hsv-normalize (tran-hsv tr))
  (if *checks* 
      (let ((rhs (make-hsv)))
	(copy-hsv-into-hsv (tran-hsv tr) rhs)
	(ftran-u tr (tran-hsv tr))
	(check-ftran-u tr rhs))
      (ftran-u tr (tran-hsv tr)))
  (hsv-remove-zeros (tran-hsv tr))
  (hsv-normalize (tran-hsv tr)))



