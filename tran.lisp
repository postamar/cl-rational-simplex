;;;;; BTRAN, FTRAN, etc...


(splay-tree :name hyper-sparse-vector-tree :val-type integer)

(stack :name new-non-zero-stack)


(defstruct (tran
	     (:constructor %make-tran))
  (coef      1   :type rational)
  (indices   #() :type vector)
  (values    #() :type vector)
  (non-zeros nil :type hyper-sparse-vector-tree)
  (new-nzs   nil :type new-non-zero-stack)
  (bm        nil :type basis-matrix)
  (u-file    #() :type bit-vector)
  (l-file    #() :type bit-vector))
  
  
(defun make-tran (bm)
  (let ((m (basis-matrix-size bm)))
    (%make-tran 
     :indices   (make-nvector m -1 fixnum)
     :values    (make-nvector m 0 integer)
     :non-zeros (make-hyper-sparse-vector-tree)
     :new-nzs   (make-new-non-zero-stack)
     :bm        bm
     :u-file    (make-nvector m 0 bit)
     :l-file    (make-nvector m 0 bit))))




;;;;; Hyper-sparse vector functions

(defun is-hsvt-component-non-zero (hsvt index)
  (multiple-value-bind (hsvt-index there)
      (hyper-sparse-vector-tree-find-key hsvt index)
    (declare (ignore hsvt-index))
    there))
	     

(defmacro do-hsv (pivot-fun other-fun scale-fun v-inds eta-inds)
  (let ((k1 (gensym))
	(k2 (gensym))
	(i1 (gensym))
	(i2 (gensym))
	(n1 (gensym))
	(n2 (gensym))
	(pivot-i (gensym)))
    `(let ((,k2 1)
	   (,k1 0)
	   (,n1 (length ,v-inds))
	   (,n2 (length ,eta-inds))
	   (,i1 -1)
	   (,i2 -1)
	   (,pivot-i (aref ,eta-inds 0)))
       (unless (= 0 ,n1)
	 (setf ,i1 (aref ,v-inds 0))
	 (when (or (= 1 ,n2)
		   (progn 
		     (setf ,i2 (aref ,eta-inds 1))
		     (loop
			(cond
			  ((= ,i1 ,pivot-i)
			   (funcall ,pivot-fun ,k1)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref ,v-inds ,k1))))
			  ((= ,i1 ,i2)
			   (funcall ,other-fun ,k1 ,k2)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref ,v-inds ,k1)))
			   (if (= ,n2 (incf ,k2))
			       (return t)
			       (setf ,i2 (aref ,eta-inds ,k2))))
			  ((< ,i1 ,i2)
			 (funcall ,scale-fun ,k1)
			   (if (= ,n1 (incf ,k1))
			       (return nil)
			       (setf ,i1 (aref ,v-inds ,k1))))
			  ((> ,i1 ,i2)
			   (if (= ,n2 (incf ,k2))
			       (return t)
			       (setf ,i2 (aref ,eta-inds ,k2))))))))
	   (loop 
	      (if (= (aref ,v-inds ,k1) ,pivot-i)
		     (funcall ,pivot-fun ,k1)
		     (funcall ,scale-fun ,k1))
	      (when (= ,n1 (incf ,k1)) 
		(return))))))))


(defmacro do-hsv-u (pivot-fun other-fun scale-fun v-inds eta-inds eta-ind-seq i->pi)
  (let ((k1 (gensym))
	(k2 (gensym))
	(i1 (gensym))
	(i2 (gensym))
	(n1 (gensym))
	(n2 (gensym)))
    `(let ((,k2 0)
	   (,k1 0)
	   (,n1 (length ,v-inds))
	   (,n2 (length ,eta-inds))
	   (,i1 -1)
	   (,i2 -1))
       (unless (= 0 ,n1)
	 (setf ,i1 (aref ,v-inds 0))
	 (when (progn 
		 (setf ,i2 (aref ,i->pi (aref ,eta-inds (aref ,eta-ind-seq 0))))
		 (loop
		    (cond
		      ((and (= ,i1 ,i2) (= ,k2 (- ,n2 1)))
		       (funcall ,pivot-fun ,k1 (aref ,eta-ind-seq ,k2))
		       (if (= ,n1 (incf ,k1))
			   (return nil)
			   (setf ,i1 (aref ,v-inds ,k1)))
		       (return t))
		      ((= ,i1 ,i2)
		       (funcall ,other-fun ,k1 (aref ,eta-ind-seq ,k2))
		       (if (= ,n1 (incf ,k1))
			   (return nil)
			   (setf ,i1 (aref ,v-inds ,k1)))
		       (if (= ,n2 (incf ,k2))
			   (return t)
			   (setf ,i2 (aref ,i->pi (aref ,eta-inds (aref ,eta-ind-seq ,k2))))))
		      ((< ,i1 ,i2)
		       (funcall ,scale-fun ,k1)
		       (if (= ,n1 (incf ,k1))
			   (return nil)
			   (setf ,i1 (aref ,v-inds ,k1))))
		      ((> ,i1 ,i2)
		       (if (= ,n2 (incf ,k2))
			   (return t)
			   (setf ,i2 (aref ,i->pi (aref ,eta-inds (aref ,eta-ind-seq ,k2)))))))))
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
  (setf (fill-pointer (tran-indices tr)) 0
	(fill-pointer (tran-values tr)) 0)
  ;; add to result vector
  (map-hyper-sparse-vector-tree
   #'(lambda (v-ind v-val)
       (vector-push v-ind (tran-indices tr))
       (vector-push v-val (tran-values tr)))
   (tran-non-zeros tr)))


;;;;
(defun tran-permutation (tr perm)
  (dotimes (k (length (tran-indices tr)))
    (let ((perm-ind (aref perm (aref (tran-indices tr) k))))
      (setf (aref (tran-indices tr) k) perm-ind))))


;;;;; BTRAN functions 

;;;;
(defun is-btran-u-residual-non-zero (tr eta-col-indices eta-col-indices-seq i->pi)
  (let ((eta-k 0)
	(max-eta-k (length eta-col-indices)))
    (unless (<= max-eta-k 1)
      (map-hyper-sparse-vector-tree
       #'(lambda (v-ind v-val)
	   (declare (ignore v-val))
	   (loop
	      (let ((eta-ind (aref i->pi (aref eta-col-indices (aref eta-col-indices-seq eta-k)))))
		(cond ((= eta-ind v-ind)
		       (return-from is-btran-u-residual-non-zero t))
		      ((> eta-ind v-ind)
		       (return))
		      ((< eta-ind v-ind)
		       (incf eta-k)
		       (when (>= eta-k (- max-eta-k 1))
			 (return-from is-btran-u-residual-non-zero nil)))))))
       (tran-non-zeros tr)))
    nil)) 
  


(defun is-btran-l-residual-non-zero (tr eta-col-indices)
  (let ((eta-k 1)
	(max-eta-k (length eta-col-indices)))
    (unless (<= max-eta-k 1)
      (map-hyper-sparse-vector-tree
       #'(lambda (v-ind v-val)
	   (declare (ignore v-val))
	   (loop
	      (let ((eta-ind (aref eta-col-indices eta-k)))
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
	(if (is-btran-u-residual-non-zero tr (hsv-is u) u-seq i->pi)
	    (progn 
	      (setf (bit (tran-u-file tr) j) 1)
	      (unless (is-hsvt-component-non-zero hsv-nz k)
		(hyper-sparse-vector-tree-set hsv-nz k 0)))
	    (when (is-hsvt-component-non-zero hsv-nz k)
	      (setf (bit (tran-u-file tr) j) 1)))))))



;;;; 
(defun btran-l-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    (loop for k from (- (basis-matrix-n-l-file bm) 1) downto 0
       do (let* ((l (aref (basis-matrix-l-file bm) k))
		 (pivot-i (hsv-j l)))
	    (if (is-hsvt-component-non-zero hsv-nz pivot-i)
		(setf (bit (tran-l-file tr) k) 1)
		(when (is-btran-l-residual-non-zero tr (hsv-is l))
		  (setf (bit (tran-l-file tr) k) 1)
		  (hyper-sparse-vector-tree-set hsv-nz pivot-i 0)))))))



(defun btran-solve-eta (tr eta u-seq)
  (let ((residue 0)
	(new-d-fact (numerator (hsv-coef eta))))
    ;; compute residue
    (do-hsv-u
	#'(lambda (result-pivot-k eta-pivot-k)
	    (declare (ignore result-pivot-k))
	    (mulf new-d-fact (aref (hsv-vis eta) eta-pivot-k)))
	#'(lambda (result-k eta-k)
	    (incf residue (* (aref (hsv-vis eta) eta-k)
			     (aref (tran-values tr) result-k))))
	#'identity
	(tran-indices tr)
	(hsv-is eta)
	u-seq
	(basis-matrix-i->pi (tran-bm tr)))
    ;; update coef
    (divf (tran-coef tr) new-d-fact)
    ;; update result values
    (do-hsv-u
	#'(lambda (result-pivot-k eta-pivot-k)
	    (declare (ignore eta-pivot-k))
	    (setf (aref (tran-values tr) result-pivot-k)
		  (- (* (denominator (hsv-coef eta))
			(aref (tran-values tr) result-pivot-k))
		     (* (numerator (hsv-coef eta)) 
			residue))))
      #'(lambda (result-k eta-k)
	  (declare (ignore eta-k))
	  (mulf (aref (tran-values tr) result-k) new-d-fact))
      #'(lambda (result-k)
	  (mulf (aref (tran-values tr) result-k) new-d-fact))
      (tran-indices tr)
      (hsv-is eta)
      u-seq
      (basis-matrix-i->pi (tran-bm tr)))))
	     


(defun btran-multiply-eta (tr eta)
  (let ((residue 0))
    ;; compute residue
    (do-hsv 
	#'(lambda (pivot-result-k)
	    (incf residue (* (aref (hsv-vis eta) 0)
			     (aref (tran-values tr) pivot-result-k))))
	#'(lambda (result-k eta-k)
	    (incf residue (* (aref (hsv-vis eta) eta-k)
			     (aref (tran-values tr) result-k))))
	#'identity
	(tran-indices tr)
	(hsv-is eta))
    ;; update coef
    (divf (tran-coef tr) (denominator (hsv-coef eta)))
    ;; update result values
    (do-hsv 
	#'(lambda (pivot-result-k)
	    (setf (aref (tran-values tr) pivot-result-k)
		  (* (numerator (hsv-coef eta)) residue)))
	#'(lambda (result-k eta-k)
	    (declare (ignore eta-k))
	    (mulf (aref (tran-values tr) result-k)
		  (denominator (hsv-coef eta))))
	#'(lambda (result-k)
	    (mulf (aref (tran-values tr) result-k)
		  (denominator (hsv-coef eta))))
	(tran-indices tr)
	(hsv-is eta))))




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
(defun btran-l-u-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (loop for k from (- (basis-matrix-n-l-file bm) 1) downto 0
       do (unless (zerop (bit (tran-l-file tr) k))
	    (btran-multiply-eta tr (aref (basis-matrix-l-file bm) k))))))


;;;;
(defun btran-u (tr perm-row perm-col inv-perm-row inv-perm-col coef indices values)
  (declare (ignore perm-row inv-perm-col))
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (length indices))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref perm-col (aref indices k))
				    (aref values k)))
    (setf (tran-coef tr) coef))
  (btran-u-non-zero tr)
  (tran-prepare-fill-in tr)
  (btran-u-fill-in tr)
  (tran-permutation tr inv-perm-row) ;; for sure
  (in-place-sort-keys-increasing (tran-indices tr) (tran-values tr)))


;;;;
(defun btran-l-u (tr coef indices values)
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (length indices))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref indices k)
				    (aref values k)))
    (setf (tran-coef tr) coef))
  (btran-l-u-non-zero tr)
  (tran-prepare-fill-in tr)
  (btran-l-u-fill-in tr))



(defun check-btran-u (tr inv-perm-row inv-perm-col coef indices values)
  (when *checks*
    (let* ((bm (tran-bm tr))
	   (m (basis-matrix-size bm))
	   (du (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (vs (make-nvector m 0 rational))
	   (vr (make-nvector m 0 rational))
	   (vt (make-nvector m 0 rational)))
      ;; make dense u
      (dotimes (k m)
	(setf (aref du (aref inv-perm-row k) (aref inv-perm-col k)) 1))
      (dotimes (k (length (basis-matrix-u-columns bm)))
	(let* ((u (aref (basis-matrix-u-columns bm) k))
	       (uj k))
	  (dotimes (l (length (hsv-is u)))
	    (setf (aref du (aref (hsv-is u) l) uj)
		(* (hsv-coef u) (aref (hsv-vis u) l))))))
      ;; make dense vectors
      (dotimes (k (length indices))
	(let ((i (aref indices k)))
	  (setf (aref vs i) (* coef (aref values k)))))
      (dotimes (k (length (tran-indices tr)))
      (setf (aref vr (aref (tran-indices tr) k))
	    (* (tran-coef tr) (aref (tran-values tr) k))))
      (dotimes (j m)
      (dotimes (i m)
	(incf (aref vt j)
	      (* (aref vr i)
		 (aref du i j)))))
      (unless (dotimes (i m t)
		(unless (= (aref vt i) (aref vs i))
		  (return nil)))
      (print vr)
      (print du)
      (print vt)
      (print vs)
      (print '---))
      (dotimes (i m t)
	(assert (= (aref vt i) (aref vs i)))))))
  
  
;;;;
(defun btran (tr perm-row perm-col inv-perm-row inv-perm-col coef indices values)
  (btran-u tr perm-row perm-col inv-perm-row inv-perm-col coef indices values)
  (check-btran-u tr inv-perm-row inv-perm-col coef indices values)
  (btran-l-u tr (tran-coef tr) (tran-indices tr) (tran-values tr)))





;;;;; FTRAN functions	  

;;;;
(defun find-ftran-u-non-zeros (tr eta-col-indices eta-col-indices-seq)
  (let ((hsv-nz (tran-non-zeros tr))
	(i->pi (basis-matrix-i->pi (tran-bm tr)))
	(eta-k 0)
	(max-eta-k (length eta-col-indices)))
    (unless (<= max-eta-k 1)
      (block find-ftran-u-non-zeros-block
	(map-hyper-sparse-vector-tree
	 #'(lambda (v-ind v-val)
	     (declare (ignore v-val))
	     (loop
		(when (>= eta-k max-eta-k)
		  (return-from find-ftran-u-non-zeros-block))
		(let ((eta-ind (aref i->pi (aref eta-col-indices (aref eta-col-indices-seq eta-k)))))
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
	 do (let ((eta-ind (aref i->pi (aref eta-col-indices (aref eta-col-indices-seq eta-last-k)))))
	      (new-non-zero-stack-push (tran-new-nzs tr) eta-ind 0)))
      ;; add new non-zeros
      (loop 
	 (when (= -1 (new-non-zero-stack-header (tran-new-nzs tr)))
	   (return))
	 (multiple-value-bind (ind val)
	     (new-non-zero-stack-pop (tran-new-nzs tr))
	   (hyper-sparse-vector-tree-set hsv-nz ind val))))))



;;;;
(defun find-ftran-l-non-zeros (tr eta-col-indices)
  (let ((hsv-nz (tran-non-zeros tr))
	(eta-k 1)
	(max-eta-k (length eta-col-indices)))
    (unless (<= max-eta-k 1)
      (block find-ftran-l-non-zeros-block
	(map-hyper-sparse-vector-tree
	 #'(lambda (v-ind v-val)
	     (declare (ignore v-val))
	     (loop
		(when (>= eta-k max-eta-k)
		  (return-from find-ftran-l-non-zeros-block))
		(let ((eta-ind (aref eta-col-indices eta-k)))
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
	 do (let ((eta-ind (aref eta-col-indices eta-last-k)))
	      (new-non-zero-stack-push (tran-new-nzs tr) eta-ind 0)))
      ;; add new non-zeros
      (loop 
	 (when (= -1 (new-non-zero-stack-header (tran-new-nzs tr)))
	   (return))
	 (multiple-value-bind (ind val)
	     (new-non-zero-stack-pop (tran-new-nzs tr))
	   (hyper-sparse-vector-tree-set hsv-nz ind val))))))



;;;;
(defun ftran-l-u-non-zero (tr)
  (let* ((bm (tran-bm tr))
	 (hsv-nz (tran-non-zeros tr)))
    ;; go through PL-file
    (dotimes (k (basis-matrix-n-l-file bm))
      (let* ((l (aref (basis-matrix-l-file bm) k))
	     (pivot-i (hsv-j l)))
	(when (is-hsvt-component-non-zero hsv-nz pivot-i)
	  (find-ftran-l-non-zeros tr (hsv-is l))
	  (setf (bit (tran-l-file tr) k) 1))))))



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
	      (find-ftran-u-non-zeros tr (hsv-is u) u-seq)
	      (setf (bit (tran-u-file tr) j) 1))))))


;;;;
(defun ftran-multiply-eta (tr eta)
  (let* ((pivot-i (aref (hsv-is eta) 0))
	 (pivot-k (find-index (tran-indices tr) pivot-i))
	 (pivot-v (aref (tran-values tr) pivot-k)))
    ;; update coef
    (divf (tran-coef tr) (denominator (hsv-coef eta)))
    ;; update result values
    (do-hsv
	#'(lambda (result-pivot-k)
	    (mulf (aref (tran-values tr) result-pivot-k)
		  (* (aref (hsv-vis eta) 0)
		     (numerator (hsv-coef eta)))))
      #'(lambda (result-k eta-k)
	  (setf (aref (tran-values tr) result-k)
		(+ (* (denominator (hsv-coef eta))
		      (aref (tran-values tr) result-k))
		   (* pivot-v
		      (numerator (hsv-coef eta))
		      (aref (hsv-vis eta) eta-k)))))
      #'(lambda (result-k)
	  (mulf (aref (tran-values tr) result-k)
		(denominator (hsv-coef eta))))
      (tran-indices tr)
      (hsv-is eta))))
			    



;;;;
(defun ftran-solve-eta (tr eta u-seq)
  (let* ((i->pi (basis-matrix-i->pi (tran-bm tr)))
	 (eta-pivot-k (aref u-seq (- (length u-seq) 1)))
	 (pivot-i (aref i->pi (aref (hsv-is eta) eta-pivot-k)))
	 (pivot-k (find-index (tran-indices tr) pivot-i))
	 (pivot-v (aref (tran-values tr) pivot-k)))
    ;; update coef
    (divf (tran-coef tr) (numerator (hsv-coef eta)))
    ;; update result values
    (do-hsv-u
	#'(lambda (result-pivot-k eta-pivot-k)
	    (divf (tran-coef tr) (aref (hsv-vis eta) eta-pivot-k))
	    (mulf (aref (tran-values tr) result-pivot-k)
		  (denominator (hsv-coef eta))))
      #'(lambda (result-k eta-k)
	  (setf (aref (tran-values tr) result-k)
		(* (numerator (hsv-coef eta))
		   (- (* (aref (hsv-vis eta) eta-pivot-k)
			 (aref (tran-values tr) result-k))
		      (* pivot-v
			 (aref (hsv-vis eta) eta-k))))))
      #'(lambda (result-k)
	  (mulf (aref (tran-values tr) result-k)
		(* (numerator (hsv-coef eta))
		   (aref (hsv-vis eta) eta-pivot-k))))
      (tran-indices tr)
      (hsv-is eta)
      u-seq
      i->pi)))



;;;;
(defun ftran-l-u-fill-in (tr)
  (let ((bm (tran-bm tr)))
    (dotimes (k (basis-matrix-n-l-file bm))
      do (unless (zerop (bit (tran-l-file tr) k))
	   (ftran-multiply-eta tr (aref (basis-matrix-l-file bm) k))))))


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
(defun ftran-u (tr perm-row perm-col inv-perm-row inv-perm-col coef indices values)
  (declare (ignore perm-col inv-perm-row))
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (length indices))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref perm-row (aref indices k))
				    (aref values k)))
    (setf (tran-coef tr) coef))
  (ftran-u-non-zero tr)
  (tran-prepare-fill-in tr)
  (ftran-u-fill-in tr)
  (tran-permutation tr inv-perm-col)
  (in-place-sort-keys-increasing (tran-indices tr) (tran-values tr)))
  


;;;;
(defun ftran-l-u (tr coef indices values)
  (tran-prepare-non-zero tr)
  (let* ((init-n-nz (length indices))
	 (hsv-nz (tran-non-zeros tr)))
    ;; add to non-zero vector
    (dotimes (k init-n-nz)
      (hyper-sparse-vector-tree-set hsv-nz
				    (aref indices k)
				    (aref values k)))
    (setf (tran-coef tr) coef))
  (ftran-l-u-non-zero tr)
  (tran-prepare-fill-in tr)
  (ftran-l-u-fill-in tr))




(defun check-ftran-u (tr inv-perm-row inv-perm-col coef indices values)
  (when *checks*
  (let* ((bm (tran-bm tr))
	 (m (basis-matrix-size bm))
	 (du (make-array (list m m) :initial-element 0 :element-type 'rational))
	 (vs (make-nvector m 0 rational))
	 (vr (make-nvector m 0 rational))
	 (vt (make-nvector m 0 rational)))
    ;; make dense u
    (dotimes (k m)
      (setf (aref du (aref inv-perm-row k) (aref inv-perm-col k)) 1))
    (dotimes (k (length (basis-matrix-u-columns bm)))
      (let* ((u (aref (basis-matrix-u-columns bm) k))
	     (uj (hsv-j u)))
	(dotimes (l (length (hsv-is u)))
	  (when (zerop l)
	    (assert (= 1 (aref du (aref (hsv-is u) 0) uj))))
	  (setf (aref du (aref (hsv-is u) l) uj)
		(* (hsv-coef u) (aref (hsv-vis u) l))))))
    ;; make dense vectors
    (dotimes (k (length indices))
      (let ((i (aref indices k)))
	(setf (aref vs i) (* coef (aref values k)))))
    (dotimes (k (length (tran-indices tr)))
      (setf (aref vr (aref (tran-indices tr) k))
	    (* (tran-coef tr) (aref (tran-values tr) k))))
    (dotimes (i m)
      (dotimes (k m)
	(incf (aref vt i)
	      (* (aref vr k)
		 (aref du i k)))))
    (print (list (tran-coef tr) (tran-indices tr) (tran-values tr)))
   (print vr)
   (print du)
   (print vt)
   (print vs)
   (print '---)
    (dotimes (i m t)
      (assert (= (aref vt i) (aref vs i)))))))


;;;;
(defun ftran (tr perm-row perm-col inv-perm-row inv-perm-col coef indices values)
  (ftran-l-u tr coef indices values)
  (ftran-u tr perm-row perm-col inv-perm-row inv-perm-col 
	   (tran-coef tr) (tran-indices tr) (tran-values tr)))
