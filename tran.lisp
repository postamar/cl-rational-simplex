
;;;; Solves system E <v-lhs> = <v-rhs>
(defun solve-Ev (v eta-column-index eta-column-vector)
  (let ((si-v   -1) (si-eta -1) 
	(di-v   -1) (di-eta -1) 
	(v_i     0) (eta_i  0))
    ;; get column indices
    (let ((num-eta_column 
	   (find-numerator-from-dense-index eta-column-vector eta-column-index)))
      (multiple-value-bind (num-v_column sparse-index-v_column)
	  (find-numerator-from-dense-index v eta-column-index)
	;; update one numerator and denominator
	(setf (sparse-column-denominator v)
	      (* (sparse-column-denominator v)
		 num-eta_column))
	(setf (aref (sparse-column-numerators v) sparse-index-v_column)
	      (* num-v_column
		 (sparse-column-denominator eta-column-vector)))
	;; update the other numerators
	(iterate-through-sparse-columns
	 (v    eta-column-vector
	  si-v si-eta
	  di-v di-eta
	  v_i  eta_i)
	 (when (/= di-v eta-column-index)
	   (setf (aref (sparse-column-numerators v) si-v)
		 (- (* v_i num-eta_column)
		    (* eta_i num-v_column)))))
	(simplify-sparse-column v)))))



;;;; Solves system E^t <v-lhs> = <v-rhs>
(defun solve-vE (v eta-column-index eta-column-vector)
  (let ((si-v   -1) (si-eta -1) 
	(di-v   -1) (di-eta -1) 
	(v_i     0) (eta_i  0))
    ;; get column indices
    (let ((num-eta_column 
	   (find-numerator-from-dense-index eta-column-vector eta-column-index)))
      (multiple-value-bind (num-v_column sparse-index-v_column)
	  (find-numerator-from-dense-index v eta-column-index)
	;; update denominator
	(setf (sparse-column-denominator v)
	      (* (sparse-column-denominator v)
		 num-eta_column))
	;; update one numerator
	(setf (aref (sparse-column-numerators v) sparse-index-v_column)
	      (* num-v_column
		 (sparse-column-denominator eta-column-vector)))
	(iterate-through-sparse-columns
	 (v    eta-column-vector
	  si-v si-eta
          di-v di-eta
          v_i  eta_i)
	 (when (/= di-v eta-column-index)
	   (decf (aref (sparse-column-numerators v) sparse-index-v_column)
		 (* v_i eta_i))))
	;; update the other numerators
	(iterate-through-sparse-columns
	 (v    eta-column-vector
          si-v si-eta
	  di-v di-eta
	  v_i  eta_i)
	 (when (/= di-v eta-column-index)
	   (setf (aref (sparse-column-numerators v) si-v)
		 (* v_i num-eta_column))))
	(simplify-sparse-column v)))))



;;;; Right-multiplies v with k-th lower eta matrix
(defun v.L_k (v k)
  (let ((L_k (aref *eta-basis-sparse-lower-columns* k))
	(denom-L_k (sparse-column-denominator (aref *eta-basis-sparse-lower-columns* k)))
	(si-v   -1) (si-eta -1) 
	(di-v   -1) (di-eta -1) 
	(v_i     0) (eta_i   0)
	(v_k     0) (si-v_k   -1))
    (iterate-through-sparse-columns
     (v    L_k
      si-v si-eta
      di-v di-eta
      v_i  eta_i)
     (incf v_k
	   (* eta_i v_i))
     (if (= k di-v)
	 (setf si-v_k si-v)
	 (setf (aref (sparse-column-numerators v) si-v)
	       (* v_i denom-L_k))))
    (setf (aref (sparse-column-numerators v) si-v_k) v_k
	  (sparse-column-denominator v) (* (sparse-column-denominator v)
					   denom-L_k))
    (simplify-sparse-column v)))
	  


;;;; Left-multiplies v with k-th lower eta matrix
(defun L_k.v (v k)
  (let ((L_k (aref *eta-basis-sparse-lower-columns* k))
	(denom-L_k (sparse-column-denominator (aref *eta-basis-sparse-lower-columns* k)))
	(v_k (find-numerator-from-dense-index v k))
	(si-v   -1) (si-eta -1) 
	(di-v   -1) (di-eta -1) 
	(v_i     0) (eta_i   0))
    ;; update numerators
    (iterate-through-sparse-columns
     (v    L_k
      si-v si-eta
      di-v di-eta
      v_i  eta_i)
     (setf (aref (sparse-column-numerators v) si-v)
	   (if (= k di-v)
	       (* v_i eta_i)
	       (+ (* v_i denom-L_k)
		  (* v_k eta_i)))))
    ;; update denominator
    (setf (sparse-column-denominator v) (* (sparse-column-denominator v)
					   denom-L_k))
    (simplify-sparse-column v)))
       


;;;; Performs backward transformation
(defun BTRAN (v)
  (let ((eta-file-size (fill-pointer *eta-column-index-file*))
	(m (length *basis-header*)))
    ;; go through eta file
    (dotimes (eta-index-aux eta-file-size)
      (let ((eta-index (- eta-file-size 1 eta-index-aux)))
	(solve-vE v 
		  (aref *eta-column-index-file* eta-index)
		  (aref *eta-column-vector-file* eta-index))))
    ;; go through factored basis
    (dotimes (k m)
      (solve-vE v k (aref *eta-basis-sparse-upper-columns* k)))
    (dotimes (aux-k m)
      (let ((k (- m 1 aux-k)))
	(v.L_k v k)
	(sparse-column-permutate v k (aref *eta-basis-permutation-row-indices* k))))))


	   
;;;; Performs forward transformation
(defun FTRAN (v)
  (let ((m (length *basis-header*)))
    ;; go through factored basis
    (dotimes (k m)
      (sparse-column-permutate v k (aref *eta-basis-permutation-row-indices* k))
      (L_k.v v k))
    (dotimes (aux-k m)
      (let ((k (- m 1 aux-k)))
	(solve-Ev v k (aref *eta-basis-sparse-upper-columns* k))))
    ;; go through eta file
    (dotimes (eta-index (fill-pointer *eta-column-index-file*))
      (let ((eta-column-index (aref *eta-column-index-file* eta-index))
	    (eta-column-vector (aref *eta-column-vector-file* eta-index)))
	(solve-Ev v eta-column-index eta-column-vector)))))
  
  

