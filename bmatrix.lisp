(in-package :rationalsimplex)

;;;;; Basis matrix implementation
;;;;;
;;;;; The basis matrix B is a m.m square nonsingular matrix
;;;;; In the course of the simplex algorithm, it is factorized
;;;;; from scratch or the factorization is updated, in order to
;;;;; solve the systems Bx=v and xB=v by forward and backward
;;;;; transformations, respectively
;;;;;
;;;;; Factorizing B from scratch produces a set of eta-matrices 
;;;;; Lf0, Lf1, ..., Lf(nfactor-1) and U0, U1, ..., U(m-1), and
;;;;; permutation matrices P(i->pi), P(pi->i), P(j->pj), P(pj->j) with
;;;;; Lf(n-factor)...Lf1.Lf0.B = U = U0.U1...U(m-1)
;;;;; and P(i->pi).U.P(j->pj) being upper-triangular

;;;;; Changing a column j of B giving a matrix B' is updated
;;;;; by producing an eta matrix Lu and modifying Uj.
;;;;; Hence after k more updates we have, in the L file:
;;;;; Lu(nfile).Lu(nfile-1)...Lu(nfactor).Lf(nfactor-1)...Lf1.Lf
;;;;; with  nfile = nfactor + k + 1 
;;;;;


;;;; Data structure definition
(symbol-macrolet
    ((err (error "basis-matrix constructor")))
  (defstruct (basis-matrix
	       (:constructor %make-basis-matrix))
    (size             0   :type fixnum) ; size of the basis matrix, i.e. m
    (refactorization-period 1 :type fixnum) ; how often the basis is refactorized
    (row-col-max      4   :type fixnum) ; pivot selection criterion (higher->finer)
    (is-singular      nil :type symbol) ; T if basis is singular
    (singular-ref     -1  :type fixnum) ; reference number of row/column causing sing.
    (refs             err :type (simple-array fixnum 1)) ; auxilliary array for fact.
    (flags            err :type (simple-array boolean 1)) ; same
    (l-file           err :type (simple-array hsv 1)) ; contents of the L-file
    ;; indices of pivots in each L-eta matrix
    (l-pivot-file     err :type (simple-array fixnum 1)) 
    (lf-file          err :type (simple-array hsv-float 1)) ; L-file in floats
    (n-l-file         0   :type fixnum) ; total L-file length
    (n-l-factor-file  0   :type fixnum) ; factorization L-file length
    (u-columns        err :type (simple-array hsv 1)) ; U-file
    (uf-columns       err :type (simple-array hsv-float 1)) ; U-file in floats
    ;; sequence in which each element Uj is in order 
    ;; of increasing row indices in the upper triangular matrix
    (u-seqs           err :type (simple-array (simple-array fixnum 1) 1))
    (update-row-vals  err :type (simple-array rational 1)) ; aux. array for update
    (fill-ins         err :type (simple-array fixnum 1)) ; aux. array for fact.
    ;; sparse representation arrays for pivoting during factorization
    (col-is           err :type (simple-array (simple-array fixnum 1) 1))
    (col-nnz          err :type (simple-array fixnum 1))
    (row-js           err :type (simple-array (simple-array fixnum 1) 1))
    (row-cis          err :type (simple-array (simple-array fixnum 1) 1)) 
    (row-nnz          err :type (simple-array fixnum 1))
    ;; buckets of rows/columns by density
    (col-buckets      err :type (simple-array (simple-array fixnum 1) 1))
    (row-buckets      err :type (simple-array (simple-array fixnum 1) 1))
    (col-bucket-sizes err :type (simple-array fixnum 1))
    (row-bucket-sizes err :type (simple-array fixnum 1))
    ;; permutation arrays
    (i->pi            err :type (simple-array fixnum 1))
    (pi->i            err :type (simple-array fixnum 1))
    (j->pj            err :type (simple-array fixnum 1))
    (pj->j            err :type (simple-array fixnum 1))))



;;;; Basis matrix constructor
(defun make-basis-matrix (&key 
			  (m -1) 
			  (lp nil) 
			  (max-l-file-size 2000) 
			  (refactorization-period 100)
			  (row-col-max 4))
  (when lp
    (setf m (adjvector-fixnum-fill-pointer (lp-active-row-refs lp))))
  (incf max-l-file-size m)
  (let* ((hsv0 (make-hsv))
	 (hsvf0 (make-hsv-float))
	 (a0 (make-array 1 :initial-element 0 :element-type 'fixnum))
	 (bm (%make-basis-matrix
	    :size             m
	    :refactorization-period refactorization-period
	    :row-col-max      row-col-max
	    :refs             (make-array m :initial-element -1 :element-type 'fixnum)
	    :flags            (make-array m :initial-element nil :element-type 'boolean)
	    :l-file           (make-array max-l-file-size :initial-element hsv0 :element-type 'hsv)
	    :lf-file          (make-array max-l-file-size :initial-element hsvf0 :element-type 'hsv-float)
	    :l-pivot-file     (make-array max-l-file-size :initial-element -1 :element-type 'fixnum)
	    :u-columns        (make-array m :initial-element hsv0 :element-type 'hsv)
	    :uf-columns       (make-array m :initial-element hsvf0 :element-type 'hsv-float)
	    :u-seqs           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :update-row-vals  (make-array m :initial-element 0 :element-type 'rational)
	    :fill-ins         (make-array m :initial-element 0 :element-type 'fixnum)
	    :col-nnz          (make-array m :initial-element 0 :element-type 'fixnum)
	    :row-nnz          (make-array m :initial-element 0 :element-type 'fixnum)
	    :col-is           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-js           (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-cis          (make-array m :initial-element a0 :element-type '(simple-array fixnum 1))
	    :col-buckets      (make-array (+ m 1) :initial-element a0 :element-type '(simple-array fixnum 1))
	    :row-buckets      (make-array (+ m 1) :initial-element a0 :element-type '(simple-array fixnum 1))
	    :col-bucket-sizes (make-array (+ m 1) :initial-element 0 :element-type 'fixnum)
	    :row-bucket-sizes (make-array (+ m 1) :initial-element 0 :element-type 'fixnum)
	    :j->pj            (make-array m :initial-element -1 :element-type 'fixnum)
	    :i->pi            (make-array m :initial-element -1 :element-type 'fixnum)
	    :pj->j            (make-array m :initial-element -1 :element-type 'fixnum)
	    :pi->i            (make-array m :initial-element -1 :element-type 'fixnum))))
    (dotimes (k max-l-file-size)
      (setf (aref (basis-matrix-l-file bm) k)      (make-hsv)
	    (aref (basis-matrix-lf-file bm) k)     (make-hsv-float)))
    (dotimes (k m bm)
      (setf (aref (basis-matrix-u-columns bm) k)   (make-hsv)
	    (aref (basis-matrix-uf-columns bm) k)  (make-hsv-float)
	    (aref (basis-matrix-u-seqs bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-col-is bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-js bm) k)      (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-cis bm) k)     (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-col-buckets bm) (+ k 1)) (make-array m :initial-element -1 :element-type 'fixnum)
	    (aref (basis-matrix-row-buckets bm) (+ k 1)) (make-array m :initial-element -1 :element-type 'fixnum)))))



;;;; Resets basis prior to factorization
(defun reset-basis-matrix (bm)
  (dotimes (k (basis-matrix-n-l-file bm))
    (reset-hsv-float (aref (basis-matrix-lf-file bm) k))
    (reset-hsv (aref (basis-matrix-l-file bm) k)))
  (setf (basis-matrix-is-singular bm) nil
	(basis-matrix-singular-ref bm) -1
	(basis-matrix-n-l-factor-file bm) 0
	(basis-matrix-n-l-file bm) 0)
  (let ((m (basis-matrix-size bm)))
    (dotimes (k m)
      (setf (aref (basis-matrix-fill-ins bm) k) 0
	    (aref (basis-matrix-refs bm) k) k
	    (aref (basis-matrix-col-nnz bm) k) 0
	    (aref (basis-matrix-row-nnz bm) k) 0
	    (aref (basis-matrix-col-bucket-sizes bm) (+ k 1)) 0
	    (aref (basis-matrix-row-bucket-sizes bm) (+ k 1)) 0
	    (aref (basis-matrix-i->pi bm) k) k
	    (aref (basis-matrix-pi->i bm) k) k
	    (aref (basis-matrix-j->pj bm) k) k
	    (aref (basis-matrix-pj->j bm) k) k)
      (let ((u_k (aref (basis-matrix-u-columns bm) k))
	    (uf_k (aref (basis-matrix-uf-columns bm) k)))
	(reset-hsv u_k)
	(reset-hsv-float uf_k)))))



;;;; Singularity declarations
(defun basis-matrix-row-is-redundant (bm i)
  (setf (basis-matrix-is-singular bm) 'redundant-row
	(basis-matrix-singular-ref bm) i))

(defun basis-matrix-column-is-redundant (bm j)
  (setf (basis-matrix-is-singular bm) 'redundant-column
	(basis-matrix-singular-ref bm) j))

	

;;;; Fill basis according to basis header prior to factorization
;;;; Basis contents are temporarily stored in U
(defun fill-basis-matrix (bm lp basis-header leaving-ref entering-ref)
  (assert (= (basis-matrix-size bm) (length basis-header)))
  (reset-basis-matrix bm)
  (dotimes (j (basis-matrix-size bm))
    (let* ((header-ref (aref basis-header j))
	   (col-ref (if (= header-ref leaving-ref) entering-ref header-ref))
	   (u_j     (aref (basis-matrix-u-columns bm) j))
	   (col          (adjvector-column-ref (lp-columns lp) col-ref))
	   (col-row-refs (hsv-is (column-hsv col)))
	   (col-vals     (hsv-vis (column-hsv col)))
	   (n-nz         (hsv-length (column-hsv col))))
      (setf (hsv-coef u_j) (hsv-coef (column-hsv col)))
      (dotimes (k n-nz)
	(let ((i (adjvector-fixnum-ref (lp-active-row-inds lp)  (aref col-row-refs k)))
	      (ci (hsv-length u_j))
	      (val (aref col-vals k)))
	  (unless (= -1 i)
	    ;; add non-zero component to U
	    (hsv-add i val u_j)
	    ;; add non-zero component to sparse representation
	    (setf (aref (aref (basis-matrix-col-is bm) j) 
			(aref (basis-matrix-col-nnz bm) j)) i)
	    (setf (aref (aref (basis-matrix-row-js bm) i)
			(aref (basis-matrix-row-nnz bm) i)) j)
	    (setf (aref (aref (basis-matrix-row-cis bm) i)
			(aref (basis-matrix-row-nnz bm) i)) ci)
	    (incf (aref (basis-matrix-col-nnz bm) j))
	    (incf (aref (basis-matrix-row-nnz bm) i)))))))
  (dotimes (k (basis-matrix-size bm))
    (let ((nnz-i (aref (basis-matrix-row-nnz bm) k))
	  (nnz-j (aref (basis-matrix-col-nnz bm) k)))
      ;; detect singularities
      (cond ((zerop nnz-j)
	     (basis-matrix-column-is-redundant bm k)
	     (return))
	    ((zerop nnz-i)
	     (basis-matrix-row-is-redundant bm k)
	     (return))
	    (t
	     ;; build buckets
	     (setf (aref (aref (basis-matrix-col-buckets bm) nnz-j)
			 (aref (basis-matrix-col-bucket-sizes bm) nnz-j))
		   k)
	     (incf (aref (basis-matrix-col-bucket-sizes bm) nnz-j))
	     (setf (aref (aref (basis-matrix-row-buckets bm) nnz-i)
			 (aref (basis-matrix-row-bucket-sizes bm) nnz-i))
		   k) 
	     (incf (aref (basis-matrix-row-bucket-sizes bm) nnz-i))))))
  ;; return T on success
  (not (basis-matrix-is-singular bm)))





;;;;; DEBUGGING FUNCTIONS

;;;; Output functions 
(defun print-2d-array (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t "   .   ")
	  (format t "~6,2F " (float (aref a i j)))))))

(defun print-2d-array-nz (a)
  (dotimes (i (array-dimension a 0) (format t "~%"))
    (dotimes (j (array-dimension a 1) (format t "~%"))
      (if (zerop (aref a i j))
	  (format t " .")
	  (format t " x")))))

(defun print-u (bm)
  (let* ((m (basis-matrix-size bm))
	 (ua (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (j m)
      (let* ((u (aref (basis-matrix-u-columns bm) j)))
	(dotimes (r (hsv-length u))
	  (let ((i (aref (hsv-is u) r)))
	    (let ((ip (aref (basis-matrix-i->pi bm) i))
		  (jp (aref (basis-matrix-j->pj bm) j)))
	      (setf (aref ua ip jp)
		    (* (hsv-coef u) (aref (hsv-vis u) r))))))))
    (print-2d-array-nz ua)))

(defun print-l-f (bm)
  (let* ((m (basis-matrix-size bm)))
    (loop for k from (basis-matrix-n-l-factor-file bm) below (basis-matrix-n-l-file bm)
       do (let ((la (make-array (list m m) :initial-element 0 :element-type 'rational))
		(l (aref (basis-matrix-l-file bm) k))
		(li (aref (basis-matrix-l-pivot-file bm) k)))
	    (dotimes (i m)
	      (setf (aref la i i) 1))
	    (dotimes (r (hsv-length l))
	      (setf (aref la li (aref (hsv-is l) r))
		    (* (hsv-coef l) (aref (hsv-vis l) r))))
	    (print-2d-array-nz la)))))



;;;; Builds an m.m array containing the basis
(defun make-dense-basis (lp bm bh)
  (let* ((m (basis-matrix-size bm))
	 (db (make-array (list m m) :initial-element 0 :element-type 'rational)))
    (dotimes (k m db)
      (let* ((col-ref (aref bh k))
	     (col (adjvector-column-ref (lp-columns lp) col-ref)))
	(dotimes (l (hsv-length (column-hsv col)))
	  (let* ((row-ref (aref (hsv-is (column-hsv col)) l))
		 (row-ind (adjvector-fixnum-ref (lp-active-row-inds lp) row-ref)))
	    (unless (= -1 row-ind)
	      (setf (aref db row-ind k)
		    (rational-in-column col l)))))))))

;;;; Checks ordered sequences
(defun check-u-seqs (bm)
  (when *checks* 
    (dotimes (k (basis-matrix-size bm))
      (let* ((j (aref (basis-matrix-pj->j bm) k))
	     (i (aref (basis-matrix-pi->i bm) k))
	     (u (aref (basis-matrix-u-columns bm) j))
	     (u-seq (aref (basis-matrix-u-seqs bm) j))
	     (lastuk (- (hsv-length u) 1)))
	(assert (<= 0 lastuk))
	(assert (= i (aref (hsv-is u) (aref u-seq lastuk))))
	(assert (not (zerop (aref (hsv-vis u) (aref u-seq lastuk)))))
	(dotimes (uk lastuk)
	  (assert (or (zerop (aref (hsv-vis u) (aref u-seq uk)))
		      (< (aref (hsv-is u) uk) (aref (hsv-is u) (+ uk 1))))))))))



;;;; Checks the LU factorization
(defun check-lu (lp bm bh)
  (when *checks*
    (let* ((m (basis-matrix-size bm))
	   (ua (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (la (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (ta (make-array (list m m) :initial-element 0 :element-type 'rational))
	   (oa (make-dense-basis lp bm bh))
	   (da (make-dense-basis lp bm bh)))
      ;; fill u
      (dotimes (j m)
	(let* ((u (aref (basis-matrix-u-columns bm) j)))
	  (dotimes (r (hsv-length u))
	    (let ((i (aref (hsv-is u) r)))
	      (setf (aref ua i j)
		    (* (hsv-coef u) (aref (hsv-vis u) r)))))))
      ;; compute for factorizations
      (dotimes (k (basis-matrix-n-l-factor-file bm))
	;; reset l
	(let ((l (aref (basis-matrix-l-file bm) k))
	      (lj (aref (basis-matrix-l-pivot-file bm) k)))
	  (dotimes (i m)
	    (dotimes (j m)
	      (setf (aref la i j) (if (= i j) 1 0))))
	  (dotimes (kl (hsv-length l))
	    (setf (aref la (aref (hsv-is l) kl) lj)
		  (* (hsv-coef l) (aref (hsv-vis l) kl)))))
	;; matrix multiplication
	(dotimes (i m)
	  (dotimes (j m)
	    (let ((v 0))
	      (dotimes (km m)
		(incf v (* (aref la i km) (aref da km j))))
	      (setf (aref ta i j) v))))
	;; matrix copy
	(dotimes (i m)
	  (dotimes (j m)
	    (setf (aref oa i j) (aref ta i j))
	    (setf (aref da i j) (aref ta i j)))))
      ;; compute for updates
      (loop for k from (basis-matrix-n-l-factor-file bm) below (basis-matrix-n-l-file bm)
	 do (let ((l (aref (basis-matrix-l-file bm) k))
		  (lj (aref (basis-matrix-l-pivot-file bm) k)))
	      ;; reset l
	      (dotimes (i m)
		(dotimes (j m)
		  (setf (aref oa i j) (aref da i j))
		  (setf (aref la i j) (if (= i j) 1 0))))
	      (dotimes (kl (hsv-length l))
		(setf (aref la lj (aref (hsv-is l) kl))
		      (* (hsv-coef l) (aref (hsv-vis l) kl))))
	      ;; matrix multiplication
	      (dotimes (i m)
		(dotimes (j m)
		  (let ((v 0))
		    (dotimes (km m)
		      (incf v (* (aref la i km) (aref da km j))))
		    (setf (aref ta i j) v))))
	      ;; matrix copy
	      (dotimes (i m)
		(dotimes (j m)
		  (setf (aref da i j) (aref ta i j))))))
      ;; check-
      (unless (equalp da ua)
	(format t "~%")
	(dotimes (i m)
	  (dotimes (j m)
	    (let ((ip (aref (basis-matrix-i->pi bm) i))
		  (jp (aref (basis-matrix-j->pj bm) j)))
	      (setf (aref ta ip jp) (aref oa i j))
	      (unless (= (aref ua i j) (aref da i j))
		(format t "discrepancy in (~A ~A), in u (~A ~A), is ~A in LB and is ~A in U~%"
			i j ip jp (aref da i j) (aref ua i j))))))
	(print-2d-array ta)
	(error "bad lu decomposition"))
      ;; check orders
      (dotimes (j m)
	(let ((u (aref (basis-matrix-u-columns bm) j))
	      (u-seq (aref (basis-matrix-u-seqs bm) j)))
	  (dotimes (k (- (hsv-length u) 1))
	    (let ((test 
		   (or (< (aref (basis-matrix-i->pi bm) (aref (hsv-is u) (aref u-seq k)))
			  (aref (basis-matrix-i->pi bm) (aref (hsv-is u) (aref u-seq (+ k 1)))))
		       (zerop (aref (hsv-vis u) (aref u-seq k)))
		       (zerop (aref (hsv-vis u) (aref u-seq (+ k 1)))))))
	      (unless test
		(error "bad row sequences in u")))))))))
