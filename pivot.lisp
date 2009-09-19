;;;;; Dynamic Markowitz pivot implementation
;;;;; potential pivots are stored in a min-heap according to #'pivot-test
;;;;;



;;;; Misc pivot heap definitions 
(defun pivot-heap-size (b)
  (length (basis-heap-js b)))

(defun check-pivot-heap-elements (b)
  (let ((n 0))
    (dotimes (index (length (basis-href->hi b)))
      (unless (= -1 (aref (basis-href->hi b) index))
	(incf n)))
    (assert (= n (pivot-heap-size b))))
  (dotimes (hi (pivot-heap-size b))
    (let ((i (aref (basis-heap-is b) hi))
	  (j (aref (basis-heap-js b) hi))
	  (ci (aref (basis-heap-cis b) hi)))
      (assert (= i (aref (lu-eta-matrix-is (aref (basis-l-columns b) j)) ci))))))

(defun check-pivot-heap (b)
  (dotimes (hi (pivot-heap-size b))
    (let ((lc-hi (+ (* 2 hi) 1))
	  (rc-hi (* 2 (+ hi 1)))
	  (j (aref (basis-heap-js b) hi))
	  (i (aref (basis-heap-is b) hi))
	  (ci (aref (basis-heap-cis b) hi)))
      (when (< lc-hi (pivot-heap-size b))
	(let ((lc-j (aref (basis-heap-js b) lc-hi))
	      (lc-i (aref (basis-heap-is b) lc-hi))
	      (lc-ci (aref (basis-heap-cis b) lc-hi)))
	  (assert (pivot-test b j i ci lc-j lc-i lc-ci))))
      (when (< rc-hi (pivot-heap-size b))
	(let ((rc-j (aref (basis-heap-js b) rc-hi))
	      (rc-i (aref (basis-heap-is b) rc-hi))
	      (rc-ci (aref (basis-heap-cis b) rc-hi)))
	  (assert (pivot-test b j i ci rc-j rc-i rc-ci)))))))

(defun print-pivot-heap (b)
  (let ((lhi 0)
	(rhi 0))
    (loop 
       (unless (< lhi (pivot-heap-size b))
	 (return))
       (loop for hi from lhi to rhi
	  do (let ((j (aref (basis-heap-js b) hi))
		   (i (aref (basis-heap-is b) hi))
		   (ci (aref (basis-heap-cis b) hi)))
	       (let* ((cs (* (- (aref (basis-row-nnz b) i) 1) (- (aref (basis-col-nnz b) j) 1)))
		      (l (aref (basis-l-columns b) j))
		      (fs (float (* (lu-eta-matrix-coef l) (aref (lu-eta-matrix-is l) ci)))))
		 (format t "(~A, ~A | ~A, ~3,1F)   " i j cs fs))))
       (format t "~%")
       (setf lhi (+ (* 2 lhi) 1))
       (setf rhi (* 2 (+ rhi 1))))
    (format t "~%")))

(defun check-nz-counts (b)
  (dotimes (k (basis-size b))
    (assert (not (zerop (aref (basis-col-nnz b) k))))
    (assert (not (zerop (aref (basis-row-nnz b) k))))))
    
(defmacro pivot-heap-swap (b hi1 j1 i1 ci1 hi2 j2 i2 ci2)
  (let ((href1 (gensym))
	(href2 (gensym)))
    `(let ((,href1 (aref (basis-hi->href ,b) ,hi1))
	   (,href2 (aref (basis-hi->href ,b) ,hi2)))
       (rotatef (aref (basis-hi->href ,b) ,hi1) (aref (basis-hi->href ,b) ,hi2))
       (rotatef (aref (basis-href->hi ,b) ,href1) (aref (basis-href->hi ,b) ,href2))
       (rotatef ,j1 ,j2)
       (rotatef ,i1 ,i2)
       (rotatef ,ci1 ,ci2)
       (rotatef ,hi1 ,hi2))))



;;;; Pivot criterion
;;;; returns t if heap element hi1 is better than heap element hi2
(defun pivot-test (b j1 i1 ci1 j2 i2 ci2)
  (let ((cc1 (aref (basis-col-nnz b) j1))
	(rc1 (aref (basis-row-nnz b) i1))
	(cc2 (aref (basis-col-nnz b) j2))
	(rc2 (aref (basis-row-nnz b) i2))
	(l1 (aref (basis-l-columns b) j1))
	(l2 (aref (basis-l-columns b) j2)))
    (let ((vi1 (aref (lu-eta-matrix-vis l1) ci1))
	  (vi2 (aref (lu-eta-matrix-vis l2) ci2)))
      (cond ((or (zerop cc1) (zerop rc1) (zerop cc2) (zerop rc2)) (error "bad nnzs"))
	    ((zerop vi2) t)
	    ((zerop vi1) nil)
	    ((or (= most-positive-fixnum cc2) (= most-positive-fixnum rc2)) t)
	    ((or (= most-positive-fixnum cc1) (= most-positive-fixnum rc1)) nil)
	    ((= 1 rc1) t)
	    ((= 1 rc2) nil)
	    ((= 1 cc1) t)
	    ((= 1 cc2) nil)
	    ((< (* (- cc1 1) (- rc1 1)) (* (- cc2 1) (- rc2 1))) t)
	    ((> (* (- cc1 1) (- rc1 1)) (* (- cc2 1) (- rc2 1))) nil)
	    ((>= (abs (* (lu-eta-matrix-coef l1) vi1))
		 (abs (* (lu-eta-matrix-coef l2) vi2)))
	     t)))))


;;;; Sift-down in pivot heap
(defun pivot-heap-sift-down (b hi himax)
  (let ((heap-size (pivot-heap-size b))
	(heap-js (basis-heap-js b))
	(heap-is (basis-heap-is b))
	(heap-cis (basis-heap-cis b)))
    (loop
       (assert (/= -1 hi))
       (when (> hi himax)
	 (return))
       (let ((lc-hi (+ (ash hi 1) 1))
	     (rc-hi (ash (+ hi 1) 1)))
	 (assert (/= -1 lc-hi))
	 (assert (/= -1 rc-hi))
	 (assert (/= -1 (aref (basis-hi->href b) lc-hi)))
	 (assert (< lc-hi heap-size))
	 (assert (<= rc-hi heap-size))
	 (unless (= rc-hi heap-size)
	   (assert (/= -1 (aref (basis-hi->href b) rc-hi))))
	 (symbol-macrolet 
	     ((j (aref heap-js hi))
	      (i (aref heap-is hi))
	      (ci (aref heap-cis hi))
	      (lc-j (aref heap-js lc-hi))
	      (lc-i (aref heap-is lc-hi))
	      (lc-ci (aref heap-cis lc-hi))
	      (rc-j (aref heap-js rc-hi))
	      (rc-i (aref heap-is rc-hi))
	      (rc-ci (aref heap-cis rc-hi)))
	   (if (or (= rc-hi heap-size)
		   (pivot-test b lc-j lc-i lc-ci rc-j rc-i rc-ci))
	       (if (pivot-test b j i ci lc-j lc-i lc-ci)
		   (return)
		   (pivot-heap-swap b hi j i ci lc-hi lc-j lc-i lc-ci))
	       (if (pivot-test b j i ci rc-j rc-i rc-ci)
		   (return)
		   (pivot-heap-swap b hi j i ci rc-hi rc-j rc-i rc-ci))))))))



;;;; Sift-up in pivot heap
(defun pivot-heap-sift-up (b hi)
  (let ((heap-js (basis-heap-js b))
	(heap-is (basis-heap-is b))
	(heap-cis (basis-heap-cis b)))
    (loop
       (when (zerop hi)
	 (return))
       (let ((p-hi (ash (- hi 1) -1)))
	 (symbol-macrolet 
	     ((j (aref heap-js hi))
	      (i (aref heap-is hi))
	      (ci (aref heap-cis hi))
	      (p-j (aref heap-js p-hi))
	      (p-i (aref heap-is p-hi))
	      (p-ci (aref heap-cis p-hi)))
	   (if (pivot-test b p-j p-i p-ci j i ci)
	       (return)
	       (pivot-heap-swap b hi j i ci p-hi p-j p-i p-ci)))))))



;;;; Add new nonzero at bottom of pivot heap
(defun pivot-heap-append (b j i ci)
  (let ((hi (pivot-heap-size b)))
    (vector-push-extend j (basis-heap-js b))
    (vector-push-extend i (basis-heap-is b))
    (vector-push-extend ci (basis-heap-cis b))
    (let ((href (length (basis-href->hi b))))
      (vector-push-extend href (basis-hi->href b))
      (vector-push-extend hi (basis-href->hi b))
      href)))



;;;; Remove a pivot off the heap
(defun pivot-heap-remove (b href)
  (let ((hi (aref (basis-href->hi b) href)))
    (unless (= hi -1)
      (assert (not (zerop (pivot-heap-size b))))
      (let ((last-href (vector-pop (basis-hi->href b)))
	    (last-j (vector-pop (basis-heap-js b)))
	    (last-i (vector-pop (basis-heap-is b)))
	    (last-ci (vector-pop (basis-heap-cis b))))
	(setf (aref (basis-href->hi b) href) -1
	      (aref (basis-hi->href b) hi) -1)
	(unless (= href last-href)
	  (setf (aref (basis-href->hi b) last-href) hi
		(aref (basis-hi->href b) hi) last-href
		(aref (basis-heap-js b) hi) last-j
		(aref (basis-heap-is b) hi) last-i
		(aref (basis-heap-cis b) hi) last-ci)
	  (let ((himax (ash (- (pivot-heap-size b) 2) -1)))
	    (unless (> hi himax)
	      (pivot-heap-sift-down b hi himax)))))
      (check-pivot-heap-elements b))))
       


;;;; Remove the best pivot off the heap
(defun pivot-heap-pop (b)
  (let ((j (aref (basis-heap-js b) 0))
	(i (aref (basis-heap-is b) 0))
	(ci (aref (basis-heap-cis b) 0)))
    (pivot-heap-remove b (aref (basis-hi->href b) 0))
    (check-pivot-heap-elements b)
    (values j i ci)))



;;;; Update a pivot in the heap
(defun pivot-heap-update (b href)
  (let ((hi (aref (basis-href->hi b) href))
	(himax (ash (- (pivot-heap-size b) 2) -1)))
    (unless (= hi -1)
      (pivot-heap-sift-up b hi)
      (setf hi (aref (basis-href->hi b) href))
      (pivot-heap-sift-down b hi himax))))


	     
;;;; Fill basis according to basis header
(defun fill-basis (b lp basis-header)
  (assert (= (basis-size b) (length basis-header)))
  (reset-basis b)
  (dotimes (j (basis-size b))
    (let ((col-ref (aref basis-header j))
	  (l_j     (aref (basis-l-columns b) j)))
      (let* ((col          (aref (lp-columns lp) col-ref))
	     (col-row-refs (column-row-refs col))
	     (col-vals     (column-values col))
	     (n-nz         (length col-row-refs)))
	(setf (lu-eta-matrix-coef l_j)    (column-coef col)
	      (lu-eta-matrix-col-ref l_j) col-ref)
	(dotimes (k n-nz)
	  (let ((i (aref (lp-active-row-inds lp) (aref col-row-refs k)))
		(ci (length (lu-eta-matrix-is l_j)))
		(hi (pivot-heap-size b))
		(val     (aref col-vals k)))
	    (unless (= -1 i)
	      (vector-push-extend i   (lu-eta-matrix-is  l_j))
	      (vector-push-extend val (lu-eta-matrix-vis l_j))
	      (vector-push-extend 0.0 (lu-eta-matrix-vfs l_j))
	      (incf (aref (basis-col-nnz b) j))
	      (incf (aref (basis-row-nnz b) i))
	      (vector-push-extend hi  (basis-href->hi b))
	      (vector-push-extend hi  (basis-hi->href b))
	      (vector-push-extend j   (basis-heap-js b))
	      (vector-push-extend i   (basis-heap-is b))
	      (vector-push-extend ci  (basis-heap-cis b))
	      (vector-push-extend hi  (aref (basis-col-hrefs b) j))
	      (vector-push-extend i   (aref (basis-col-is b) j))
	      (vector-push-extend hi  (aref (basis-row-hrefs b) i))
	      (vector-push-extend j   (aref (basis-row-js b) i))
	      (vector-push-extend ci  (aref (basis-row-cis b) i))))))))
  (dotimes (k (basis-size b))
    (cond ((zerop (aref (basis-col-nnz b) k))
	   (basis-column-is-redundant b k)
	   (return))
	  ((zerop (aref (basis-row-nnz b) k))
	   (basis-row-is-redundant b k)
	   (return))))
  (unless (basis-is-singular b)
    (let ((himax (ash (- (pivot-heap-size b) 2) -1)))
      (loop for hi from himax downto 0 
	 do (pivot-heap-sift-down b hi himax)))
    (check-pivot-heap-elements b)
    t))



;;;; Update the pivot counts and the pivot heap
(defun pivot-count-update (b pivot-j pivot-i)
  (check-pivot-heap b)
  (format t "pivot  = (~A, ~A |~%" pivot-i pivot-j)
  (print-pivot-heap b)
  (check-pivot-heap-elements b)
  ;; remove nzs in pivot column
  (dotimes (index (length (aref (basis-col-is b) pivot-j)))
    (let ((i    (aref (aref (basis-col-is b) pivot-j) index))
	  (href (aref (aref (basis-col-hrefs b) pivot-j) index)))
      (when (and (zerop (decf (aref (basis-row-nnz b) i)))
		 (/= i pivot-i))
	(basis-row-is-redundant b i)
	(return-from pivot-count-update))
      (pivot-heap-remove b href)
      (check-pivot-heap-elements b)))
  ;; remove nzs in pivot row
  (dotimes (index (length (aref (basis-row-js b) pivot-i)))
    (let ((j    (aref (aref (basis-row-js b) pivot-i) index))
	  (href (aref (aref (basis-row-hrefs b) pivot-i) index)))
      (when (and (zerop (decf (aref (basis-col-nnz b) j)))
		 (/= j pivot-j))
	(basis-column-is-redundant b j)
	(return-from pivot-count-update))
      (pivot-heap-remove b href)
      (check-pivot-heap-elements b)))
  ;; remove row/column counts
  (setf (aref (basis-col-nnz b) pivot-j) most-positive-fixnum
	(aref (basis-row-nnz b) pivot-i) most-positive-fixnum)
  (check-nz-counts b)
  (print-pivot-heap b)
  (check-pivot-heap b)
  ;; update modifications via pivot column
  (dotimes (index (length (aref (basis-col-is b) pivot-j)))
    (let ((i (aref (aref (basis-col-is b) pivot-j) index))
	  (pivot-index -1))
      (assert (< 0 (length (aref (basis-row-js b) i))))
      (unless (= i pivot-i)
	(dotimes (row-index (length (aref (basis-row-js b) i)))
	  (let ((j (aref (aref (basis-row-js b) i) row-index))
		(href (aref (aref (basis-row-hrefs b) i) row-index)))
	    (if (= j pivot-j)
		(setf pivot-index row-index)
		(pivot-heap-update b href))))
	(let ((last-j    (vector-pop (aref (basis-row-js b) i)))
	      (last-href (vector-pop (aref (basis-row-hrefs b) i)))
	      (last-cis  (vector-pop (aref (basis-row-cis b) i))))
	  (unless (= pivot-index (length (aref (basis-row-js b) i)))
	    (setf (aref (aref (basis-row-js b) i) pivot-index)    last-j
		  (aref (aref (basis-row-cis b) i) pivot-index)   last-cis
		  (aref (aref (basis-row-hrefs b) i) pivot-index) last-href))))))
  ;; update modifications via pivot row 
  (dotimes (index (length (aref (basis-row-js b) pivot-i)))
    (let ((j (aref (aref (basis-row-js b) pivot-i) index))
	  (pivot-index -1))
      (assert (< 0 (length (aref (basis-col-is b) j))))
      (unless (= j pivot-j)
	(dotimes (col-index (length (aref (basis-col-is b) j)))
	  (let ((i (aref (aref (basis-col-is b) j) col-index))
		(href (aref (aref (basis-col-hrefs b) j) col-index)))
	    (if (= i pivot-i)
		(setf pivot-index col-index)
		(pivot-heap-update b href))))
	(let ((last-i    (vector-pop (aref (basis-col-is b) j)))
	      (last-href (vector-pop (aref (basis-col-hrefs b) j))))
	  (unless (= pivot-index (length (aref (basis-col-is b) j)))
	    (setf (aref (aref (basis-col-is b) j) pivot-index) last-i
		  (aref (aref (basis-col-hrefs b) j) pivot-index) last-href)))))))


;;;; Select the best pivot in the matrix
(defun basis-select-pivot (b)
  (if (zerop (pivot-heap-size b))
      (values -1 -1 -1)
      (multiple-value-bind (j i ci) 
	  (pivot-heap-pop b)
	(let* ((l (aref (basis-l-columns b) j))
	       (vi (aref (lu-eta-matrix-vis l) ci)))
	  (cond ((not (zerop vi))
		 (pivot-count-update b j i)
		 (if (basis-is-singular b)
		     (values -1 -1 -1)
		     (values j i ci)))
		((zerop (pivot-heap-size b))
		 (basis-row-is-redundant b i)
		 (values -1 -1 -1))
		((= 1 (aref (basis-row-nnz b) i))
		 (basis-row-is-redundant b i)
		 (values -1 -1 -1))
		((= 1 (aref (basis-col-nnz b) j))
		 (basis-column-is-redundant b j)
		 (values -1 -1 -1))
		(t 
		 (error "error in heap")))))))
		   
  
#|
		 (dotimes (hi (pivot-heap-size b))
		   (let* ((l (aref (basis-l-columns b) (aref (basis-heap-js b) hi)))
			  (vi (aref (lu-eta-matrix-vis l) (aref (basis-heap-cis b) hi))))
		     (format t "~A " (abs (* (lu-eta-matrix-coef l) vi)))))
		 (format t "~%")
		 (dotimes (hi (pivot-heap-size b))
		   (let* ((j (aref (basis-heap-js b) hi))
			  (i (aref (basis-heap-is b) hi))
			  (rc (aref (basis-row-nnz b) i))
			  (cc (aref (basis-col-nnz b) j)))
		     (format t "~A ~A | " rc cc)))
		 (format t "~%")
|#

;;;; Perform a pivot in the matrix
(defun basis-perform-pivot (b pk)
  (multiple-value-bind (j i ci)
      (basis-select-pivot b)
    (if (= -1 j)
	(values -1 -1 -1)
	(let ((perm-i (aref (basis-i->pi b) i))
	      (perm-j (aref (basis-j->pj b) j))
	      (swap-i (aref (basis-pi->i b) pk))
	      (swap-j (aref (basis-pj->j b) pk)))
	  (setf (aref (basis-i->pi b) i) pk
		(aref (basis-j->pj b) j) pk
		(aref (basis-pi->i b) pk) i
		(aref (basis-pj->j b) pk) j
		(aref (basis-i->pi b) swap-i) perm-i
		(aref (basis-j->pj b) swap-j) perm-j
		(aref (basis-pi->i b) perm-i) swap-i
		(aref (basis-pj->j b) perm-j) swap-j
		(aref (basis-lu-ppivots1 b) pk) j
		(aref (basis-lu-ppivots2 b) pk) swap-j)
	  (values j i ci)))))
