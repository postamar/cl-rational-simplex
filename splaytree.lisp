
(defstruct splay-tree-node
  (key (error "splay-tree-node constructor") :type fixnum)
  (val (error "splay-tree-node constructor") :type t)
  (lc  nil :type t)
  (rc  nil :type t))

(defstruct splay-tree
  (root-node nil :type t))


(defun splay-tree-access (splay-tree key)
  (assert (splay-tree-p splay-tree))
  (let ((n (splay-tree-root-node splay-tree)))
    (when n
      (loop
	 (let ((lc (splay-tree-node-lc n))
	       (rc (splay-tree-node-rc n))
	       (k  (splay-tree-node-key n)))
	   (cond ((= key k)
		  ;; stop splaying - success
		  (setf (splay-tree-root-node splay-tree) n)
		  (return t))
		 ((and lc (< key k))
		  ;; splay left
		  (let ((llc (splay-tree-node-lc lc))
			(lrc (splay-tree-node-rc lc))
			(lk  (splay-tree-node-key lc)))
		    (cond ((= key lk)
			   ;; zig
			   (let ((new-n lc))
			     (setf (splay-tree-node-rc new-n) n
				   (splay-tree-node-lc n)     lrc
				   n new-n)))
			  ((and llc (< key lk))
			   ;; zig zig
			   (let ((new-n llc))
			     (setf (splay-tree-node-lc lc)    (splay-tree-node-rc llc)
				   (splay-tree-node-rc lc)    n
				   (splay-tree-node-lc n)     lrc
				   (splay-tree-node-rc new-n) lc
				   n new-n)))
			  ((and lrc (> key lk))
			   ;; zig zag
			   (let ((new-n lrc))
			     (setf (splay-tree-node-lc n)     (splay-tree-node-rc lrc)
				   (splay-tree-node-rc lc)    (splay-tree-node-lc lrc)
				   (splay-tree-node-lc new-n) lc
				   (splay-tree-node-rc new-n) n
				   n new-n)))
			  (t
			   ;; failure
			   (setf (splay-tree-root-node splay-tree) n)
			   (return)))))
		 ((and rc (> key k))
		  ;; splay right
		  (let ((rlc (splay-tree-node-lc rc))
			(rrc (splay-tree-node-rc rc))
			(rk  (splay-tree-node-key rc)))
		    (cond ((= key rk)
			   ;; zig
			   (let ((new-n rc))
			     (setf (splay-tree-node-lc new-n) n
				   (splay-tree-node-rc n)     rlc
				   n new-n)))
			  ((and rlc (< key rk))
			   ;; zig zag
			   (let ((new-n rlc))
			     (setf (splay-tree-node-rc n)     (splay-tree-node-lc rlc)
				   (splay-tree-node-lc rc)    (splay-tree-node-rc rlc)
				   (splay-tree-node-lc new-n) n
				   (splay-tree-node-rc new-n) rc
				   n new-n)))
			  ((and rrc (> key rk))
			   ;; zig zig
			   (let ((new-n rrc))
			     (setf (splay-tree-node-rc rc)    (splay-tree-node-lc rrc)
				   (splay-tree-node-lc rc)    n
				   (splay-tree-node-rc n)     rlc
				   (splay-tree-node-lc new-n) rc
				   n new-n)))
			  (t
			   ;; failure
			   (setf (splay-tree-root-node splay-tree) n)
			   (return)))))
		 (t
		  ;; stop splaying - failure
		  (setf (splay-tree-root-node splay-tree) n)
		  (return))))))))



(defun set-in-splay-tree (splay-tree key value)
  ;; TODO..
  (if (splay-tree-access splay-tree key)
      (setf (splay-tree-node-val (splay-tree-root-node splay-tree)) value)
      (let ((new-node (make-splay-tree-node :key key :val value))
	    (root-node (splay-tree-root-node splay-tree)))
	(if (not root-node)
	    (setf (splay-tree-root-node splay-tree) new-node)
	    (let ((lc (splay-tree-node-lc  root-node))
		  (rc (splay-tree-node-rc  root-node))
		  (k  (splay-tree-node-key root-node)))
	      (cond ((< key k)
		     (if (not lc)
			 (setf (splay-tree-node-lc root-node) new-node)
			 (let ((llc (splay-tree-node-lc  lc))
			       (lrc (splay-tree-node-rc  lc))
			       (lk  (splay-tree-node-key lc)))
			   (cond ((and (not llc) (< key lk))
				  (setf (splay-tree-node-lc lc) new-node))
				 ((and (not lrc) (> key lk))
				  (setf (splay-tree-node-rc lc) new-node))
				 (t 
				  (error "splay tree error during set"))))))
		    ((> key k)
		     (if (not rc)
			 (setf (splay-tree-node-rc root-node) new-node)
			 (let ((rlc (splay-tree-node-lc  rc))
			       (rrc (splay-tree-node-rc  rc))
			       (rk  (splay-tree-node-key rc)))
			   (cond ((and (not rlc) (< key rk))
				  (setf (splay-tree-node-lc rc) new-node))
				 ((and (not rrc) (> key rk))
				  (setf (splay-tree-node-rc rc) new-node))
				 (t 
				  (error "splay tree error during set"))))))
		    (t 
		     (error "splay tree error during set")))))
	value)))
  

(defsetf splay-tree-access set-in-splay-tree)
  
  
  
  

