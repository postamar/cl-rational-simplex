(defpackage #:ratiosimpl-asd
  (:use :cl :asdf))

(in-package :ratiosimpl-asd)

(defsystem ratiosimpl
  :name "rational simplex"
  :version "wip 0.0.7"
  :components ((:file "common")
	       (:file "datastruct" :depends-on ("common"))
	       (:file "mps" :depends-on ("common"))
	       (:file "lp" :depends-on ("mps"))
	       (:file "preproc" :depends-on ("lp"))
	       (:file "eta" :depends-on ("common"))
	       (:file "bmatrix" :depends-on ("lp" "eta" "datastruct"))
	       (:file "pivot" :depends-on ("bmatrix"))
	       (:file "basis" :depends-on ("bmatrix" "pivot"))
	       (:file "lu" :depends-on ("bmatrix" "pivot"))
	       (:file "debug" :depends-on ("lu"))))
  
  
