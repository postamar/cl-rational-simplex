(defpackage #:rationalsimplex
  (:use :cl :asdf :sb-profile))

(in-package :rationalsimplex)


(defsystem rationalsimplex
  :name "rational simplex"
  :version "0.9"
  :maintainer "Marius Posta"
  :author "Marius Posta"
  :licence "BSD sans advertising clause (see file COPYING for details)"
  :description "An exact dual-simplex algorithm implementation using rational numbers"
  :serial t
  :components ((:file "misc")
	       (:file "datastruct")
	       (:file "mps")
	       (:file "hsv")
	       (:file "lp"
		      :depends-on ("misc" "datastruct" "mps" "hsv"))
	       (:file "preproc"
		      :depends-on ("lp"))
	       (:file "bmatrix"
		      :depends-on ("lp"))
	       (:file "luupdt"
		      :depends-on ("bmatrix"))
	       (:file "pivot"
		      :depends-on ("bmatrix"))
	       (:file "lufact"
		      :depends-on ("pivot"))
	       (:file "basis"
		      :depends-on ("bmatrix"))
	       (:file "tran"
		      :depends-on ("basis"))
	       (:file "simplex"
		      :depends-on ("basis"))
	       (:file "exitenter"
		      :depends-on ("simplex"))
	       (:file "update" 
		      :depends-on ("simplex"))
	       (:file "iteration" 
		      :depends-on ("simplex" "exitenter" "update"))
	       (:file "main" 
		      :depends-on ("iteration"))))
		      
	       
