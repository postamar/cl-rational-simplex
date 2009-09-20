(in-package :rationalsimplex)

;;;; Common functions and macros that don't belong anywhere else
;;;;

;;;; 
(defparameter *checks* nil
  "if set to T, activates SLOW result-verifying checks")



;;;; Modify macros
(defun absmax (a b)
  (max (abs a) (abs b)))
(define-modify-macro absmaxf (value) absmax)
(define-modify-macro maxf (value) max)
(define-modify-macro minf (value) min)
(define-modify-macro mulf (value) *)
(define-modify-macro divf (value) /)
(define-modify-macro ashf (value) ash)
(defmacro orf (place &rest rest)
  `(setf ,place (or ,@rest ,place)))


;;;; Performs linear search and replace in vector
;*** Not sure if this is actually used
(defun linear-modify (vector old-value new-value)
  (let ((len (length vector)))
    (dotimes (k len)
      (when (= old-value (aref vector k))
	(setf (aref vector k) new-value)
	(return t)))))



;;;; Binary search in increasing-ordered vector
;;;; returns index on success, -1 on failure
(defun find-index (vector value)
  (declare ((array fixnum *) vector)
	   (fixnum value))
  (let ((l 0)
	(u (length vector))
	(p 0)
	(d 0))
    (declare (fixnum l u p d))
    (loop
       (setf p (ash (- u l) -1))
       (when (zerop p)
	 (if (eql value (aref vector l))
	     (return l)
	     (return -1)))
       (incf p l)
       (setf d (aref vector p))
       (if (eql value d)
	   (return p)
	   (if (< d value)
	       (setf l p)
	       (setf u p))))))


;;;; Binary search in increasing-ordered vector in range 0-nelts
;;;; returns index on success, -1 on failure
(defun find-index-bounded (vector nelts value) 
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (declare ((simple-array fixnum 1) vector)
	   (fixnum nelts value))
  (let ((l 0)
	(u nelts)
	(p 0)
	(d 0))
    (declare (fixnum l u p d))
    (loop
       (setf p (ash (- u l) -1))
       (when (zerop p)
	 (if (eql value (aref vector l))
	     (return l)
	     (return -1)))
       (incf p l)
       (setf d (aref vector p))
       (if (eql value d)
	   (return p)
	   (if (< d value)
	       (setf l p)
	       (setf u p))))))



;;;; Binary search in increasing-ordered vector in range 0-nelts
;;;; returns index on success, index of previous element on failure
(defun find-closest-index-bounded (vector nelts value)
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (declare ((simple-array fixnum 1) vector)
	   (fixnum nelts value))
  (let ((l 0)
	(u nelts)
	(p 0)
	(d 0))
    (declare (fixnum l u p d))
    (cond ((zerop nelts)
	   -1)
	  ((< value (aref vector 0))
	   -1)
	  (t 
	   (loop
	      (setf p (ash (- u l) -1))
	      (when (zerop p)
		(return l))
	      (incf p l)
	      (setf d (aref vector p))
	      (if (eql value d)
		  (return p)
		  (if (< d value)
		      (setf l p)
		      (setf u p))))))))
    
    
;;;; Binary log rounded down of n
(defun floor-log2 (n)
  (let ((p 0))
    (unless (< n 65536)
      (ashf n -16)
      (incf p 16))
    (unless (< n 256)
      (ashf n -8)
      (incf p 8))
    (unless (< n 16)
      (ashf n -4)
      (incf p 4))
    (unless (< n 4)
      (ashf n -2)
      (incf p 2))
    (unless (< n 2)
      (ashf n -1)
      (incf p))
    (if (zerop n)
	-1
	p)))


;;;; Heap sort, increasing order, of vector in range 0-length
(defun sort-increasing-bounded (vector length key)
  (declare ((simple-array fixnum 1) vector)
	   ((function (fixnum) fixnum) key))
  (flet ((sift-down (root end)
	   (loop
	      (let ((child (+ 1 (* 2 root))))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (funcall key (aref vector child)) 
			      (funcall key (aref vector (+ child 1)))))
		  (incf child))
		(when (>= (funcall key (aref vector root)) 
			  (funcall key (aref vector child)))
		  (return))
		(rotatef (aref vector root) (aref vector child))
		(setf root child)))))
    ;; max-heapify
    (loop for heapify-start from (floor (- length 2) 2) downto 0
       do (sift-down heapify-start (- length 1)))
    ;; sort
    (loop for end from (- length 1) downto 1
       do (progn (rotatef (aref vector 0) (aref vector end))
		 (sift-down 0 (- end 1))))))


;;;; In-place heap sort of key-value vector pair, increasing order
(defun in-place-sort-keys-increasing (keys values)
  (flet ((sift-down (root end)
	   (loop
	      (let ((child (+ 1 (* 2 root))))
		(when (> child end)
		  (return))
		(when (and (< child end)
			   (< (aref keys child) (aref keys (+ child 1))))
		  (incf child))
		(when (>= (aref keys root) (aref keys child))
		  (return))
		(rotatef (aref keys root) (aref keys child))
		(rotatef (aref values root) (aref values child))
		(setf root child)))))
    ;; max-heapify
    (let ((count (length keys)))
      (loop for heapify-start from (floor (- count 2) 2) downto 0
	 do (sift-down heapify-start (- count 1)))
      ;; sort
      (loop for end from (- count 1) downto 1
	 do (progn (rotatef (aref keys 0) (aref keys end))
		   (rotatef (aref values 0) (aref values end))
		   (sift-down 0 (- end 1)))))))


;;;; Common denominator of rationals in vector in range 0-count
(defun common-denominator (vector count)
  (assert (not (zerop count)))
  (let ((cd (denominator (aref vector 0))))
    (loop for k from 1 below count
       unless (zerop (aref vector k))
       do (let ((dk (denominator (aref vector k))))
	    (mulf cd (/ dk (gcd cd dk)))))
    cd))
  

      
;;;; Initial random state
(defparameter *simplex-random-state*
  #S(RANDOM-STATE :STATE #.(MAKE-ARRAY 627 :ELEMENT-TYPE '(UNSIGNED-BYTE 32)
				       :INITIAL-CONTENTS
				       '(0 2567483615 130 2392588325 471108515
					 3173724518 1152824458 2132340551
					 2244806021 1313474216 1780589052
					 2314565241 3206115047 2193879442
					 1689714782 4230326843 763475449
					 3184399428 1604139008 130236397
					 2707501403 1882178574 3257938498
					 473441615 3091798989 2644522112
					 2351373380 3555194049 1849441759
					 3436181722 2414457878 2134260483
					 3715007841 2159999292 1735170440
					 2054536437 681931667 2487578070
					 1309884058 2304845079 4075201141
					 2622763192 2869442284 4293411497
					 1927133943 4080136802 1003450158
					 112879339 1883767721 1671134356
					 3660268176 59942333 3292658891
					 3081481086 1447998738 2451097631
					 3911568701 2843963088 4216691060
					 195455153 3821825711 980839722 574180326
					 2143714355 43887121 1485205132 284553560
					 2075371525 2631579075 3622429382
					 17781738 1763421735 4040067941
					 3957866248 4110279836 1031565209
					 4116008391 137498034 4001762942
					 3785169947 3618436825 993369572
					 1799656160 271754189 2866718011
					 2459152174 231392354 3478747567
					 3976374509 2751852576 2286309860
					 3278272289 1053975615 2625422330
					 3281287286 2348451 2349894529 2555187100
					 1609032360 4082289237 2874621875
					 2817778230 1333987322 1969439223
					 1682023125 1712836632 2981373836
					 878542537 4151409495 30547586 2135419214
					 3286096203 2800784521 776288948
					 3316923504 3452033309 789872555
					 2373485109 2680359799 1572694880
					 1802409796 3318639283 1157762585
					 2830249190 3053429998 2073617585
					 2482471571 280086780 947902392
					 1010699839 4219797285 874346594
					 582742962 2080508973 406459247 73512632
					 1904828348 2059990283 872704017 17642350
					 1721597926 2374599305 1374284971
					 2199083060 2543463760 2837723895
					 2622236605 3672148426 2387902506
					 3502725541 1245469319 2432150736
					 2039707092 2038744387 1222957769
					 1678547254 996003390 521502881
					 2966944963 662604332 1083153128
					 4140820719 573918933 3619556082
					 2967582914 4150528733 512761535
					 3329439144 2767801356 2608080347
					 543031745 1642126846 4119147190
					 2909575481 2644472987 3798547236
					 2314329408 1820966119 2397727085
					 576320346 4119854330 3378078485
					 728933143 3462323264 2841745636
					 629644115 2816473017 1889025286
					 2386332686 3679564689 1942315251
					 1458412252 3129570136 1283420191
					 4221938885 1793772738 953702610
					 1066042893 2650034703 266265560
					 3018322012 487461675 1188369009
					 799922574 1686216774 2726019625
					 2602604043 1423338900 1585351664
					 3951724247 1922610013 874369386
					 3190793546 21212165 1023596327
					 3202247984 1565333364 520947555
					 1865457129 519547862 2275726046
					 1290768129 3944961187 1281541644
					 1606036232 2049659087 2093319285
					 559500754 2733373863 4273050268
					 3523281719 3545218462 2254215763
					 2564727560 1519738931 1281933462
					 1768265831 2755061340 3781070367
					 686106678 2995070843 2408858016
					 2219998627 2583419742 4096815095
					 3104747756 1313890679 1144637230
					 90352275 2086355272 1375354627
					 1448749542 2645517783 1076382028
					 2637838223 4220103318 4144020619
					 777347600 946031859 1569742334
					 3988993191 3881453948 1709339159
					 2178810974 3357203827 2436298600
					 3333698579 1185540950 3037764487
					 858778076 3013836575 3316757846
					 3539357755 3737531104 907205123
					 3008051006 947976151 866940012
					 1392943895 618569486 2187427507
					 3743922216 3432124419 480039110
					 3449404631 1951494764 2545996133
					 3830161035 3725784699 1256221411
					 3722816333 749464903 3719085123
					 300853491 2564520405 3700530795
					 181578715 1366764115 3917321541
					 2376551207 830215419 3167290475
					 2344357637 2895458555 2152901691
					 2348698179 1108478861 2493902807
					 2764162115 429175315 2371313349
					 2506940971 313169451 306578275 110923557
					 3799377607 1110786747 2185707083
					 1069089445 406077387 203446459
					 3868416675 2875653677 4090169799
					 890516835 246160179 1672579861
					 4242444427 2062603291 793557843
					 3444933573 1802323335 3472899707
					 3188655691 3646873253 4149889019
					 3061854939 1888359715 503333261
					 2888002519 2659609347 2907644147
					 515861061 579734827 1918458496
					 2571339494 2764548314 642383838
					 269898520 597333446 1605001362
					 2699722410 3559800000 2585428982
					 3644288354 1949764494 3421994896
					 2750764078 3587445010 184081082
					 3405450064 1486346422 4023652090
					 50636894 2160348184 1976828486
					 2152286162 614069306 1583799648
					 428612214 1377190530 3944872510
					 549998704 3885312878 1586437026
					 2683794458 1038075520 1043973222
					 1056088186 1163017758 4113197624 7683142
					 3132872018 1636986346 2243070336
					 2666314166 2800321218 576835086
					 2122636144 2076647278 867898194
					 3404814106 194058128 1793529782
					 1972501498 3253420862 288319640
					 1884860454 2914337458 2684746170
					 1707666240 4060973782 485736834
					 3307952318 4235449776 2290953678
					 2533258146 2957288346 824121856
					 323905958 691856986 3641892126
					 3482776088 4274629958 2281365010
					 3387417066 1606109056 1479189686
					 3223033058 3111470798 3304567568
					 3566889646 3087629138 288533946
					 4203777872 3642609782 3886609786
					 2584116062 2164404248 1225860422
					 3906075602 978005818 406600544
					 3645555830 3180779458 282790846
					 593340336 1187811438 3406062626
					 2585587738 3286364992 1901695206
					 1813719162 691833118 3151560376
					 144251782 3207818130 2378676202
					 3859029760 3619050550 2564564290
					 402223246 2545718896 1818039342
					 3848190162 2420371738 2898471376
					 2293433014 2400076986 2516551102
					 4067211480 1299752614 719580850
					 1230907514 3944975296 2539376086
					 900099330 3158928382 4163742256
					 3924404174 2774796258 2420761562
					 3600180352 988294118 836425690
					 1529705694 1304226456 550753478
					 724443410 2851630506 882366400
					 3908295414 3305960418 3522158990
					 568811280 2006333742 2035246098
					 2303708218 3595490768 1305458614
					 2449318010 3140409566 3530306584
					 2853685830 1895944402 2691659322
					 615162464 3768611446 3924580226
					 3730995518 2521661552 640786286
					 3693885602 1827683482 1956778112
					 613201510 1771375482 2831891870
					 4087991480 3113514054 1880748114
					 3158756074 1507820288 133906870
					 1092947906 975458830 2812910448
					 1837782894 2667330514 1383521946
					 2997398928 1499105718 1034617082
					 3305875902 2305898136 1707886246
					 2014401586 3922637242 1176156224
					 338372566 1885102850 67464638 2143678512
					 3831501134 1459418274 178517658
					 3758970240 11588262 661874906 2231307294
					 2424601112 4229409990 2046931730
					 2609821418 714208640 3393407926
					 2752616162 3187378894 1886219792
					 851486254 142960466 2277713082
					 3584251728 3280757366 372419962
					 2518412638 3564815000 1231035846
					 3783671634 1826266810 755713504
					 922259446 1685010626 642708286
					 3779283120 4198421486 581130402
					 931331866 1324547392 3640239206
					 3041883130 2341342901 351573757
					 1870192537 3361900299 3912787753
					 3246114829 2803006785 64818723
					 2113948277 224985637 393417985
					 3206735835 890291145 1968293261 49537585
					 2587894507 1084089269 3020516733
					 858637337 3382211915 1819218073
					 2413262957 463238017 3851413955
					 2250001573 3838135813 3049789185
					 834272907 2987448489 420631293
					 4264408065 52894059 3611125813
					 3274014621 2215562137 3467538635
					 1956943913 4251593869 1926001985
					 3050418243 209119669 651068677 143505217
					 3561852251 3631033257 2351166669
					 4269369841 149556971 3766515989
					 2446468093 3289922489 188744811
					 2234069273 31079088 3970535662
					 1044597111))))

