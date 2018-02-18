(defun fan-order (list-a list-b)
  (let ((fan-order (sort (mapcar #'cons list-a list-b)
                         #'<
                         :key (lambda (cons)
                                (atan (cdr cons) (car cons))))))
    (values (mapcar #'car fan-order)
            (mapcar #'cdr fan-order))))

(defun mappair (fun seq)
  (typecase seq
    (list
     (let ((firsts seq)
           (seconds (rest seq)))
       (mapcar fun firsts seconds)))
    (vector
     (let ((firsts seq)
           (seconds (make-array (1- (length seq))
                                :displaced-to seq
                                :displaced-index-offset 1)))
       (map 'vector fun firsts seconds)))))

(defun shifts (basis-points exponents-left exponents-right)
  (flet ((shifts% (basis-point)
           (let ((basis-left (aref basis-point 0))
                 (basis-right (aref basis-point 1)))
             (list* basis-left
                    basis-right
                    (mapcar (lambda (i j)
                              (max (* basis-left i)
                                   (* basis-right j)))
                            exponents-left
                            exponents-right))))
         (default-shifts ()
           (let ((shift-length (+ 2 (list-length exponents-left))))
             (loop
                with as-vector = (make-array shift-length
                                             :element-type 'integer
                                             :initial-element 0)
                for i from 2 below shift-length
                do
                  (setf (aref as-vector (1- i)) 0
                        (aref as-vector i) 1)
                collect (coerce as-vector 'list)))))
    (nconc (default-shifts)
           (mapcar #'shifts% basis-points))))

(defstruct (inequality (:constructor make-inequality (left-variable right-variable right-coefficient right-constant)))
  (left-variable "" :type string)
  (right-variable "" :type string)
  (right-coefficient 0 :type number)
  (right-constant 0 :type number))

(defun make-axis-names (number)
  (loop
     for i from 1 to number
     collect (format nil "z~d" i)))

(defun inequalities (exponents-left exponents-right axis-names)
  (flet ((inequalities% (x-shift y-shift x-coefficients y-coefficients z-shifts)
           (delete-if
            #'null
            (list*
             (and (> x-shift 0)
                  (make-inequality "x" "" 0 x-shift))
             (and (> y-shift 0)
                  (make-inequality "y" "" 0 y-shift))
             (mapcan (lambda (z-shift x-coefficient y-coefficient axis-name)
                       (list (and (> z-shift (* x-coefficient x-shift))
                                  (make-inequality axis-name
                                                   "x"
                                                   x-coefficient
                                                   (- z-shift (* x-coefficient
                                                                 x-shift))))
                             (and (> z-shift (* y-coefficient y-shift))
                                  (make-inequality axis-name
                                                   "y"
                                                   y-coefficient
                                                   (- z-shift (* y-coefficient
                                                                 y-shift))))))
                     z-shifts
                     x-coefficients
                     y-coefficients
                     axis-names)))))
    (let* ((hilbert-basis (hilbert-basis exponents-left exponents-right))
           (shifts (shifts hilbert-basis exponents-left exponents-right)))
      (mapcar (lambda (shift)
                (destructuring-bind (x y &rest zs) shift
                  (inequalities% x y exponents-left exponents-right zs)))
              shifts))))

(define-condition hilbert-basis-parse-error (error)
  ((level :initarg :level
          :type symbol
          :reader level)
   (expected :initarg :expected
             :type list
             :reader expected)
   (datum :initarg :datum
          :reader datum))
  (:report (lambda (condition stream)
             (format stream
                     "While parsing ~s: Expected ~s but got ~s"
                     (level condition)
                     (expected condition)
                     (datum condition)))))

(defun hilbert-basis (exponents-left exponents-right)
  ;; The arguments are included for the future, in case I ever find out how to
  ;; get the Hilbert basis without reading it in from another program
  (declare (ignore exponents-left exponents-right))
  (let ((buffer (make-array 255 :element-type 'character :fill-pointer 0)))
    (labels ((parse-hilbert-basis (start)
               (unless (char= #\{ (char buffer start))
                 (error 'hilbert-basis-parse-error
                        :level 'parse-hilbert-basis
                        :expected '(#\{)
                        :datum (char buffer start)))
               (multiple-value-bind (points end) (parse-points (1+ start))
                 (unless (char= #\} (char buffer end))
                   (error 'hilbert-basis-parse-error
                          :level 'parse-hilbert-basis
                          :expected '(#\})
                          :datum (char buffer end)))
                 points))
             (parse-points (start)
               (loop
                  with point
                  with end
                  do
                    (multiple-value-setq (point end) (parse-point start))
                  collect point into points
                  until (char= #\} (char buffer end))
                  unless (char= #\, (char buffer end))
                  do
                    (error 'hilbert-basis-parse-error
                           :level 'parse-points
                           :expected '(#\, #\})
                           :datum (char buffer end))
                  end
                  do (setq start (1+ end))
                  finally (return (values points end))))
             (parse-point (start)
               (setq start (position-if-not #'ccl:whitespacep
                                            buffer
                                            :start start))
               (unless (char= #\{ (char buffer start))
                 (error 'hilbert-basis-parse-error
                        :level 'parse-point
                        :expected '(#\{)
                        :datum (char buffer start)))
               (multiple-value-bind (x end)
                   (read-from-string buffer t nil :start (1+ start))
                 (unless (typep x '(integer 0 *))
                   (error 'type-error
                          :expected '(integer 0 *)
                          :datum x))
                 (unless (char= #\, (char buffer end))
                   (error 'hilbert-basis-parse-error
                          :level 'parse-point
                          :expected '(#\,)
                          :datum (char buffer (1+ end))))
                 (multiple-value-bind (y end)
                     (read-from-string buffer t nil :start (1+ end) :end (position #\} buffer :start (1+ end)))
                   (unless (typep y '(integer 0 *))
                     (error 'type-error
                            :expected '(integer 0 *)
                            :datum y))
                   (unless (char= #\} (char buffer end))
                     (error 'hilbert-basis-parse-error
                            :level 'parse-point
                            :expected '(#\})
                            :datum (char buffer end)))
                   (values (vector x y) (1+ end))))))
      (loop       
         do
           (loop
              for char = (read-char)
              until (char= #\newline char)
              do (vector-push char buffer))
         until (string= "o5 = " buffer :end2 (min 5 (length buffer)))
         do (setf (fill-pointer buffer) 0)
         finally (return (parse-hilbert-basis 5))))))

(defun write-inequality (inequality &optional (stream *standard-output*))
  (let ((*standard-output* stream)) 
    (write-string (inequality-left-variable inequality))
    (write-string " < ")
    (cond
      ((zerop (inequality-right-coefficient inequality))
       (princ (inequality-right-constant inequality)))
      ((zerop (inequality-right-constant inequality))
       (princ (inequality-right-coefficient inequality))
       (write-string (inequality-right-variable inequality)))
      (t
       (princ (inequality-right-coefficient inequality))
       (write-string (inequality-right-variable inequality))
       (write-string (if (> (inequality-right-constant inequality) 0)
                         " + "
                         " - "))
       (princ (abs (inequality-right-constant inequality))))) 
    inequality))

(defun write-compound-inequality (inequalities)
  (loop       
     for (first-conjunction . remaining-conjunctions) on inequalities
     do
       (loop
          initially (write-string "((")
          for (first-disjunction . remaining-disjunctions) on first-conjunction
          do (write-inequality first-disjunction)
          unless (endp remaining-disjunctions)
          do (write-string " || ")
          finally (write-string "))"))
     unless (endp remaining-conjunctions)
     do (write-line " &&")))

(defun write-top-inequality (exponents-left exponents-right axis-names)
  (loop
     initially (format t "0 <= x && 0 <= y && ")
     for left in exponents-left
     for right in exponents-right
     for axis-name in axis-names
     do (format t "~dx <= ~a && ~dy <= ~a && " left axis-name right axis-name)
     finally (terpri)))

(defun write-integration-range (axis-names)
  (format t "], {x, 0, 2000}, {y, 0, 2000}, ~{{~a, 0, 2000}~^, ~}]" axis-names))

(defun hkm-main (inputs)
  (flet ((half (list)
           (loop
              for slow on list
              for fast on list by #'cddr
              finally (return slow)))
         (amputate (list sublist)
           (loop
              for last on list
              when (eql (rest last) sublist)
              do (setf (cdr last) nil)
              and return list)))
    (let* ((exponents-right (half inputs))
           (exponents-left (amputate inputs exponents-right))
           (axis-names (make-axis-names (list-length exponents-left))))
      (write-line "Integrate[Boole[")
      (write-top-inequality exponents-left exponents-right axis-names)
      (write-compound-inequality (inequalities exponents-left
                                               exponents-right
                                               axis-names))
      (write-integration-range axis-names))))

(defun f-signature (exponents-left exponents-right)
  (loop
     initially (write-line "Integrate[Boole[")
     for remaining on exponents-left
     for coefficient-left = (car remaining)
     for coefficient-right in exponents-right
     for i from 1
     for max = (max coefficient-left coefficient-right)
     then (max max coefficient-left coefficient-right)
     do
       (format t
               "~dx <= z~d < 1 + ~@*~dx &&~%"
               coefficient-left
               i)
       (format t
               "~dy <= z~d < 1 + ~@*~dy"
               coefficient-right
               i)
     unless (endp (cdr remaining))
     do (write-line " && ")
     finally
       (let ((upper-bound (* 10 max)))
         (format t "], {x, 0, ~d}, {y, 0, ~:*~d}, " upper-bound)
         (dotimes (j i)
           (format t
                   "{z~d, 0, ~d}~:[~;, ~]"
                   (1+ j)
                   upper-bound
                   (< (1+ j) i)))
         (write-string "]"))))
