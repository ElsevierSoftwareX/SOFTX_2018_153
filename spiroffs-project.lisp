
(defun mappairn (fun list)
  "(mappairn fun list) => list

Accepts a two-argument function FUN that returns a list, and a list LIST.
Applies the function to each pair of adjacent elements in LIST (i.e., each
element shows up in two pairs, one with the element before and one with the
element after) and concatenates the results with NCONC."
  (let ((firsts list)
        (seconds (rest list)))
    (mapcan fun firsts seconds)))

(defun determinant (vector1 vector2)
  "(determinant vector1 vector2) => number

Treats VECTOR1 and VECTOR2 (which should each have 2 elements) as rows in a 2x2
matrix and finds the determinant."
  ;; Treats vector1 and vector2 as rows in a 2x2 matrix and finds the
  ;; determinant.
  (- (* (aref vector1 0) (aref vector2 1))
     (* (aref vector1 1) (aref vector2 0))))

(defun vector-scale (scalar vector)
  "(vector-scale scalar vector) => vector

Multiplies each element of VECTOR by SCALAR."
  (map 'vector
       (lambda (coordinate) (* scalar coordinate))
       vector))

(defun vector+ (&rest sequences)
  "(vector+ &rest sequences) => vector

Adds the corresponding elements of each sequence in SEQUENCES, and returns a
vector containing the results."
  (apply #'map 'vector #'+ sequences))

(defun vector= (vector-1 vector-2)
  "(vector= vector-1 vector-2) => boolean

Tests whether the corresponding elements of VECTOR-1 and VECTOR-2 (should each
have 2 elements) are equal."
  (and (eql (aref vector-1 0) (aref vector-2 0))
       (eql (aref vector-1 1) (aref vector-2 1))))

(defun shifts (basis-points exponents-left exponents-right)
  "(shifts basis-points exponents-left exponents-right) => list of lists

Accepts a Hilbert set BASIS-POINTS and two lists of exponents EXPONENTS-LEFT and
EXPONENTS-RIGHT. Returns the Hilbert set augmented with the intersections of
corresponding ideals (or, the amount by which each variable in the initial
inequalities has to be shifted to obtain the new inequalities when finding the
Hilbert-Kunz multiplicity)."
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
  "Represents a bounding inequality for the solid whose volume equals the
Hilbert-Kunz multiplicity of a given intersection ideal."
  (left-variable "" :type string)
  (right-variable "" :type string)
  (right-coefficient 0 :type number)
  (right-constant 0 :type number))

(defun make-axis-names (number)
  "(make-axis-names number) => list of axis names

Returns NUMBER names for axes beyond the x and y axis. Names will be of the form
z1, z2, z3, etc."
  (loop
     for i from 1 to number
     collect (format nil "z~d" i)))

(defun inequalities (exponents-left exponents-right axis-names)
  "(inequalities exponents-left exponents-right axis-names) => list of inequalities

Returns a list of the inequalities bounding the solid whose volume is equal to
the Hilbert-Kunz multiplicity of the intersection algebra of the ideals with
exponents EXPONENTS-LEFT and EXPONENTS-RIGHT."
  (flet ((inequalities% (x-coefficients
                         y-coefficients
                         x-shift
                         y-shift
                         &rest
                         z-shifts)
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
    (let* ((hilbert-set (hilbert-set exponents-left exponents-right))
           (shifts (shifts hilbert-set exponents-left exponents-right)))
      (mapcar (lambda (shift)
                (apply #'inequalities% exponents-left exponents-right shift))
              shifts))))

(defun hilbert-set (exponents-left exponents-right)
  "(hilbert-set exponents-left exponents-right) => list of vectors

Finds the Hilbert set for the division of the plane into cones corresponding to
the given lists of ideal exponents."
  (delete-duplicates (mappairn #'hilbert-basis
                               (append (list #(0 1))
                                       (mapcar #'vector
                                               exponents-right
                                               exponents-left)
                                       (list #(1 0))))
                     :test #'vector=))

(defun canonicalize-initial-points (point-left point-right)
  "(canonicalize-initial-points point-left point-right) => point-left, point-right

Ensures that the given points are suitable for use in the Hilbert basis
algorithm. Returns two points such that the coordinates of each are relatively
prime and (determinant point-left point-right) is greater than 0."
  (flet ((ensure-coprime (point)
           (vector-scale (/ (gcd (aref point 0)
                                 (aref point 1)))
                         point)))
    (let ((point-left (ensure-coprime point-left))
          (point-right (ensure-coprime point-right)))
      (if (< 0 (determinant point-left point-right))
          (values point-left point-right)
          (values point-right point-left)))))

(defun next-basis-point-guess (point)
  "(next-basis-point-guess point) => new-point

Given POINT, returns NEW-POINT such that (determinant point new-point) equals
1. Used as the initial \"guess\" for a new Hilbert basis point."
  (let ((b (aref point 0))
        (a (aref point 1)))
    (if (zerop b)
        (vector (/ a) 1)
        (loop
           for k from 1
           for ka-1 = (1- a) then (+ a ka-1)
           when (zerop (mod ka-1 b))
           return (vector k (/ ka-1 b))))))

(defun next-basis-point (point-left point-right)
  "(next-basis-point point-left point-right) => point

Returns a new point in the Hilbert basis that already contains POINT-LEFT and
POINT-RIGHT."
  (let* ((guess (next-basis-point-guess point-right))
         (scalar (ceiling (determinant guess point-left)
                          (determinant point-left point-right))))
    (vector+ guess (vector-scale scalar point-right))))

(defun hilbert-basis (point-left point-right)
  "(hilbert-basis point-left point-right) => list of points

Returns the Hilbert basis of the cone delineated by rays from the origin through
POINT-LEFT and POINT-RIGHT."
  (multiple-value-bind (point-left point-right)
      (canonicalize-initial-points point-left point-right)
    (if (vector= point-left point-right)
        '()
        (loop
           with so-far = (list point-left point-right)
           for first = (first so-far)
           for second = (second so-far)
           if (eql 1 (determinant first second))
           return so-far
           else
           do
             (push (next-basis-point first second) (cdr so-far))))))

;;; The following functions beginning with "write-" are all for writing the
;;; inequalities for the Hilbert-Kunz multiplicity in a format compatible with
;;; Mathematica.

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
  (format t
          "], {x, 0, 2000}, {y, 0, 2000}, ~{{~a, 0, 2000}~^, ~}]~%"
          axis-names))

(defun hkm-main (inputs)
  "(hkm-main inputs) => nil

Given a list of exponents of ideals (the first half of the list is assumed to be
the exponents of one ideal and the second half those of the other), outputs
Mathematica code to calculate the Hilbert-Kunz multiplicity and F-signature of
the corresponding intersection algebra. Destructively modifies the given list."
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
      (write-integration-range axis-names)
      (terpri)
      (f-signature exponents-left exponents-right))))

(defun f-signature (exponents-left exponents-right)
  "(f-signature exponents-left exponents-right) => nil

Prints out the inequalities bounding the solid whose volume is equal to the
f-signature of the intersection algebra of the ideals with exponents
EXPONENTS-LEFT and EXPONENTS-RIGHT in a format compatible with Mathematica."
  (loop
     initially (write-line "Integrate[Boole[")
     for remaining on exponents-left
     for coefficient-left = (car remaining)
     for coefficient-right in exponents-right
     for i from 1
     maximizing (max coefficient-left coefficient-right) into max
     do
       (format t
               "~dx <= z~d <= 1 + ~@*~dx &&~%"
               coefficient-left
               i)
       (format t
               "~dy <= z~d <= 1 + ~@*~dy"
               coefficient-right
               i)
     unless (endp (cdr remaining))
     do (write-line " && ")
     finally
       (let ((upper-bound (* 10 max)))
         (write-line "], {x, 0, 1}, {y, 0, 1}, ")
         (dotimes (j i)
           (format t
                   "{z~d, 0, ~d}~:[~;, ~]"
                   (1+ j)
                   upper-bound
                   (< (1+ j) i)))
         (write-string "]"))))
