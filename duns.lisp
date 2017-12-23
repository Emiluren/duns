(defun invert-vec (v)
  (destructuring-bind (x y z) v
	(list (- x) (- y) (- z))))

(defun magnitude (v)
  (destructuring-bind (x y z) v
	(sqrt (+ (* x x) (* y y) (* z z)))))

(defun square-magnitude (v)
  (destructuring-bind (x y z) v
	(+ (* x x) (* y y) (* z z))))

(defun normalize (v)
  (let ((l (magnitude v)))
	(destructuring-bind (x y z) v
	  (list (/ x l) (/ y l) (/ z l)))))

(defun scale-vec (v s)
  (destructuring-bind (x y z) v
	(list (* x s) (* y s) (* z s))))

(defun add-vec (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
	(destructuring-bind (x2 y2 z2) v2
	  (list (+ x1 x2) (+ y1 y2) (+ z1 z2)))))

(defun sub-vec (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
	(destructuring-bind (x2 y2 z2) v2
	  (list (- x1 x2) (- y1 y2) (- z1 z2)))))

(defun add-scaled-vec (v1 v2 s)
  (add-vec v1 (scale-vec v2 s)))

(defun component-product (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
	(destructuring-bind (x2 y2 z2) v2
	  (list (* x1 x2) (* y1 y2) (* z1 z2)))))

(defun scalar-product (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
	(destructuring-bind (x2 y2 z2) v2
	  (+ (* x1 x2) (* y1 y2) (* z1 z2)))))

(defun vector-product (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
	(destructuring-bind (x2 y2 z2) v2
	  (list (- (* y1 z2) (* z1 y2))
			(- (* z1 x2) (* x1 z2))
			(- (* x1 y2) (* y1 x2))))))

(defclass physical-particle ()
  ((particle-position
	:initarg :position
	:initform '(0 0 0)
	:accessor particle-position)
   (velocity
	:initarg :velocity
	:initform '(0 0 0)
	:accessor velocity)
   (acceleration
	:initarg :acceleration
	:initform '(0 0 0)
	:accessor acceleration)
   (damping
	:initarg :damping
	:initform 0.999
	:accessor damping)
   inverse-mass))

(defgeneric mass (particle))

(defmethod mass ((particle physical-particle))
  (with-slots (inverse-mass) particle
	(if (zerop inverse-mass)
		:positive-infinity
		(/ 1 inverse-mass))))

(defgeneric (setf mass) (value particle))

(defmethod (setf mass) (value (particle physical-particle))
  (setf (slot-value particle 'inverse-mass) (/ 1 value)))

(defmethod initialize-instance :after ((particle physical-particle)
									   &key mass)
  (if mass
	  (setf (mass particle) mass)
	  (setf (slot-value particle 'inverse-mass) 0)))

(defgeneric clear-accumulator (particle)
  (:documentation "Remove acting forces from a particle."))

(defmethod clear-accumulator ((particle physical-particle))
  )

(defgeneric integrate-particle (particle duration))

(defmethod integrate-particle ((particle physical-particle) duration)
  ;; Don't integrate things with infinite mass.
  (when (plusp (slot-value particle 'inverse-mass))
	(assert (plusp duration))
	(with-accessors ((pos particle-position)
					 (velocity velocity)
					 (acceleration acceleration)
					 (damping damping)) particle
	  ;; Update linear position.
	  (setf pos (add-scaled-vec pos velocity duration))
	  ;; Work out the acceleration from the force.
	  ;; (We'll add to this when we come to generate forces...)
	  (let ((resulting-acc acceleration))
		;; Update linear velocity from the acceleration.
		(setf velocity (add-scaled-vec velocity resulting-acc duration)))
	  ;; Impose drag.
	  (setf velocity (scale-vec velocity (expt damping duration)))
	  ;; Clear the forces
	  (clear-accumulator particle))))