from motion_constraint import MotionConstraint
from motion import Motion

#constraint = MotionConstraint(8.33, 100, 500, 10000, 25, 0.7139)
constraint = MotionConstraint(8.33, 100, 500, 10000, 26)
motion = Motion(constraint, 400)

motion.simulate(True, False)
