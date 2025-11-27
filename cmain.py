from motion_constraint import MotionConstraint
from motion import Motion

constraint = MotionConstraint(8.33, 100, 500, 10000, 0.0025 * 5)

motion = Motion(constraint, 400)

motion.simulate()
