from motion_constraint import MotionConstraint
from s_curve_partial import SCurvePartial
from s_curve_full import SCurveFull

#constraint = MotionConstraint(8.33, 100, 500, 10000, 0.0625)
#constraint = MotionConstraint(8.33, 100, 500, 10000, 0.0650)

constraint = MotionConstraint(8.33, 100, 500, 10000, 10, 0.4)

curve = SCurvePartial(constraint, 400)
#curve = SCurveFull(constraint, 400)

#fit_done = part.fit()
#print(fit_done)

fit_done = curve.solve()
