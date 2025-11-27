from s_curve_full import SCurveFull
from s_curve_partial import SCurvePartial

from motion_constraint import MotionConstraint

class Motion:
    def __init__(self, c:MotionConstraint, alpha:float):
        self.c = c
        self.alpha = alpha
        self.curve = None

    def simulate(self):
        solved = False

        # First choice: full s-curve
        self.curve = SCurveFull(self.c, self.alpha)
        if (not self.curve.solve()):
            # Second choice: partial s-curve
            self.curve = SCurvePartial(self.c, self.alpha)
            if (self.curve.solve()):
                solved = True
        else:
            solved = True

        if not solved:
            self.curve = None

        return solved
    
    def getCurve(self):
        return self.curve


