from s_curve_full import SCurveFull
from motion_constraint import MotionConstraint

class Motion:
    def __init__(self, c:MotionConstraint, alpha:float):
        self.c = c
        self.alpha = alpha
        self.curve = None

    def simulate(self, full = True, partial = True):
        solved = False

        # First choice: full s-curve
        if (full):
            self.curve = SCurveFull(self.c, self.alpha)
            solved = self.curve.solve()

        if not solved:
            self.curve = None

        return solved
    
    def getCurve(self):
        return self.curve


