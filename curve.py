import matplotlib.pyplot as plt
from matplotlib.pyplot import axes

import numpy as np
import sys
import math
from typing import List

from curve_math import solve_third_order_newton, solve_second_order_pos

class CurveSegment:    
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float):
        """
        Class constructor.

        Args:
            id (int): id of the segment
            t (float): duration of the segment
            vi (float): entry velocity of the segment
            ve (float): exit velocity of the segment
            ai (float): entry acceleration of the segment
            ae (float): exit acceleration of the segment
            j (float): jerk of the segment
            si (float): entry accumulated displacement of the segment
            se (float): exit accummulated displacement of the segment
            stpi (float): entry accumulated steps
            stpe (float): exit accummulated steps

        """
        self.id = id
        self.vi = vi
        self.ve = ve
        self.ai = ai
        self.ae = ae
        self.j  = j
        self.t  = t
        self.si = si
        self.se = se
        self.s = self.se - self.si

    def fa(self, t):
        return self.ai + self.j * (t)
    
    def fv(self, t):
        return self.vi + self.ai * t + (1 / 2) * self.j * (t**2)
    
    def fs(self,t):
        return self.si + self.vi * t + (1 / 2) * self.ai * (t**2) + (1 / 6) * self.j * (t**3)

class ContinuousCurveSegment(CurveSegment):
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float):
        super().__init__(id, t, vi, ve, ai, ae, j, si, se)

    def dump(self):
        print("|%10.4f|%10.4f| %10.4f| %10.4f| %10.4f| %10.4f|" % (self.t * 1000, self.vi, self.ve, self.ai, self.ae, self.s))

class DiscreteCurveSegment(CurveSegment):
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float, stpi: int, stpe: int, td: float = 0, t0: float = 0, s0: float = 0):
        super().__init__(id, t, vi, ve, ai, ae, j, si, se)

        self.stpi = stpi
        self.stpe = stpe

        self.td = td
        self.t0 = t0
        self.s0 = s0

    def dump(self):
        print("|%10.4f|%10.4f| %10.4f| %10.4f| %10.4f| %10.4f| %8d|" % (self.t * 1000, self.vi, self.ve, self.ai, self.ae, self.s, self.stpe - self.stpi + 1))
        
class Curve:
    def __init__(self, alpha: float):
        self.segments: List[ContinuousCurveSegment] = []
        self.discreteSegments: List[DiscreteCurveSegment] = []
        self.deltas = []

        self.epsilon = 1E-9
        self.alpha = alpha
        self.beta = 1 / self.alpha

        self.min_segment_steps = 2
        self.solve_error = 0.01
        self.debugBounds = True
        self.debugDiscretize = True
        self.hasSegment4 = False

    def __check_min_displacement__(self):
        # Get the minimum half displacement required for the curve
        min_s = self.__get_min_displacement__()

        for s in min_s:
            if (s < self.beta):
                return False

        return sum(min_s) <= (self.c.s / 2)

    def __get_min_steps__(self):
        # Get the minimum half displacement required for the curve
        min_s = self.__get_min_displacement__()

        # Create a list to store the required steps
        min_st = [sys.maxsize] * len(min_s)

        # Convert to steps
        for idx, s in enumerate(min_s):
            if (s > 0):
                min_st[idx] = math.floor(s * self.alpha)

        return min_st

    def __get_totals(self):
        t = 0
        stp = 0
        for segment in self.segments:
            t = t + segment.t
            stp = stp + segment.stp

        return [t, stp]

    def __discretize__(self):
        self.deltas = []

        # First step is generated at t = 0
        t = 0
        acc_s = 0
        step = 0

        v = 0

        beta = self.beta

        for segment in self.discreteSegments:
            # Get revelant data from segment
            vi = segment.vi
            ai = segment.ai
            j = segment.j

            td = segment.td
            t0 = segment.t0
            s0 = segment.s0

            if (self.debugDiscretize and (segment.id > 1)):
                print("------------------------------------------------------------------")

            if (td > 0) or (segment.t == 0):
                delta = td

                # Append delta
                self.deltas.append(delta)

                t = t0

                if (self.debugDiscretize):
                    print("step %4d, delta %f, t %f, s %f, v %f" % (step + 1, delta * 1000000000, t, segment.fs(t), segment.fv(t)))

                acc_s = s0
                step = step + 1
            else:
                t = 0
                acc_s = 0

            # For each, segment the inverse function of s(t) is defined using a lambda function fsInv(s), which 
            # is used later to discretize the velocity profile.
            if (segment.id == 1) or (segment.id == 3) or (segment.id == 5) or (segment.id == 7):
                fsInv = lambda s: solve_third_order_newton(j, 3 * ai, 6 * vi, -6 * (s + beta), x0, 10**-9)
                x0 = self.beta / vi

                def getDelta(s: float) -> float:
                    nonlocal t
                    nonlocal x0

                    prev_t = t
                    t = solve_third_order_newton(j, 3 * ai, 6 * vi, -6 * (s + beta), x0, 10**-9)
                    delta = t - prev_t

                    x0 = (beta / segment.fv(t)) + t                

                    return delta

            elif (segment.id == 2) or (segment.id == 6):
                fsInv = lambda s: solve_second_order_pos(0.5 * ai, vi, - (s + beta))

                def getDelta(s: float) -> float:
                    nonlocal t

                    prev_t = t
                    t = solve_second_order_pos(0.5 * ai, vi, - (s + beta))
                    delta = t - prev_t

                    return delta
            else:
                fsInv = lambda s: (s + beta) / vi

                def getDelta(s: float) -> float:
                    return beta / vi
            
            segment_step = 0

            while (step < segment.stpe):
                # Once a step is generated, we must find an increment of t (delta) which increases
                # the total displacement by beta.
                delta = getDelta(acc_s)

                # Append delta
                self.deltas.append(delta)

                if (self.debugDiscretize):
                    print("step %4d, delta %f, t %f, s %f, v %f" % (step + 1, delta * 1000000000, t, segment.fs(t), segment.fv(t)))


                step = step + 1
                segment_step = segment_step + 1
                acc_s = s0 + segment_step * beta

        print(sum(self.deltas))


    def addSegment(self, segment: CurveSegment):
        if isinstance(segment, ContinuousCurveSegment):
            self.segments.append(segment)
        elif isinstance(segment, DiscreteCurveSegment):
            self.discreteSegments.append(segment)

        if (segment.id == 4):
            self.hasSegment4 = True

    def solve(self):
        solved = self.__solve_motion_constraints__()

        if (self.c.t > 0):
            solved = self.__solve_time_and_motion_constraints__()

        if solved:
            if self.__bounds__():
                self.__discretize__()
            else:
                solved = False

        return solved
    
    def getTotalTime(self):
        t = 0

        for segment in self.segments:
            t = t + segment.t

        return t

    def getMaxVelocity(self):
        max_v = 0

        for segment in self.segments:
            if (segment.ve > max_v):
                max_v = segment.ve

        return max_v

    def getMaxAcceleration(self):
        max_a = 0

        for segment in self.segments:
            if (segment.ae > max_a):
                max_a = segment.ae

        return max_a
    
    def getProfile(self):
        return self.profile
    
    def continuousDump(self):
        t = 0
        s = 0

        print("-----------------------------------------------------------------------")
        print("| Continuous form                                                     |")
        print("-----------------------------------------------------------------------")
        print("| t        | vi       | ve        | ai        | ae        | s         |")
        print("-----------------------------------------------------------------------")

        for segment in self.segments:
            t = t + segment.t
            s = s + segment.s
            segment.dump()

        print("-----------------------------------------------------------------------")
        print("|%10.4f|          |           |           |           | %10.4f|" % (t, s))

    def discreteDump(self):
        t = 0
        s = 0
        stp = 0

        print("---------------------------------------------------------------------------------")
        print("| Discrete form                                                                 |")
        print("---------------------------------------------------------------------------------")
        print("| t        | vi       | ve        | ai        | ae        | s         | stp     |")
        print("---------------------------------------------------------------------------------")

        for segment in self.discreteSegments:
            t = t + segment.t
            s = s + segment.s
            stp = stp + segment.stpe - segment.stpi + 1
            segment.dump()

        print("---------------------------------------------------------------------------------")
        print("|%10.4f|          |           |           |           | %10.4f| %8d|" % (t, s, stp))

    def plotV(self, ax, samplingPoints = False):
        acc_t = 0
        for segment in self.discreteSegments:
            x_axis = np.linspace(0, segment.t)
            y_axis = segment.fv(x_axis)
            ax.plot(x_axis + acc_t, y_axis)
            ax.set_xlabel("time (s)", loc = "right")
            ax.set_ylabel("velocity (units/s)", loc = "center")
            acc_t = acc_t + segment.t      

        if samplingPoints:
            t = 0
            x_axis = [0]
            y_axis = [self.discreteSegments[0].fv(0)]

            for segment in self.discreteSegments:  
                segment_deltas = self.deltas[segment.stpi - 1:segment.stpe]

                if (segment.t0 > 0):
                    segment_deltas[0] = segment.t0

                segment_x = np.array(segment_deltas).cumsum()
                segment_y = segment.fv(segment_x)
                
                x_axis = x_axis + (segment_x + t).tolist()
                y_axis = y_axis + segment_y.tolist()

                t = t + segment.t

            ax.plot(x_axis, y_axis, marker=".",linestyle="")

    def plotA(self, ax, samplingPoints = False):
        acc_t = 0
        for segment in self.discreteSegments:
            x_axis = np.linspace(0, segment.t)
            y_axis = segment.fa(x_axis)
            ax.plot(x_axis + acc_t, y_axis)
            ax.set_xlabel("time (s)", loc = "right")
            ax.set_ylabel("acceleration (units/s^2)", loc = "center")
            acc_t = acc_t + segment.t

        if samplingPoints:
            t = 0
            x_axis = [0]
            y_axis = [self.discreteSegments[0].fa(0)]

            for segment in self.discreteSegments:  
                segment_deltas = self.deltas[segment.stpi - 1:segment.stpe]

                if (segment.t0 > 0):
                    segment_deltas[0] = segment.t0

                segment_x = np.array(segment_deltas).cumsum()
                segment_y = segment.fa(segment_x)
                
                x_axis = x_axis + (segment_x + t).tolist()
                y_axis = y_axis + segment_y.tolist()

                t = t + segment.t

            ax.plot(x_axis, y_axis, marker=".",linestyle="")

    def plotS(self, ax, samplingPoints = False):
        acc_t = 0
        for segment in self.discreteSegments:
            x_axis = np.linspace(0, segment.t)
            y_axis = segment.fs(x_axis)
            ax.plot(x_axis + acc_t, y_axis)            
            ax.set_xlabel("time (s)", loc = "right")
            ax.set_ylabel("displacement (units)", loc = "center")
            acc_t = acc_t + segment.t

        if samplingPoints:
            t = 0
            x_axis = [0]
            y_axis = [self.discreteSegments[0].fs(0)]

            for segment in self.discreteSegments:
                segment_deltas = self.deltas[segment.stpi - 1:segment.stpe]

                if (segment.t0 > 0):
                    segment_deltas[0] = segment.t0

                segment_x = np.array(segment_deltas).cumsum()
                segment_y = segment.fs(segment_x)
                
                x_axis = x_axis + (segment_x + t).tolist()
                y_axis = y_axis + segment_y.tolist()

                t = t + segment.t

            ax.plot(x_axis, y_axis, marker=".",linestyle="")
