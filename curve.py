import matplotlib.pyplot as plt
from matplotlib.pyplot import axes

import numpy as np
import sys
import math
from typing import List

from curve_math import solve_third_order_newton, solve_second_order_pos

class CurvePhase:    
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float):
        """
        Class constructor.

        Args:
            id (int): id of the phase
            t (float): duration of the phase
            vi (float): entry velocity of the phase
            ve (float): exit velocity of the phase
            ai (float): entry acceleration of the phase
            ae (float): exit acceleration of the phase
            j (float): jerk of the phase
            si (float): entry accumulated displacement of the phase
            se (float): exit accummulated displacement of the phase
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

class ContinuousCurvePhase(CurvePhase):
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float):
        super().__init__(id, t, vi, ve, ai, ae, j, si, se)

    def dump(self):
        print("|%10.4f|%10.4f| %10.4f| %10.4f| %10.4f| %10.4f|" % (self.t * 1000, self.vi, self.ve, self.ai, self.ae, self.s))

class DiscreteCurvePhase(CurvePhase):
    def __init__(self, id: int, t:float, vi:float, ve:float, ai:float, ae:float, j:float, si:float, se:float, stpi: int, stpe: int):
        super().__init__(id, t, vi, ve, ai, ae, j, si, se)

        self.stpi = stpi
        self.stpe = stpe

    def dump(self):
        print("|%10.4f|%10.4f| %10.4f| %10.4f| %10.4f| %10.4f| %8d|" % (self.t * 1000, self.vi, self.ve, self.ai, self.ae, self.s, self.stpe - self.stpi + 1))
        
class Curve:
    def __init__(self, alpha: float):
        self.phases: List[ContinuousCurvePhase] = []
        self.discretePhases: List[DiscreteCurvePhase] = []
        self.deltas = []

        self.alpha = alpha
        self.beta = 1 / self.alpha

        self.debugBounds = False
        self.debugDiscretize = False

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

    def __check_min_steps__(self):
        # Get displacement for phases 1 & 3
        min_st = self.__get_min_steps__()

        min_s = 0

        for s in min_st:
            if (s != sys.maxsize) and (s > 0):
                min_s = min_s + s
            else:
                return False

        return (min_s <= (self.steps / 2))
    
    def __get_totals(self):
        t = 0
        stp = 0
        for phase in self.phases:
            t = t + phase.t
            stp = stp + phase.stp

        return [t, stp]
    
    def __discretize__(self):
        self.deltas = []

        # First step is generated at t = 0
        t = 0
        acc_s = 0
        beta = self.beta
        step = 0

        v = 0

        for phase in self.discretePhases:
            # Get revelant data from phase
            vi = phase.vi
            ai = phase.ai
            j = phase.j

            # For each, phase the inverse function of s(t) is defined using a lambda function fsInv(s), which 
            # is used later to discretize the velocity profile.
            if (phase.id == 1) or (phase.id == 3) or (phase.id == 5) or (phase.id == 7):
                fsInv = lambda s: solve_third_order_newton(j / 6, ai / 2, vi, - s - beta, x0, 0.000001)
                x0 = self.beta / vi
            elif (phase.id == 2) or (phase.id == 6):
                fsInv = lambda s: solve_second_order_pos(ai / 2, vi, - s - beta)
            else:
                fsInv = lambda s: (s + beta) / vi
            
            t = 0
            acc_s = 0
            phase_step = 0

            if (self.debugDiscretize and (phase.id > 1)):
                print("------------------------------------------------------------------")
            
            while (step < phase.stpe):
                # Once a step is generated, we must find an increment of t (delta) which increases
                # the total displacement by beta.
                prev_t = t
                t = fsInv(acc_s)
                delta = t - prev_t

                # Append delta
                self.deltas.append(delta)

                if (self.debugDiscretize):
                    print("step %4d, delta %f, t %f, s %f, v %f" % (step + 1, delta * 1000000000, t, phase.fs(t), phase.fv(t)))

                x0 = t
                step = step + 1
                phase_step = phase_step + 1
                acc_s = phase_step * beta

    def addPhase(self, phase: CurvePhase):
        if isinstance(phase, ContinuousCurvePhase):
            self.phases.append(phase)
        elif isinstance(phase, DiscreteCurvePhase):
            self.discretePhases.append(phase)

    def solve(self):
        if (self.c.t > 0):
            solved = self.__solve_time_and_motion_constraints__()
        else:
            solved = self.__solve_motion_constraints__()

        if solved:
            if self.__bounds__():
                self.__discretize__()
            else:
                solved = False

        return solved
    
    def getTotalTime(self):
        t = 0

        for phase in self.phases:
            t = t + phase.t

        return t

    def getMaxVelocity(self):
        max_v = 0

        for phase in self.phases:
            if (phase.ve > max_v):
                max_v = phase.ve

        return max_v

    def getMaxAcceleration(self):
        max_a = 0

        for phase in self.phases:
            if (phase.ae > max_a):
                max_a = phase.ae

        return max_a
    
    def continuousDump(self):
        t = 0
        s = 0

        print("-----------------------------------------------------------------------")
        print("| Continuous form                                                     |")
        print("-----------------------------------------------------------------------")
        print("| t        | vi       | ve        | ai        | ae        | s         |")
        print("-----------------------------------------------------------------------")

        for phase in self.phases:
            t = t + phase.t
            s = s + phase.s
            phase.dump()

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

        for phase in self.discretePhases:
            t = t + phase.t
            s = s + phase.s
            stp = stp + phase.stpe - phase.stpi + 1
            phase.dump()

        print("---------------------------------------------------------------------------------")
        print("|%10.4f|          |           |           |           | %10.4f| %8d|" % (t, s, stp))

    def plotV(self, ax, title):
        acc_t = 0
        for phase in self.discretePhases:
            x_axis = np.linspace(0, phase.t)
            y_axis = phase.fv(x_axis)
            ax.plot(x_axis + acc_t, y_axis)
            ax.title.set_text(title)
            ax.set_xlabel("t", loc = "right")
            ax.set_ylabel("v", loc = "center")
            acc_t = acc_t + phase.t

    def plotA(self, ax, title):
        acc_t = 0
        for phase in self.discretePhases:
            x_axis = np.linspace(0, phase.t)
            y_axis = phase.fa(x_axis)
            ax.plot(x_axis + acc_t, y_axis)
            ax.title.set_text(title)
            ax.set_xlabel("t", loc = "right")
            ax.set_ylabel("a", loc = "center")
            acc_t = acc_t + phase.t

    def plotS(self, ax, title):
        acc_t = 0
        for phase in self.discretePhases:
            x_axis = np.linspace(0, phase.t)
            y_axis = phase.fs(x_axis)
            ax.plot(x_axis + acc_t, y_axis)            
            ax.title.set_text(title)
            ax.set_xlabel("t", loc = "right")
            ax.set_ylabel("s", loc = "center")
            acc_t = acc_t + phase.t

        x_axis = []
        y_axis = []

        t = 0
        s = 0.0
        step = 0

        for delta in self.deltas:
            x_axis.append(t)
            y_axis.append(s)

            t = t + delta
            step = step + 1
            s = self.beta * step

        x_axis.append(t)
        y_axis.append(s)

        #ax.plot(x_axis, y_axis, marker="o")
