import sys
import math
import numpy as np

from motion_constraint import MotionConstraint
from curve import Curve, CurvePhase, ContinuousCurvePhase, DiscreteCurvePhase
from curve_math import solve_third_order_newton, solve_second_order_newton, solve_second_order_pos, solve_second_order_neg

class SCurveFull(Curve):
    def __init__(self, c: MotionConstraint, alpha: float):
        super().__init__(alpha)

        self.c = c
        self.steps = math.floor(c.s * self.alpha)

    def __get_min_displacement__(self):   
        """
        Get the minimum half displacement required for a full s-curve applying the displacement
        constraints. For a full s-curve, this is the displacement for phases 1, 2 & 3.

        Returns:
            list: a list with the displacement required for phase 1, 2 & 3.
        """

        # In phase 1 displacement depends on initial velocity, acceleration and jerk
        s1 = (((6 * self.c.a * self.c.j * self.c.v0) + self.c.a ** 3) / (6 * self.c.j ** 2))

        # In phase 2 displacement depends on initial velocity, target velocity, acceleration and jerk
        s2 = ((self.c.v ** 2 * self.c.j - self.c.a ** 2 * self.c.v - self.c.v0 ** 2 * self.c.j - self.c.a ** 2 * self.c.v0) / (2 * self.c.a * self.c.j))

        # In phase 3 displacement depends on target velocity, acceleration and jerk
        s3 = ((6 * self.c.v * self.c.a * self.c.j - self.c.a ** 3) / (6 * self.c.j ** 2))

        return [s1, s2, s3]
    
    def __bounds(self):
        """
        Compute the bounds for each phase of the full s-curve. The bounds are computed for the
        continuous and the discrete form of the curve.
         
        For the continuous form, the following is computed for each phase:

        - Duration of the phase
        - Entry velocity of the phase
        - Exit velocity of the phase
        - Entry acceleration of the phase
        - Exit acceleration of the phase
        - Jerk of the phase
        - Entry accumulated displacement of the phase 
        - Exit accumulated displacement of the phase 

        For the discrete form, the following is computed for each phase:

        - Duration of the phase
        - Entry velocity of the phase
        - Exit velocity of the phase
        - Entry acceleration of the phase
        - Exit acceleration of the phase
        - Jerk of the phase
        - Entry accumulated displacement of the phase 
        - Exit accumulated displacement of the phase 
        - Entry accumulated steps of the phase 
        - Exit accumulated steps of the phase 
        """
        # Get total steps
        total_stp = math.floor(self.c.s * self.alpha)

        # Compute the displacement made in each phase, then discretize it in steps, and
        # finaly adjust the displacement to the real displacement made
        s1, s2, s3 = self.__get_min_displacement__()
        s1_stp, s2_stp, s3_stp = [math.floor(s * self.alpha) for s in [s1, s2, s3]]
        s1_d, s2_d, s3_d = [s * self.beta for s in [s1_stp, s2_stp, s3_stp]]

        # Phase 1
        # ------------------------------------------------------------------------------

        # Phase characterization (continuous form)
        t = self.c.a / self.c.j
        ai = 0
        ae = self.c.a
        j = self.c.j
        vi = self.c.v0
        ve = self.c.v0 + ((self.c.a**2) / (2 * self.c.j))
        s  = s1

        si = 0
        se = s

        phase1 = ContinuousCurvePhase(1, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase1)

        # Phase characterization (discrete form)
        aid = ai
        vid = vi

        # The discrete displacement is given by:
        #
        #   s(td) = (j / 6) * td^3 + (aid / 2) * td^2 + vid * td
        #
        # To know the discrete phase duration we need to solve s(td) = s1_d 
        td = solve_third_order_newton(j / 6, aid / 2, vid, - s1_d, t, 0.000001)

        # The discrete velocity is the first derivate of s(td)
        ved = vid + aid * td + (j / 2) * td ** 2

        # The discrete acceleration is the second derivate of s(td)
        aed = aid + j * td

        sd = s1_d

        sid = 0
        sed = sd

        stpi = 1
        stpe = s1_stp

        phase1d = DiscreteCurvePhase(1, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase1d)

        # Phase 2
        # ------------------------------------------------------------------------------

        # Phase characterization (continuous form)
        t = (self.c.v - self.c.v0 - (self.c.a ** 2 / self.c.j)) / self.c.a
        ai = ae
        ae = ai
        j = 0
        vi = ve
        ve = self.c.v - (self.c.a ** 2 / (2 * self.c.j))
        s  = s2

        si = se
        se = si + s

        phase2 = ContinuousCurvePhase(2, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase2)

        # Phase characterization (discrete form)
        aid = aed
        vid = ved

        # The discrete displacement is given by:
        #
        #   s(td) = (aid / 2) * td^2 + vid * td
        #
        # To know the discrete phase duration we need to solve s(td) = s2_d 
        td = solve_second_order_pos(aid / 2, vid, - s2_d)

        # The discrete velocity is the first derivate of s(td)
        ved = vid + aid * td

        # The discrete acceleration is the second derivate of s(td)
        aed = aid

        sd = s2_d

        sid = sed
        sed = sid + sd

        stpi = stpe + 1
        stpe = stpe + s2_stp

        phase2d = DiscreteCurvePhase(2, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase2d)

        # Phase 3
        # ------------------------------------------------------------------------------

        # Phase characterization (continuous form)
        t = self.c.a / self.c.j
        ai = ae
        ae = 0
        j = -self.c.j
        vi = phase2.ve
        ve = self.c.v
        s = s3

        si = se
        se = si + s

        phase3 = ContinuousCurvePhase(3, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase3)

        # Phase characterization (discrete form)
        aid = aed
        vid = ved

        # The discrete displacement is given by:
        #
        #   s(td) = (j / 6) * td^3 + (aid / 2) * td^2 + vid * td
        #
        # To know the discrete phase duration we need to solve s(td) = s3_d 
        td = solve_third_order_newton(j / 6, aid / 2, vid, - s3_d, t, 0.000001)

        # The discrete acceleration is the second derivate of s(td)
        aed = aid + j * td

        # In this phase, since the jerk is < 0, and due to the discretization, it is possible that
        # the exit acceleration is < 0, which is not possible because we are in the acceleration phase
        # of the s-curve. To solve this, the number of steps for this phase (which will later be added to phase 4)
        # are decremented until the acceleration is >= 0.
        while (aed < 0) and (s3_stp > 1):
            s3_stp = s3_stp - 1
            s3_d = s3_stp * self.beta
            td = solve_third_order_newton(j / 6, aid / 2, vid, - s3_d, t, 0.000001)
            aed = aid + j * td

        if (aed < 0):
            return False

        # The discrete velocity is the first derivate of s(td)
        ved = vid + aid * td + (j / 2) * td ** 2
        
        sd = s3_d

        sid = sed
        sed = sid + sd

        stpi = stpe + 1
        stpe = stpe + s3_stp

        phase3d = DiscreteCurvePhase(3, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase3d)

        # Phase 4
        # ------------------------------------------------------------------------------

        # Phase characterization (continuous form)
        s4 = self.c.s - 2 * (s1 + s2 + s3)

        if (s4 > 0):
            # Continuous form
            ai = 0
            ae = 0
            j = 0
            vi = ve
            ve = vi
            s = s4
            t = s / vi

            si = se
            se = si + s

            phase4 = ContinuousCurvePhase(4, t, vi, ve, ai, ae, j, si, se,)
            self.addPhase(phase4)

        # Phase characterization (discrete form)
        stp = total_stp - 2 * (s1_stp + s2_stp + s3_stp)

        if (stp > 0):
            aid = 0
            aed = 0
            j = 0
            vid = ved
            ved = vid
            sd = stp * self.beta
            t = sd / vid
            td = t

            sid = sed
            sed = sid + sd

            stpi = stpe + 1
            stpe = stpe + stp

            phase4d = DiscreteCurvePhase(4, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
            self.addPhase(phase4d)

        # Phase 5 (symmetrical to phase 3)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = phase3.t
        ai = ae
        ae = -phase3.ai
        j  = phase3.j
        vi = ve
        ve = phase3.vi        
        s  = phase3.s

        si = se
        se = si + s

        phase5 = ContinuousCurvePhase(5, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase5)

        # Discrete form
        td  = phase3d.t
        aid = -phase3d.ae
        aed = -phase3d.ai
        vid = phase3d.ve
        ved = phase3d.vi        
        sd  = phase3d.s

        sid = sed
        sed = sid + sd

        stpi = stpe + 1
        stpe = stpe + s3_stp

        phase5d = DiscreteCurvePhase(5, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase5d)

        # Phase 6 (symmetrical to phase 2)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = phase2.t
        ai = -phase2.ae
        ae = -phase2.ae
        j  = phase2.j
        vi = phase2.ve
        ve = phase2.vi        
        s  = phase2.s

        si = se
        se = si + s

        phase6 = ContinuousCurvePhase(6, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase6)

        # Discrete form
        td  = phase2d.t
        aid = -phase2d.ae
        aed = -phase2d.ae
        vid = phase2d.ve
        ved = phase2d.vi        
        sd  = phase2d.s

        sid = sed
        sed = sid + sd

        stpi = stpe + 1
        stpe = stpe + s2_stp

        phase6d = DiscreteCurvePhase(6, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase6d)

        # Phase 7 (symmetrical to phase 1)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = phase1.t
        ai = -phase1.ae
        ae = -phase1.ai
        j  = phase1.j
        vi = phase1.ve
        ve = phase1.vi        
        s  = phase1.s

        si = se
        se = si + s

        phase7 = ContinuousCurvePhase(7, t, vi, ve, ai, ae, j, si, se)
        self.addPhase(phase7)

        # Discrete form
        td  = phase1d.t
        aid = -phase1d.ae
        aed = -phase1d.ai
        vid = phase1d.ve
        ved = phase1d.vi        
        sd  = phase1d.s

        sid = sed
        sed = sid + sd

        stpi = stpe + 1
        stpe = stpe + s1_stp

        phase7d = DiscreteCurvePhase(7, td, vid, ved, aid, aed, j, sid, sed, stpi, stpe)
        self.addPhase(phase7d)

        if (self.debugBounds):
            self.continuousDump()
            print()
            self.discreteDump()

        return True

    def __solve_time_and_motion_constraints(self):
        """
        Solve the time and motion constraints for a full s-curve, tuning-down the initial acceleration and velocity
        if necessary to satisfy the time and motion constraints.

        Returns:
            boolean: True if the time and motion constraints have been solved, False otherwise.
        """

        # If we have a time constraint, we need to tune-down acceleration and velocity to find an acceleration
        # and a velocity that satisfy the motion constraints, and mades the motion duration equal to the time
        # constraint. To do this, we can apply a bisect algorithm on velocity.

        solved = False
        
        # Init the tune-down velocity interval
        min_v = self.c.v0
        max_v = self.c.v 

        while (abs(max_v - min_v) > 1):
            # Set a candidate velocity witch is in half of the current interval
            self.c.update_v(min_v + ((max_v - min_v) / 2))

            # Get a candidate acceleration that corresponds to the candidate velocity
            a = solve_second_order_newton(
                self.c.v - self.c.v0,
                self.c.j * (self.c.s - self.c.v * self.c.t),
                self.c.j * (self.c.v0 ** 2 - 2 * self.c.v * self.c.v0 + self.c.v ** 2),
                self.c.a,
                0.000001)                
                        
            if (not np.isnan(a) and (a <= self.c.a) and (a > 0)):
                # The candidate acceleration and velocity satisfy the acceleration and velocity constraints

                # Check if the candidate acceleration and velocity satisfies the displacement constraints
                self.c.update_a(a)

                if (super().__check_min_steps__()):
                    # Acceleration, velocity, and displacement constraints are satisfied
                    # Update velocity interval
                    min_v = self.c.v
                    min_a = a
                else:   
                    # Acceleration, velocity, and displacement constraints are not satisfied
                    # Update velocity interval
                    max_v = self.c.v

                self.c.restore_a()

            else:
                # The candidate acceleration and velocity don't satisfy the acceleration and velocity constraints.
                # Update velocity interval
                min_v = self.c.v

            self.c.restore_v()

        if (min_v > self.c.v0):
            self.c.update_a(min_a)
            self.c.update_v(min_v)
            solved = True

        return solved
    
    def __solve_motion_constraints(self):
        """
        Solve the motion constraints to a full s-curve, tuning-down the initial acceleration and velocity
        if necessary to satisfy the motion constraints.

        Returns:
            boolean: True if the motion constraints have been solved, False otherwise.
        """
        solved = False

        if (not super().__check_min_steps__()):
            # As there are not enough space with the current constraints, we need to tune-down acceleration and
            # velocity to find an acceleration and a velocity that satisfy the constaints. To do this, we can
            # apply a bisect algorithm on acceleration.

            # Init the tune-down acceleration interval
            min_a = 0
            max_a = self.c.a  

            while (abs(max_a - min_a) > 0.01):
                # Set a candidate acceleration witch is in half of the current interval
                self.c.update_a(min_a + ((max_a - min_a) / 2))

                # Get a candidate velocity that corresponds to the candidate acceleration
                v = solve_second_order_pos(1, self.c.a ** 2 / self.c.j, ((self.c.a ** 2 * self.c.v0) / self.c.j) - self.c.v0 ** 2 - self.c.s * self.c.a)
                if (not np.isnan(v) and (v > self.c.v0) and (v < self.c.v)):
                    # The candidate acceleration and velocity satisfy the acceleration and velocity constraints

                    # Check if the candidate acceleration and velocity satisfies the displacement constraints
                    self.c.update_v(v)

                    if (super().__check_min_steps__()):
                        # Acceleration, velocity, and displacement constraints are satisfied
                        # Update acceleration interval
                        min_a = self.c.a
                        min_v = v
                    else:
                        # Acceleration, velocity, and displacement constraints are not satisfied
                        # Update acceleration interval
                        max_a = self.c.a

                    self.c.restore_v()                
                else:
                    # The candidate acceleration and velocity don't satisfy the acceleration and velocity constraints.
                    # Update acceleration interval
                    max_a = self.c.a

                self.c.restore_a()                

            if (min_a > 0):
                self.c.update_a(min_a)
                self.c.update_v(min_v)
                solved = True
        else:
            solved = True
            
        return solved
    
    def solve(self):
        if (self.c.t > 0):
            solved = self.__solve_time_and_motion_constraints()
        else:
            solved = self.__solve_motion_constraints()

        if solved:
            if self.__bounds():
                self.__discretize__()
            else:
                solved = False

        return solved
