import sys
import math
import numpy as np

from motion_constraint import MotionConstraint
from curve import Curve, CurvePhase, ContinuousCurvePhase, DiscreteCurvePhase
from curve_math import solve_third_order_newton, solve_second_order_pos

class SCurvePartial(Curve):
    def __init__(self, c: MotionConstraint, alpha: float):
        self.profile = "Partial S-Curve"
        
        super().__init__(alpha)

        self.c = c
        self.steps = math.floor(c.s * self.alpha)

    def __get_min_displacement__(self):   
        """
        Get the minimum half displacement required for a partial s-curve applying the displacement
        constraints. For a partial s-curve, this is the displacement for phases 1 & 3.

        Returns:
            list: a list with the displacement required for phase 1 & 3.
        """

        # In phase 1 & 3 displacement depends on initial velocity, acceleration and jerk
        s1 = (((6 * self.c.a * self.c.j * self.c.v0) + self.c.a ** 3) / (6 * self.c.j ** 2))
        s3 = (((6 * self.c.a * self.c.j * self.c.v0) + (5 * self.c.a ** 3)) / (6 * self.c.j ** 2))

        return [s1, s3]
    
    def __bounds__(self):
        """
        Compute the bounds for each phase of the partial s-curve. The bounds are computed for the
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
        s1, s3 = self.__get_min_displacement__()
        s1_stp, s3_stp = [math.floor(s * self.alpha) for s in [s1, s3]]
        s1_d, s3_d = [s * self.beta for s in [s1_stp, s3_stp]]

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

        # Phase 3
        # ------------------------------------------------------------------------------

        # Phase characterization (continuous form)
        t = self.c.a / self.c.j
        ai = ae
        ae = 0
        j = -self.c.j
        vi = phase1.ve
        ve = self.c.v0 + ((self.c.a ** 2) / self.c.j)
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
        s4 = self.c.s - 2 * (s1 + s3)        

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

            phase4 = ContinuousCurvePhase(4, t, vi, ve, ai, ae, j, si, se)
            self.addPhase(phase4)

        # Phase characterization (discrete form)
        stp = total_stp - 2 * (s1_stp + s3_stp)

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

    def __solve_time_and_motion_constraints__(self):
        """
        Solve the time and motion constraints for a partial s-curve, tuning-down the initial acceleration
        if necessary to satisfy the time and motion constraints.

        Returns:
            boolean: True if the time and motion constraints have been solved, False otherwise.
        """

        solved = False

        # Get a candidate acceleration that satisfy the time constraint
        a = solve_third_order_newton(
            2, 
            - self.c.t * self.c.j,
            0,
            self.c.j ** 2 * (self.c.s - self.c.t * self.c.v0),
            self.c.a  ,
            0.000001
        )

        if ((not np.isnan(a)) and (a < self.c.a) and (a > 0)):
            # Get a candidate velocity that corresponds to the candidate acceleration
            v = self.c.v0 + ((a ** 2) / self.c.j)
    
            if ((v > self.c.v0) and (v < self.c.v)):
                # Check if the candidate acceleration and velocity satisfies the displacement constraints
                self.c.update_a(a)
                self.c.update_v(v)

                if (super().__check_min_steps__()):
                    # Displacement constraints are satisfied
                    solved = True
                else:   
                    # Displacement constraints are not satisfied
                    self.c.restore_a()
                    self.c.restore_v()

        return solved

    def __solve_motion_constraints__(self):
        """
        Solve the motion constraints to a partial s-curve, tuning-down the initial acceleration and velocity
        if necessary to satisfy the motion constraints.

        Returns:
            boolean: True if the motion constraints have been solved, False otherwise.
        """

        solved = False

        if (not super().__check_min_steps__()):
            # Not enough space for the current constraints. As required space for a partial s-curve
            # depends only on initial velocity, acceleration and jerk, we can try to tune-down acceleration
            # maintaining invariant the initial velocity and jerk.

            # We need to solve the following third order equation
            fa = (1 / self.c.j ** 2)
            fb = 0
            fc = (2 * self.c.v0) / self.c.j
            fd = - 1 * (self.steps / 2) * self.beta
            x0 = self.c.a

            a = solve_third_order_newton(fa, fb, fc, fd, x0, 0.000001)
            if (not np.isnan(a)):
                v = self.c.v0 + ((a**2)/self.c.j)

                if (v > self.c.v):
                    a = math.sqrt(self.c.j * (self.c.v - self.c.v0))
                    v = self.c.v

            self.c.update_a(a)
            self.c.update_v(v)

            if (not super().__check_min_steps__()):
                self.c.restore_a()
                self.c.restore_v()                
            else:
                solved = True
        else:
            v = self.c.v0 + ((self.c.a**2)/self.c.j)
            self.c.update_v(v)

            solved = True

        return solved