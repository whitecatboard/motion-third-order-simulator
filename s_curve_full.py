import sys
import math
import numpy as np

from motion_constraint import MotionConstraint
from curve import Curve, CurveSegment, ContinuousCurveSegment, DiscreteCurveSegment
from curve_math import solve_third_order_newton, solve_second_order_newton, solve_second_order_pos, solve_second_order_neg, s_to_steps

class SCurveFull(Curve):
    def __init__(self, c: MotionConstraint, alpha: float):
        self.profile = "Full S-Curve"
        
        super().__init__(alpha)

        self.c = c
        self.steps = math.floor(c.s * self.alpha)

    def __get_min_displacement__(self):   
        """
        Get the minimum half displacement required for a full s-curve applying the displacement
        constraints. For a full s-curve, this is the displacement for segments 1, 2 & 3.

        Returns:
            list: a list with the displacement required for segment 1, 2 & 3.
        """

        # In segment 1 displacement depends on initial velocity, acceleration and jerk
        s1 = (((6 * self.c.a * self.c.j * self.c.v0) + self.c.a ** 3) / (6 * self.c.j ** 2))

        # In segment 2 displacement depends on initial velocity, target velocity, acceleration and jerk
        s2 = ((self.c.v ** 2 * self.c.j - self.c.a ** 2 * self.c.v - self.c.v0 ** 2 * self.c.j - self.c.a ** 2 * self.c.v0) / (2 * self.c.a * self.c.j))

        # In segment 3 displacement depends on target velocity, acceleration and jerk
        s3 = ((6 * self.c.v * self.c.a * self.c.j - self.c.a ** 3) / (6 * self.c.j ** 2))

        return [s1, s2, s3]
    
    def __bounds__(self):
        """
        Compute the bounds for each segment of the full s-curve. The bounds are computed for the
        continuous and the discrete form of the curve.
         
        For the continuous form, the following is computed for each segment:

        - Duration of the segment
        - Entry velocity of the segment
        - Exit velocity of the segment
        - Entry acceleration of the segment
        - Exit acceleration of the segment
        - Jerk of the segment
        - Entry accumulated displacement of the segment 
        - Exit accumulated displacement of the segment 

        For the discrete form, the following is computed for each segment:

        - Duration of the segment
        - Entry velocity of the segment
        - Exit velocity of the segment
        - Entry acceleration of the segment
        - Exit acceleration of the segment
        - Jerk of the segment
        - Entry accumulated displacement of the segment 
        - Exit accumulated displacement of the segment 
        - Entry accumulated steps of the segment 
        - Exit accumulated steps of the segment 
        """
        
        # Get total steps
        total_stp = math.floor(self.c.s * self.alpha)

        # Compute the displacement made in each segment, then discretize it in steps, and
        # finaly adjust the displacement to the real displacement made
        s1, s2, s3 = self.__get_min_displacement__()

        s1_stp = s_to_steps(s1, self.alpha, self.epsilon)
        s2_stp = s_to_steps(s2, self.alpha, self.epsilon)
        s3_stp = s_to_steps(s3, self.alpha, self.epsilon)

        s1_d = s1_stp * self.beta
        s2_d = s2_stp * self.beta
        s2_d = s2_stp * self.beta

        # Segment 1
        # ------------------------------------------------------------------------------

        # Segment characterization (continuous form)
        t = self.c.a / self.c.j
        ai = 0
        ae = self.c.a
        j = self.c.j
        vi = self.c.v0
        ve = self.c.v0 + ((self.c.a**2) / (2 * self.c.j))
        s  = s1

        si = 0
        se = s

        segment1 = ContinuousCurveSegment(1, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment1)

        # Segment characterization (discrete form)
        t_err = 0
        s_err = 0

        t1_err = 0
        s1_err = 0

        if (s1 - s1_d > self.epsilon):
            t1_err = t - solve_third_order_newton(j, 3 * ai, 6 * vi, - 6 * s1_d, t, 0.000001)
            s1_err = s1 - s1_d

            t_err = t1_err
            s_err = s1_err
            stp_err = math.ceil(s_err * self.alpha)

        stpi = 1
        stpe = s1_stp

        segment1d = DiscreteCurveSegment(1, t, vi, ve, ai, ae, j, si, se, stpi, stpe, 0, 0, 0)
        self.addSegment(segment1d)

        # Segment 2
        # ------------------------------------------------------------------------------

        # Segment characterization (continuous form)
        t = (self.c.v - self.c.v0 - (self.c.a ** 2 / self.c.j)) / self.c.a
        ai = ae
        ae = ai
        j = 0
        vi = ve
        ve = self.c.v - (self.c.a ** 2 / (2 * self.c.j))
        s  = s2

        si = se
        se = si + s

        segment2 = ContinuousCurveSegment(2, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment2)

        # Segment characterization (discrete form)
        t2_err = 0
        s2_err = 0
        s2_0 = 0
        t2_0 = 0
        t2_d = 0

        if (stp_err > 0):    
            # Check if there is enough space in segment to fix the accumulated displacement error      
            s_fix_err = stp_err * self.beta - s1_err
            if (s_fix_err <= s2):
                # Enough space
                t2_0 = solve_second_order_pos(0.5 * ai, vi, - s_fix_err)
                s2_0 = s_fix_err

                s2_stp = s_to_steps(s2 - s_fix_err, self.alpha, self.epsilon)
                s2_d = s_fix_err + s2_stp * self.beta
                s2_stp = s2_stp + stp_err

                t2_d = t_err + t2_0

                # Update accumulated error
                s_err = 0
                t_err = 0
                stp_err = 0
            else:
                # Not enough space

                # Update accumulated error
                s_err = s_err + s2
                t_err = t_err + t

        if (s2 - s2_d > self.epsilon):
            t2_err = t - solve_second_order_pos(0.5 * ai, vi, - s2_d)
            s2_err = s2 - s2_d

            t_err = t_err + t2_err
            s_err = s_err + s2_err
            stp_err = math.ceil(s_err * self.alpha)

        stpi = stpe + 1
        stpe = stpe + s2_stp

        segment2d = DiscreteCurveSegment(2, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t2_d, t2_0, s2_0)
        self.addSegment(segment2d)

        # Segment 3
        # ------------------------------------------------------------------------------

        # Segment characterization (continuous form)
        t = self.c.a / self.c.j
        ai = ae
        ae = 0
        j = -self.c.j
        vi = segment2.ve
        ve = self.c.v
        s = s3

        si = se
        se = si + s

        segment3 = ContinuousCurveSegment(3, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment3)

        # Segment characterization (discrete form)
        t3_err = 0
        s3_err = 0
        s3_0 = 0
        t3_0 = 0
        t3_d = 0

        if (stp_err > 0):    
            # Check if there is enough space in segment to fix the accumulated displacement error      
            s_fix_err = stp_err * self.beta - s_err
            if (s_fix_err <= s3):
                # Enough space
                t3_0 = solve_third_order_newton(j, 3 * ai, 6 * vi, - 6 * s_fix_err, t, 0.000001)
                s3_0 = s_fix_err
                
                s3_stp = s_to_steps(s3 - s_fix_err, self.alpha, self.epsilon)
                s3_d = s_fix_err + s3_stp * self.beta
                s3_stp = s3_stp + stp_err

                t3_d = t_err + t3_0

                # Update accumulated error
                s_err = 0
                t_err = 0
                stp_err = 0
            else:
                # Not enough space

                # Update accumulated error
                s_err = s_err + s3
                t_err = t_err + t

        if (s3 - s3_d > self.epsilon):
            t3_err = t - solve_third_order_newton(j, 3 * ai, 6 * vi, - 6 * s3_d, t, 0.000001)
            s3_err = s3 - s3_d

            t_err = t_err + t3_err
            s_err = s_err + s3_err
            stp_err = math.ceil(s_err * self.alpha)
       
        stpi = stpe + 1
        stpe = stpe + s3_stp

        segment3d = DiscreteCurveSegment(3, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t3_d, t3_0, s3_0)
        self.addSegment(segment3d)

        # Segment 3
        # ------------------------------------------------------------------------------

        # Segment characterization (continuous form)
        s4 = self.c.s - 2 * (s1 + s2 + s3)

        if (s4 > self.epsilon):
            s4_d = math.floor(s4 * self.alpha)

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

            segment4 = ContinuousCurveSegment(4, t, vi, ve, ai, ae, j, si, se)
            self.addSegment(segment4)

            # Segment characterization (discrete form)
            t4_err = 0
            s4_err = 0
            s4_0 = 0
            t4_0 = 0
            t4_d = 0

            if (stp_err > 0):    
                # Check if there is enough space in segment to fix the accumulated displacement error      
                s_fix_err = stp_err * self.beta - s_err
                if (2 * s_fix_err <= s4):
                    # Enough space
                    t4_0 = s_fix_err / vi
                    s4_0 = s_fix_err
                    
                    s4_stp = s_to_steps(s4 - 2 * s_fix_err, self.alpha, self.epsilon)
                    s4_d = s_fix_err + s4_stp * self.beta
                    s4_stp = s4_stp + stp_err

                    t4_d = t_err + t4_0

                    # Update accumulated error
                    s_err = 0
                    t_err = 0
                    stp_err = 0
                else:
                    # Not enough space

                    # Update accumulated error
                    s_err = s_err + s4
                    t_err = t_err + t

            if (s4 - s4_d > self.epsilon):
                t4_err = t - (s4_d / vi)
                s4_err = s4 - s4_d

                t_err = t_err + t4_err
                s_err = s_err + s4_err
                stp_err = math.ceil(s_err * self.alpha)
        
            stpi = stpe + 1
            stpe = stpe + s4_stp

            segment4d = DiscreteCurveSegment(4, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t4_d, t4_0, s4_0)
            self.addSegment(segment4d)

        # Segment 5 (symmetrical to segment 3)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = segment3.t
        ai = ae
        ae = -segment3.ai
        j  = segment3.j
        vi = ve
        ve = segment3.vi        
        s  = segment3.s

        si = se
        se = si + s

        segment5 = ContinuousCurveSegment(5, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment5)

        # Discrete form
        if (s3_err <= s3):
            s3_stp = s_to_steps(s3 - s3_err, self.alpha, self.epsilon)
            s3_stp = s3_stp + math.ceil(s3_err * self.alpha)

        stpi = stpe + 1
        stpe = stpe + s3_stp

        if (not self.hasSegment4):
            segment5d = DiscreteCurveSegment(5, t, vi, ve, ai, ae, j, si, se, stpi, stpe, 2 * t3_err, t3_err, s3_err)
        else:
            segment5d = DiscreteCurveSegment(5, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t3_err + t4_0, t3_err, s3_err)

        self.addSegment(segment5d)

        # Segment 6 (symmetrical to segment 2)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = segment2.t
        ai = -segment2.ae
        ae = -segment2.ae
        j  = segment2.j
        vi = segment2.ve
        ve = segment2.vi        
        s  = segment2.s

        si = se
        se = si + s

        segment6 = ContinuousCurveSegment(6, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment6)

        # Discrete form
        if (s2_err <= s2):
            s2_stp = s_to_steps(s2 - s2_err, self.alpha, self.epsilon)
            s2_stp = s2_stp + math.ceil(s2_err * self.alpha)

        stpi = stpe + 1
        stpe = stpe + s2_stp

        segment6d = DiscreteCurveSegment(6, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t3_0 + t2_err, t2_err, s2_err)
        self.addSegment(segment6d)

        # Segment 7 (symmetrical to segment 1)
        # ------------------------------------------------------------------------------

        # Continuous form
        t  = segment1.t
        ai = -segment1.ae
        ae = -segment1.ai
        j  = segment1.j
        vi = segment1.ve
        ve = segment1.vi        
        s  = segment1.s

        si = se
        se = si + s

        segment7 = ContinuousCurveSegment(7, t, vi, ve, ai, ae, j, si, se)
        self.addSegment(segment7)

        # Discrete form
        if (s1_err <= s1):
            s1_stp = s_to_steps(s1 - s1_err, self.alpha, self.epsilon)
            s1_stp = s1_stp + math.ceil(s1_err * self.alpha)

        stpi = stpe + 1
        stpe = stpe + s1_stp

        segment7d = DiscreteCurveSegment(7, t, vi, ve, ai, ae, j, si, se, stpi, stpe, t2_0 + t1_err, t1_err, s1_err)
        self.addSegment(segment7d)

        if (self.debugBounds):
            self.continuousDump()
            print()
            self.discreteDump()

        return True

    def __solve_time_and_motion_constraints_step__(self):
        solved = False

        min_v = self.c.v0 + (self.c.a ** 2 / self.c.j)

        v_ = solve_second_order_neg(
            self.c.j,
            - self.c.a * self.c.j * self.c.t - 2 * self.c.j * self.c.v0 + self.c.a**2,
            self.c.j * self.c.v0**2 - self.c.a**2 * self.c.v0 + self.c.a * self.c.j * self.c.s
        )

        if (not np.isnan(v_)) and (v_ > min_v) and (v_ <= self.c.v):
            self.c.update_v(v_)

            if (super().__check_min_displacement__()):
                solved = True

            self.c.restore_v()

        if (solved):
            self.c.update_v(v_)

        return solved

    def __solve_time_and_motion_constraints__(self):
        """
        Solve the time and motion constraints for a full s-curve, tuning-down the initial acceleration and velocity
        if necessary to satisfy the time and motion constraints.

        Returns:
            boolean: True if the time and motion constraints have been solved, False otherwise.
        """

        if (self.c.s <= self.c.v0 * self.c.t):
            return False
        
        solved = self.__solve_time_and_motion_constraints_step__()

        if (not solved):
            # Init the tune-down acceleration interval
            min_a = 0
            max_a = self.c.a  

            # Compute the iterations required for the given error
            its = math.ceil(math.log2((max_a - min_a) / self.solve_error))

            for _ in range(its):
                # Original formula:
                #   min_a + 0.5 * (max_a - min_a)
                #
                # Can be simplified to:
                #   0.5 * (min_a + max_a)
                self.c.update_a(0.5 * (min_a + max_a))

                if (self.__solve_time_and_motion_constraints_step__()):
                    min_a = self.c.a
                    min_v = self.c.v
                else:
                    max_a = self.c.a

                self.c.restore_a()

            if (min_a > 0):
                self.c.update_a(min_a)
                self.c.update_v(min_v)

                solved = True

        return solved
    
    def __solve_motion_constraints__(self):
        """
        Solve the motion constraints to a full s-curve, tuning-down the initial acceleration and velocity
        if necessary to satisfy the motion constraints.

        Returns:
            boolean: True if the motion constraints have been solved, False otherwise.
        """
        solved = False

        # Compute a_peak, which is the acceleration that makes the segment 2 not exist
        a_peak = math.sqrt(self.c.j * (self.c.v - self.c.v0))

        if (a_peak < self.c.a):
            # The current acceleration constraint makes the segment 2 not exist

            # Compute a new acceleration constraint that makes segment 2 exist with a displacement of a given steps
            a_ = solve_second_order_pos(self.c.v + self.c.v0, 2 * self.min_segment_steps * self.beta * self.c.j, -self.c.j * (self.c.v**2 - self.c.v0**2))

            self.c.update_a(a_)

        if (not super().__check_min_displacement__()):
            # As there are not enough space with the current constraints, we need to tune-down acceleration and
            # velocity to find an acceleration and a velocity that satisfy the constraints. To do this, we can
            # apply a bisect algorithm on acceleration.

            # Init the tune-down acceleration interval
            min_a = 0
            max_a = self.c.a  

            # Compute the iterations required for the given error
            its = math.ceil(math.log2((max_a - min_a) / self.solve_error))

            for _ in range(its):
                # Set a candidate acceleration witch is in half of the current interval
                #
                # Original formula:
                #   min_a + 0.5 * (max_a - min_a)
                #
                # Can be simplified to:
                #   0.5 * (min_a + max_a)
                self.c.update_a(0.5 * (min_a + max_a))

                # Get a candidate velocity that corresponds to the candidate acceleration
                v = solve_second_order_pos(1, self.c.a ** 2 / self.c.j, ((self.c.a ** 2 * self.c.v0) / self.c.j) - self.c.v0 ** 2 - self.c.s * self.c.a)

                if (not np.isnan(v) and (v > self.c.v0) and (v < self.c.v)):
                    # The candidate acceleration and velocity satisfy the acceleration and velocity constraints

                    # Check if the candidate acceleration and velocity satisfies the displacement constraints
                    self.c.update_v(v)

                    if (super().__check_min_displacement__()):
                        # Displacement constraints are satisfied, update acceleration interval
                        min_a = self.c.a
                        min_v = v
                    else:
                        # Displacement constraints are not satisfied, update acceleration interval
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