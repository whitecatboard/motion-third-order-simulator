import sys
import unittest
import math

sys.path.append(".")
sys.path.append("ui")

from motion_constraint import MotionConstraint
from curve import Curve
from s_curve_full import SCurveFull
from s_curve_partial import SCurvePartial

class TestCurve(unittest.TestCase):
    def check_bounds_continuous(self, curve: Curve, constraints: MotionConstraint):
        # For each phase:
        #
        # - Check that the entry displacement is almost equal to the exit displacement of the previous phase
        # - Check that the entry velocity is almost equal to the exit velocity of the previous phase
        # - Check that the entry acceleration is almost equal to the exit acceleration of the previous phase
        #
        # - Check that the entry velocity is less or equal to the max velocity
        # - Check that the entry velocity is greater or equal to the initial velocity
        # - Check that the exit velocity is less or equal to the max velocity
        # - Check that the exit velocity is greater or equal to the initial velocity
        se = curve.phases[0].si
        ve = curve.phases[0].vi
        ae = curve.phases[0].ai

        for phase in curve.phases:
            si = phase.si
            vi = phase.vi
            ai = phase.ai

            self.assertAlmostEqual(si, se, None, "s bound check failded in phase %d" % (phase.id), self.s_epsilon)
            self.assertAlmostEqual(vi, ve, None, "v bound check failded in phase %d" % (phase.id), self.v_epsilon)
            self.assertAlmostEqual(ai, ae, None, "a bound check failded in phase %d" % (phase.id), self.a_epsilon)

            self.assertLessEqual(vi, constraints.v, "vi <= V check failded in phase %d" % (phase.id))
            self.assertGreaterEqual(vi, constraints.v0, "vi >= v0 check failded in phase %d" % (phase.id))
            self.assertLessEqual(ve, constraints.v, "ve <= V check failded in phase %d" % (phase.id))
            self.assertGreaterEqual(ve, constraints.v0, "ve >= v0 check failded in phase %d" % (phase.id))
            self.assertLessEqual(ai, constraints.a, "ai <= A check failded in phase %d" % (phase.id))
            self.assertLessEqual(ae, constraints.a, "ae <= A check failded in phase %d" % (phase.id))

            se = phase.se
            ve = phase.ve
            ae = phase.ae

    def check_bounds_discrete(self, curve: Curve, constraints: MotionConstraint):
        # For each phase:
        #
        # - Check that the entry displacement is almost equal to the exit displacement of the previous phase
        # - Check that the entry velocity is almost equal to the exit velocity of the previous phase
        # - Check that the entry acceleration is almost equal to the exit acceleration of the previous phase
        #
        # - Check that the entry velocity is less or equal to the max velocity
        # - Check that the entry velocity is greater or equal to the initial velocity
        # - Check that the exit velocity is less or equal to the max velocity
        # - Check that the exit velocity is greater or equal to the initial velocity
        se = curve.discretePhases[0].si
        ve = curve.discretePhases[0].vi
        ae = curve.discretePhases[0].ai

        for phase in curve.discretePhases:
            si = phase.si
            vi = phase.vi

            if (phase.id != 4):
                ai = phase.ai

            self.assertAlmostEqual(si, se, None, "s bound check failded in phase %d" % (phase.id), self.s_epsilon)
            self.assertAlmostEqual(vi, ve, None, "v bound check failded in phase %d" % (phase.id), self.v_epsilon)

            if (phase.id != 4):
                self.assertAlmostEqual(ai, ae, None, "a bound check failded in phase %d" % (phase.id), self.a_epsilon)
            else:
                self.assertEqual(phase.ai, 0, "ai not zero in phase %d" % (phase.id))
                self.assertEqual(phase.ae, 0, "ae not zero in phase %d" % (phase.id))

            self.assertLessEqual(vi, constraints.v, "vi <= V check failded in phase %d" % (phase.id))
            self.assertGreaterEqual(vi, constraints.v0, "vi >= v0 check failded in phase %d" % (phase.id))
            self.assertLessEqual(ve, constraints.v, "ve <= V check failded in phase %d" % (phase.id))
            self.assertGreaterEqual(ve, constraints.v0, "ve >= v0 check failded in phase %d" % (phase.id))
            self.assertLessEqual(ai, constraints.a, "ai <= A check failded in phase %d" % (phase.id))
            self.assertLessEqual(ae, constraints.a, "ae <= A check failded in phase %d" % (phase.id))
            
            se = phase.se
            ve = phase.ve

            if (phase.id != 4):
                if (phase.id == 3):
                    ae = -phase.ae
                else:
                    ae = phase.ae

    def check_discretize(self, curve: Curve):
        # - Check that the displacement between 2 consecutive steps is almost equal to beta
        #
        # - Check acceleration / deceleration:
        #   - When in accelerating phase, deltas are growing
        #   - When in decelererating phase, deltas are growing
        #   - Otherwise, deltas are constant
        step = 1
        idx = 0

        accelerating = False
        decelerating = False

        for phase in curve.discretePhases:
            acc_s = phase.si
            acc_t = 0

            accelerating = phase.id < 4
            decelerating = phase.id > 4

            while (step <= phase.stpe):
                delta_s = phase.fs(acc_t + curve.deltas[idx]) - acc_s
                self.assertAlmostEqual(delta_s, curve.beta, None, "s delta check failed for step %d" % (step), self.s_epsilon)

                if (idx > 0):
                    delta_diff = math.floor(curve.deltas[idx - 1] * 1000000000000) - math.floor(curve.deltas[idx] * 1000000000000)

                    if (accelerating):
                        self.assertGreaterEqual(delta_diff, 1, "accelerating phase check failed in phase %d (idx %d)" % (phase.id, idx))
                    elif (decelerating):                        
                        self.assertLessEqual(delta_diff, 1, "decelerating phase check failed in phase %d (idx %d)" % (phase.id, idx))
                    else:
                        self.assertGreaterEqual(delta_diff, 0, "constant velocity check failed in phase %d (idx %d)" % (phase.id, idx))

                acc_s = step * curve.beta
                acc_t = acc_t + curve.deltas[idx]

                step = step + 1
                idx = idx + 1

        # Check that the deltas are symmetrical along the entire curve
        idx_a = 0
        idx_b = len(curve.deltas) - 1

        while idx_a < idx_b:
            delta_diff = math.floor(curve.deltas[idx_a] * 1000000000000) - math.floor(curve.deltas[idx_a] * 1000000000000)
            self.assertEqual(delta_diff, 0, "delta symmetry check failed for idx %d - %d" % (idx_a, idx_b))

            idx_a = idx_a + 1
            idx_b = idx_b - 1

    def test_fit(self):
        self.s_epsilon = 0.0001
        self.v_epsilon = 0.0001
        self.a_epsilon = 0.0001
        self.t_epsilon = 0.0001

        cases = [
            [8.33, 100, 500, 10000, 000.0625, 400, True],
            [8.33, 100, 500, 10000, 000.0650, 400, True],
            [8.33, 100, 500, 10000, 000.0675, 400, True],
            [8.33, 100, 500, 10000, 001.0000, 400, True],
            [8.33, 100, 500, 10000, 001.0025, 400, True],
            [8.33, 100, 500, 10000, 001.0050, 400, True],
            [8.33, 100, 500, 10000, 010.0000, 400, True],
            [8.33, 100, 500, 10000, 010.0025, 400, True],
            [8.33, 100, 500, 10000, 010.0050, 400, True],
            [8.33, 100, 500, 10000, 100.0000, 400, True],
            [8.33, 100, 500, 10000, 100.0025, 400, True],
            [8.33, 100, 500, 10000, 100.0050, 400, True],
        ]

        for case in cases:
            constraint = MotionConstraint(case[0], case[1], case[2], case[3], case[4])        
            curve = SCurveFull(constraint, case[5])
            fit_done = curve.fit()
            self.assertEqual(fit_done, case[6], "fit check failed for %d displacement" % (case[4]))
            if (case[6]):
                self.check_bounds_continuous(curve, constraint)
                self.check_bounds_discrete(curve, constraint)
                self.check_discretize(curve)