import sys
import unittest
import math

sys.path.append(".")
sys.path.append("ui")

from motion_constraint import MotionConstraint
from curve import Curve
from s_curve_full import SCurveFull

class TestCurve(unittest.TestCase):
    def check_bounds_continuous(self, curve: Curve, constraints: MotionConstraint):
        # For each segment:
        #
        # - Check that the entry displacement is almost equal to the exit displacement of the previous segment
        # - Check that the entry velocity is almost equal to the exit velocity of the previous segment
        # - Check that the entry acceleration is almost equal to the exit acceleration of the previous segment
        #
        # - Check that the entry velocity is less or equal to the max velocity
        # - Check that the entry velocity is greater or equal to the initial velocity
        # - Check that the exit velocity is less or equal to the max velocity
        # - Check that the exit velocity is greater or equal to the initial velocity
        se = curve.segments[0].si
        ve = curve.segments[0].vi
        ae = curve.segments[0].ai

        for segment in curve.segments:
            si = segment.si
            vi = segment.vi
            ai = segment.ai

            self.assertAlmostEqual(si, se, None, "s bound check failded in segment %d" % (segment.id), self.s_epsilon)
            self.assertAlmostEqual(vi, ve, None, "v bound check failded in segment %d" % (segment.id), self.v_epsilon)
            self.assertAlmostEqual(ai, ae, None, "a bound check failded in segment %d" % (segment.id), self.a_epsilon)

            self.assertLessEqual(vi, constraints.v, "vi <= V check failded in segment %d" % (segment.id))
            self.assertGreaterEqual(vi, constraints.v0, "vi >= v0 check failded in segment %d" % (segment.id))
            self.assertLessEqual(ve, constraints.v, "ve <= V check failded in segment %d" % (segment.id))
            self.assertGreaterEqual(ve, constraints.v0, "ve >= v0 check failded in segment %d" % (segment.id))
            self.assertLessEqual(ai, constraints.a, "ai <= A check failded in segment %d" % (segment.id))
            self.assertLessEqual(ae, constraints.a, "ae <= A check failded in segment %d" % (segment.id))

            se = segment.se
            ve = segment.ve
            ae = segment.ae

    def check_bounds_discrete(self, curve: Curve, constraints: MotionConstraint):
        # For each segment:
        #
        # - Check that the entry displacement is almost equal to the exit displacement of the previous segment
        # - Check that the entry velocity is almost equal to the exit velocity of the previous segment
        # - Check that the entry acceleration is almost equal to the exit acceleration of the previous segment
        #
        # - Check that the entry velocity is less or equal to the max velocity
        # - Check that the entry velocity is greater or equal to the initial velocity
        # - Check that the exit velocity is less or equal to the max velocity
        # - Check that the exit velocity is greater or equal to the initial velocity
        se = curve.discreteSegments[0].si
        ve = curve.discreteSegments[0].vi
        ae = curve.discreteSegments[0].ai

        for segment in curve.discreteSegments:
            si = segment.si
            vi = segment.vi

            if (segment.id != 4):
                ai = segment.ai

            self.assertAlmostEqual(si, se, None, "s bound check failded in segment %d" % (segment.id), self.s_epsilon)
            self.assertAlmostEqual(vi, ve, None, "v bound check failded in segment %d" % (segment.id), self.v_epsilon)

            if (segment.id != 4):
                self.assertAlmostEqual(ai, ae, None, "a bound check failded in segment %d" % (segment.id), self.a_epsilon)
            else:
                self.assertEqual(segment.ai, 0, "ai not zero in segment %d" % (segment.id))
                self.assertEqual(segment.ae, 0, "ae not zero in segment %d" % (segment.id))

            self.assertLessEqual(vi, constraints.v, "vi <= V check failded in segment %d" % (segment.id))
            self.assertGreaterEqual(vi, constraints.v0, "vi >= v0 check failded in segment %d" % (segment.id))
            self.assertLessEqual(ve, constraints.v, "ve <= V check failded in segment %d" % (segment.id))
            self.assertGreaterEqual(ve, constraints.v0, "ve >= v0 check failded in segment %d" % (segment.id))
            self.assertLessEqual(ai, constraints.a, "ai <= A check failded in segment %d" % (segment.id))
            self.assertLessEqual(ae, constraints.a, "ae <= A check failded in segment %d" % (segment.id))
            
            se = segment.se
            ve = segment.ve

            if (segment.id != 4):
                if (segment.id == 3):
                    ae = -segment.ae
                else:
                    ae = segment.ae

    def check_discretize(self, curve: Curve):
        # - Check that the displacement between 2 consecutive steps is almost equal to beta
        #
        # - Check acceleration / deceleration:
        #   - When in accelerating segment, deltas are growing
        #   - When in decelererating segment, deltas are growing
        #   - Otherwise, deltas are constant
        step = 1
        idx = 0

        s_err = 0
        t_err = 0

        accelerating = False
        decelerating = False

        for segment in curve.discreteSegments:
            acc_s = segment.si
            acc_t = 0
            s0 = segment.s0
            segment_stp = 1

            while (step <= segment.stpe):
                if (segment.id in [4]) and (segment.s0 > 0) and (segment_stp == 1):                    
                    accelerating = False
                    decelerating = True
                elif (segment.id in [5]) and (segment.s0 > 0) and (segment_stp == 1) and (not curve.hasSegment4):
                    accelerating = True
                    decelerating = False
                else:
                    accelerating = segment.id < 4
                    decelerating = segment.id > 4

                delta_s = segment.fs(acc_t + curve.deltas[idx] - t_err) - acc_s + s_err
                self.assertAlmostEqual(delta_s, curve.beta, None, "s delta check failed for step %d" % (step), self.s_epsilon)

                if (idx > 0):
                    delta_diff = math.floor(curve.deltas[idx - 1] * 1000000000000) - math.floor(curve.deltas[idx] * 1000000000000)

                    if (segment.id == 4) and (step > 5335):
                        print("diff %.10f, prev %.10f, curr %.10f" % (delta_diff, curve.deltas[idx - 1] * 1000000000000, curve.deltas[idx] * 1000000000000))

                    if (accelerating):
                        self.assertGreaterEqual(delta_diff, 1, "accelerating segment check failed in segment %d (idx %d)" % (segment.id, idx))
                    elif (decelerating):                        
                        self.assertLessEqual(delta_diff, 1, "decelerating segment check failed in segment %d (idx %d)" % (segment.id, idx))
                    else:
                        self.assertGreaterEqual(delta_diff, 0, "constant velocity check failed in segment %d (idx %d)" % (segment.id, idx))

                if (s0 > 0):
                    acc_s = segment.si + s0 + (segment_stp - 1) * curve.beta
                else:
                    acc_s = step * curve.beta

                acc_t = acc_t + curve.deltas[idx] - t_err

                s_err = 0
                t_err = 0

                step = step + 1
                segment_stp = segment_stp + 1
                idx = idx + 1

            s_err = segment.s - (segment.fs(acc_t) - segment.si)
            t_err = segment.t - acc_t

            print("s_err %.10f, t_err %.10f" % (s_err, t_err))

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
            [8.33, 100, 500, 10000, 000.8000, 400, True],
            [8.33, 100, 500, 10000, 001.0000, 400, True],
            [8.33, 100, 500, 10000, 001.0025, 400, True],
            [8.33, 100, 500, 10000, 001.0050, 400, True],
            [8.33, 100, 500, 10000, 010.0000, 400, True],
            [8.33, 100, 500, 10000, 010.0025, 400, True],
            [8.33, 100, 500, 10000, 010.0050, 400, True],
            [8.33, 100, 500, 10000, 026.0000, 400, True],
            [8.33, 100, 500, 10000, 100.0000, 400, True],
            [8.33, 100, 500, 10000, 100.0025, 400, True],
            [8.33, 100, 500, 10000, 100.0050, 400, True],
        ]

        for case in cases:
            constraint = MotionConstraint(case[0], case[1], case[2], case[3], case[4])        
            curve = SCurveFull(constraint, case[5])
            fit_done = curve.solve()
            self.assertEqual(fit_done, case[6], "fit check failed for %d displacement" % (case[4]))
            if (case[6]):
                self.check_bounds_continuous(curve, constraint)
                self.check_bounds_discrete(curve, constraint)
                self.check_discretize(curve)