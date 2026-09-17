import math
import unittest

import numpy as np

from pyThermoModels.activity import (
    ActivityCore,
    Margules,
    RedlichKister,
    VanLaar,
    Wilson,
)


class Phase3ClassicalActivityModelTests(unittest.TestCase):
    def assert_close(self, first, second, rel=1e-9, abs_tol=1e-12):
        self.assertTrue(
            math.isclose(first, second, rel_tol=rel, abs_tol=abs_tol),
            f"{first} != {second}",
        )

    def assert_gibbs_identity(self, x, gamma, gE_RT):
        identity = float(np.sum(np.asarray(x) * np.log(np.asarray(gamma))))
        self.assert_close(identity, gE_RT)

    def test_wilson_ideal_limit_and_gibbs_identity(self):
        model = Wilson(components=["A", "B", "C"])
        model_input = {
            "mole_fraction": {"A": 0.2, "B": 0.3, "C": 0.5},
            "lambda_ij": np.ones((3, 3)),
        }

        res, other = model.cal(model_input)

        np.testing.assert_allclose(res["value"], [1.0, 1.0, 1.0])
        self.assert_close(other["excess_gibbs_RT"], 0.0)
        self.assert_gibbs_identity(
            res["mole_fraction"],
            res["value"],
            other["excess_gibbs_RT"],
        )

    def test_wilson_binary_hand_calculation(self):
        model = Wilson(components=["A", "B"])
        lambda_ij = np.array([[1.0, 0.8], [1.4, 1.0]])
        x = np.array([0.3, 0.7])

        res, other = model.cal({
            "mole_fraction": {"A": x[0], "B": x[1]},
            "lambda_ij": lambda_ij,
        })

        sums = lambda_ij @ x
        expected_ln = np.array([
            1.0 - math.log(sums[0]) - np.sum(x * lambda_ij[:, 0] / sums),
            1.0 - math.log(sums[1]) - np.sum(x * lambda_ij[:, 1] / sums),
        ])
        expected_gamma = np.exp(expected_ln)
        expected_gE_RT = float(-np.sum(x * np.log(sums)))

        np.testing.assert_allclose(res["value"], expected_gamma)
        self.assert_close(other["excess_gibbs_RT"], expected_gE_RT)
        self.assert_gibbs_identity(x, res["value"], expected_gE_RT)

    def test_margules_two_suffix_hand_calculation(self):
        model = Margules(components=["A", "B"])
        res, other = model.cal(
            {
                "mole_fraction": {"A": 0.4, "B": 0.6},
                "A": 1.2,
            },
            calculation_mode="two_suffix",
        )

        expected_ln = np.array([1.2 * 0.6**2, 1.2 * 0.4**2])
        expected_gamma = np.exp(expected_ln)

        np.testing.assert_allclose(res["value"], expected_gamma)
        self.assert_close(other["excess_gibbs_RT"], 1.2 * 0.4 * 0.6)
        self.assert_gibbs_identity(
            res["mole_fraction"],
            res["value"],
            other["excess_gibbs_RT"],
        )

    def test_margules_three_suffix_hand_calculation(self):
        model = Margules(components=["A", "B"])
        x1 = 0.25
        x2 = 0.75
        A12 = 0.9
        A21 = 1.4

        res, other = model.cal({
            "mole_fraction": {"A": x1, "B": x2},
            "A12": A12,
            "A21": A21,
        })

        g = x1 * x2 * (A12 * x1 + A21 * x2)
        dg = (1.0 - 2.0 * x1) * (A12 * x1 + A21 * x2) + x1 * x2 * (A12 - A21)
        expected_ln = np.array([g + x2 * dg, g - x1 * dg])

        np.testing.assert_allclose(res["value"], np.exp(expected_ln))
        self.assert_close(other["excess_gibbs_RT"], g)
        self.assert_gibbs_identity([x1, x2], res["value"], g)

    def test_van_laar_hand_calculation(self):
        model = VanLaar(components=["A", "B"])
        x1 = 0.35
        x2 = 0.65
        a1 = 1.3
        a2 = 0.7

        res, other = model.cal({
            "mole_fraction": {"A": x1, "B": x2},
            "a1": a1,
            "a2": a2,
        })

        expected_gE_RT = a1 * a2 * x1 * x2 / (a1 * x1 + a2 * x2)

        self.assertTrue(all(gamma > 0 for gamma in res["value"]))
        self.assert_close(other["excess_gibbs_RT"], expected_gE_RT)
        self.assert_gibbs_identity([x1, x2], res["value"], expected_gE_RT)

    def test_redlich_kister_hand_calculation(self):
        model = RedlichKister(components=["A", "B"])
        x1 = 0.4
        x2 = 0.6
        coeffs = [1.0, 0.2, -0.1]

        res, other = model.cal({
            "mole_fraction": {"A": x1, "B": x2},
            "a": coeffs,
        })

        d = x1 - x2
        series = coeffs[0] + coeffs[1] * d + coeffs[2] * d**2
        dseries = coeffs[1] * 2.0 + coeffs[2] * 2.0 * d * 2.0
        expected_gE_RT = x1 * x2 * series
        dg = (1.0 - 2.0 * x1) * series + x1 * x2 * dseries
        expected_ln = np.array([
            expected_gE_RT + x2 * dg,
            expected_gE_RT - x1 * dg,
        ])

        np.testing.assert_allclose(res["value"], np.exp(expected_ln))
        self.assert_close(other["excess_gibbs_RT"], expected_gE_RT)
        self.assert_gibbs_identity([x1, x2], res["value"], expected_gE_RT)

    def test_validation_guards_and_activity_core_selection(self):
        with self.assertRaises(Exception):
            Margules(components=["A", "B", "C"]).cal({
                "mole_fraction": {"A": 0.2, "B": 0.3, "C": 0.5},
                "A": 1.0,
            })

        with self.assertRaises(Exception):
            VanLaar(components=["A", "B"]).cal({
                "mole_fraction": {"A": 0.4, "B": 0.7},
                "a1": 1.0,
                "a2": 1.0,
            })

        core = ActivityCore(datasource={}, equationsource={}, components=["A", "B"])
        self.assertIsInstance(core.select("WILSON"), Wilson)
        self.assertIsInstance(core.select("MARGULES"), Margules)
        self.assertIsInstance(core.select("VAN_LAAR"), VanLaar)
        self.assertIsInstance(core.select("REDLICH_KISTER"), RedlichKister)


if __name__ == "__main__":
    unittest.main()
