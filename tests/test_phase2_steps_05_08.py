import math
import unittest

import numpy as np

from pyThermoModels.eos.eosmodels import EOSModels
from pyThermoModels.eos.lee_kesler import lee_kesler_ploecker_mixture, lee_kesler_pure
from pyThermoModels.eos.mixing import classical_quadratic_mixing, wong_sandler_pr_nrtl_mixing


class Phase2LeeKeslerAndMixingTests(unittest.TestCase):
    def setUp(self):
        self.params = [
            {"a": 1.5, "b": 2.0e-5, "A": 0.12, "B": 0.01, "eos-model": "PR"},
            {"a": 0.8, "b": 3.0e-5, "A": 0.08, "B": 0.02, "eos-model": "PR"},
        ]
        self.x = np.array([0.35, 0.65])
        self.kij = np.array([[0.0, 0.07], [0.07, 0.0]])

    def assert_close(self, first, second, rel=1e-9, abs_tol=1e-12):
        self.assertTrue(
            math.isclose(first, second, rel_tol=rel, abs_tol=abs_tol),
            f"{first} != {second}",
        )

    def test_classical_mixing_matches_legacy_formula(self):
        result = classical_quadratic_mixing(self.x, self.params, self.kij)
        expected_aij = np.array(
            [
                [math.sqrt(1.5 * 1.5), (1 - 0.07) * math.sqrt(1.5 * 0.8)],
                [(1 - 0.07) * math.sqrt(0.8 * 1.5), math.sqrt(0.8 * 0.8)],
            ]
        )
        expected_amix = float(np.sum(self.x[:, None] * self.x[None, :] * expected_aij))
        expected_bmix = float(np.dot(self.x, [2.0e-5, 3.0e-5]))
        self.assert_close(result.a_mix, expected_amix)
        self.assert_close(result.b_mix, expected_bmix)
        np.testing.assert_allclose(result.a_ij, expected_aij)

        legacy = EOSModels({}, {}).eos_mixing_rule(self.x, self.params, self.kij)
        self.assert_close(legacy[0], result.a_mix)
        self.assert_close(legacy[1], result.b_mix)
        np.testing.assert_allclose(legacy[2], result.a_ij)

    def test_lee_kesler_low_pressure_limit(self):
        result = lee_kesler_pure(
            P=100.0,
            T=300.0,
            Tc=190.564,
            Pc=4.5992e6,
            acentric_factor=0.011,
        )
        self.assertLess(abs(result.compressibility_factor - 1.0), 1e-3)
        self.assertLess(abs(result.enthalpy_departure), 5.0)
        self.assertLess(abs(result.entropy_departure), 0.05)
        self.assertFalse(result.metadata["fugacity_coefficient_supported"])

    def test_lee_kesler_simple_and_reference_limits(self):
        simple = lee_kesler_pure(1.0e6, 350.0, 500.0, 5.0e6, 0.0)
        reference = lee_kesler_pure(1.0e6, 350.0, 500.0, 5.0e6, 0.3978)
        middle = lee_kesler_pure(1.0e6, 350.0, 500.0, 5.0e6, 0.1989)
        self.assert_close(
            middle.compressibility_factor,
            0.5 * (simple.compressibility_factor + reference.compressibility_factor),
            rel=1e-7,
        )

    def test_ploecker_pure_component_limit_and_metadata(self):
        pure = lee_kesler_pure(
            P=1.0e5,
            T=300.0,
            Tc=190.564,
            Pc=4.5992e6,
            acentric_factor=0.011,
        )
        mix = lee_kesler_ploecker_mixture(
            P=1.0e5,
            T=300.0,
            components=["methane"],
            x=[1.0],
            Tc=[190.564],
            Pc=[4.5992e6],
            acentric_factors=[0.011],
            Zc=[0.286],
        )
        self.assert_close(mix.result.compressibility_factor, pure.compressibility_factor, rel=2e-3)
        self.assertTrue(mix.metadata["k_ij_defaulted_to_zero"])

    def test_ploecker_permutation_invariance(self):
        first = lee_kesler_ploecker_mixture(
            P=1.0e5,
            T=330.0,
            components=["methane", "ethane"],
            x=[0.25, 0.75],
            Tc=[190.564, 305.32],
            Pc=[4.5992e6, 4.872e6],
            acentric_factors=[0.011, 0.099],
            Zc=[0.286, 0.279],
        )
        second = lee_kesler_ploecker_mixture(
            P=1.0e5,
            T=330.0,
            components=["ethane", "methane"],
            x=[0.75, 0.25],
            Tc=[305.32, 190.564],
            Pc=[4.872e6, 4.5992e6],
            acentric_factors=[0.099, 0.011],
            Zc=[0.279, 0.286],
        )
        self.assert_close(first.result.compressibility_factor, second.result.compressibility_factor)
        self.assert_close(first.pseudo_critical_temperature, second.pseudo_critical_temperature)

    def test_wong_sandler_pure_limit_and_scope_guard(self):
        pure_result = wong_sandler_pr_nrtl_mixing(
            T=300.0,
            xi=[1.0],
            params_list=[self.params[0]],
            excess_gibbs_over_rt=lambda _T, _x: 0.0,
        )
        self.assert_close(pure_result.b_mix, self.params[0]["b"], rel=1e-7)
        self.assert_close(pure_result.a_mix, self.params[0]["a"], rel=1e-7)
        self.assertEqual(pure_result.metadata["rule"], "wong_sandler")

        with self.assertRaises(ValueError):
            wong_sandler_pr_nrtl_mixing(
                T=300.0,
                xi=[1.0],
                params_list=[{**self.params[0], "eos-model": "SRK"}],
                excess_gibbs_over_rt=lambda _T, _x: 0.0,
            )

    def test_wong_sandler_responds_to_excess_gibbs_callback(self):
        ideal = wong_sandler_pr_nrtl_mixing(
            T=300.0,
            xi=self.x,
            params_list=self.params,
            excess_gibbs_over_rt=lambda _T, _x: 0.0,
            k_ij=self.kij,
        )
        nonideal = wong_sandler_pr_nrtl_mixing(
            T=300.0,
            xi=self.x,
            params_list=self.params,
            excess_gibbs_over_rt=lambda _T, _x: 0.25,
            k_ij=self.kij,
        )
        self.assertNotEqual(ideal.a_mix, nonideal.a_mix)
        self.assertNotEqual(ideal.b_mix, nonideal.b_mix)


if __name__ == "__main__":
    unittest.main()
