import math
import unittest

from pyThermoModels.configs import R_CONST
from pyThermoModels.eos.alpha import (
    alpha_value,
    attraction_parameter_derivatives,
    dalpha_dT,
    d2alpha_dT2,
)
from pyThermoModels.eos.departure import residual_properties_cubic_pure
from pyThermoModels.eos.eosmanager import EOSManager
from pyThermoModels.eos.eosmodels import EOSModels


class Phase1CubicTests(unittest.TestCase):
    def setUp(self):
        self.datasource = {
            "methane": {
                "Pc": {"value": 4.5992e6, "unit": "Pa"},
                "Tc": {"value": 190.564, "unit": "K"},
            }
        }
        self.eos = EOSModels(self.datasource, {})
        self.manager = EOSManager(self.datasource, {})

    def finite_difference(self, func, x_value, step):
        return (func(x_value + step) - func(x_value - step)) / (2 * step)

    def second_finite_difference(self, func, x_value, step):
        return (func(x_value + step) - 2 * func(x_value) + func(x_value - step)) / (step ** 2)

    def assert_close(self, first, second, rel=1e-9, abs_tol=1e-12):
        self.assertTrue(
            math.isclose(first, second, rel_tol=rel, abs_tol=abs_tol),
            f"{first} != {second}",
        )

    def test_alpha_values_match_legacy_formulas(self):
        T = 300.0
        Tc = 500.0
        Tr = T / Tc
        for omega in [0.0, 0.07780, 0.2]:
            pr_legacy = (
                1
                + (0.37464 + 1.54226 * omega - 0.26992 * omega ** 2)
                * (1 - Tr ** 0.5)
            ) ** 2
            srk_legacy = (
                1
                + (0.480 + 1.574 * omega - 0.176 * omega ** 2)
                * (1 - Tr ** 0.5)
            ) ** 2
            self.assert_close(alpha_value(T, Tc, omega, model="PR"), pr_legacy)
            self.assert_close(alpha_value(T, Tc, omega, model="SRK"), srk_legacy)

        self.assert_close(alpha_value(T, Tc, model="RK"), Tr ** -0.5)
        self.assert_close(alpha_value(T, Tc, model="vdW"), 1.0)

    def test_alpha_derivatives_match_finite_difference(self):
        for model, omega in [("PR", 0.07780), ("SRK", 0.08664), ("RK", None)]:
            T = 320.0
            Tc = 500.0
            step = 1e-3
            func = lambda temp: alpha_value(temp, Tc, omega, model=model)
            self.assert_close(
                dalpha_dT(T, Tc, omega, model=model),
                self.finite_difference(func, T, step),
                rel=2e-7,
            )
            self.assert_close(
                d2alpha_dT2(T, Tc, omega, model=model),
                self.second_finite_difference(func, T, step),
                rel=2e-4,
                abs_tol=1e-9,
            )

    def test_alpha_rejects_invalid_inputs(self):
        with self.assertRaises(ValueError):
            alpha_value(0.0, 500.0, model="PR")
        with self.assertRaises(ValueError):
            alpha_value(300.0, 0.0, model="PR")
        with self.assertRaises(ValueError):
            alpha_value(300.0, 500.0, model="Twu")

    def test_eos_parameters_preserve_legacy_alpha_outputs(self):
        P = 1.0e6
        T = 300.0
        Tc = self.datasource["methane"]["Tc"]["value"]
        for method in ["vdW", "RK", "SRK", "PR"]:
            params = self.eos.eos_parameters(P, T, "methane", method=method)
            omega = params["omega"]
            Tr = T / Tc
            if method == "vdW":
                expected_alpha = 1.0
            elif method == "RK":
                expected_alpha = Tr ** -0.5
            elif method == "SRK":
                expected_alpha = (
                    1
                    + (0.480 + 1.574 * omega - 0.176 * omega ** 2)
                    * (1 - Tr ** 0.5)
                ) ** 2
            else:
                expected_alpha = (
                    1
                    + (0.37464 + 1.54226 * omega - 0.26992 * omega ** 2)
                    * (1 - Tr ** 0.5)
                ) ** 2
            self.assert_close(params["alpha"], expected_alpha)

    def residual_case(self, model, phase, P=1.0e6, T=300.0):
        root_id = {"LIQUID": 2, "VAPOR": 3, "SUPERCRITICAL": 4}[phase]
        roots, params_list, _ = self.manager.eos_roots(
            P,
            T,
            ["methane"],
            {"root": [root_id]},
            eos_model=model,
            solver_method="root",
            mode="single",
        )
        params = params_list[0]
        Tc = self.datasource["methane"]["Tc"]["value"]
        Pc = self.datasource["methane"]["Pc"]["value"]


        return residual_properties_cubic_pure(
            P=P,
            T=T,
            Z=float(roots[0]),
            Tc=Tc,
            Pc=Pc,
            alpha_acentric_input=params["omega"],
            eos_model=model,
            sigma=params["sigma"],
            epsilon=params["epsilon"],
            psi=params["psi"],
            omega=params["omega"],
            phase=phase,
        )

    def test_residual_property_identities_for_pr_and_srk(self):
        for model in ["PR", "SRK"]:
            residual = self.residual_case(model, "VAPOR")
            self.assert_close(
                residual.residual_gibbs,
                residual.residual_enthalpy - 300.0 * residual.residual_entropy,
                rel=1e-10,
            )
            self.assert_close(
                residual.residual_internal_energy,
                residual.residual_enthalpy
                - R_CONST * 300.0 * (residual.compressibility_factor - 1),
                rel=1e-10,
            )
            self.assertTrue(math.isfinite(residual.residual_isobaric_heat_capacity))
            self.assertTrue(math.isfinite(residual.residual_isochoric_heat_capacity))

    def test_residual_cp_matches_finite_difference_of_enthalpy(self):
        for model in ["PR", "SRK"]:
            T = 300.0
            P = 1.0e6
            step = 1e-3
            residual = self.residual_case(model, "VAPOR", P=P, T=T)
            h_func = lambda temp: self.residual_case(model, "VAPOR", P=P, T=temp).residual_enthalpy
            finite_cp = self.finite_difference(h_func, T, step)
            self.assert_close(
                residual.residual_isobaric_heat_capacity,
                finite_cp,
                rel=2e-5,
                abs_tol=1e-5,
            )
    def test_residual_properties_approach_zero_at_low_pressure(self):
        residual = self.residual_case("PR", "VAPOR", P=10.0, T=300.0)
        self.assertLess(abs(residual.residual_enthalpy), 1e-1)
        self.assertLess(abs(residual.residual_entropy), 1e-3)
        self.assertLess(abs(residual.residual_gibbs), 1e-1)

    def test_attraction_parameter_derivatives_are_scaled_alpha_derivatives(self):
        ac = 2.5
        T = 350.0
        Tc = 500.0
        da, d2a = attraction_parameter_derivatives(ac, T, Tc, 0.1, model="PR")
        self.assert_close(da, ac * dalpha_dT(T, Tc, 0.1, model="PR"))
        self.assert_close(d2a, ac * d2alpha_dT2(T, Tc, 0.1, model="PR"))


if __name__ == "__main__":
    unittest.main()


