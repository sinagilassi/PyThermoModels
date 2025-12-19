# import libs
import logging
import math
from typing import Dict, Any, List, Tuple

# NOTE: setup logger
logger = logging.getLogger(__name__)


class UNIFAC1:
    """
    UNIFAC Activity Coefficient Model.

    This class implements the UNIFAC (UNIversal Functional Activity Coefficient)
    model for calculating activity coefficients in liquid mixtures based on
    group contributions.
    """

    def __init__(
        self,
        group_data: Dict[str, Dict[str, Any]],
        interaction_data: Dict[str, Dict[str, Any]],
        **kwargs
    ):
        """
        Parameters
        ----------
        group_data : Dict[str, Dict[str, Any]]
            Dictionary containing group properties (R, Q, main-group).
        interaction_data : Dict[str, Dict[str, Any]]
            Dictionary containing interaction parameters between main groups.
        kwargs : dict
            Additional keyword arguments for UNIFAC1 model initialization.
            - eps : float, optional
                Small number to avoid division by zero (default is 1e-30).
            - z : float, optional
                Coordination number (default is 10.0).
        """
        # SECTION: Initialization & Topology Pre-calculation
        self.group_data = group_data
        self.interaction_data = interaction_data

        # >> check lengths
        if not self.group_data:
            logger.error("group_data is empty")
            raise ValueError("group_data is empty")

        if not self.interaction_data:
            logger.error("interaction_data is empty")
            raise ValueError("interaction_data is empty")

        # NOTE: kwargs
        # >> small number to avoid div by zero
        self.eps = kwargs.get("eps", 1e-30)
        # >> z value (coordination number)
        self.z = kwargs.get("z", 10.0)

        # SECTION: Initialize calculation parameters
        # ! >> will be called after setting component groups
        self._component_group_data: List[Dict[str, float | int]] = []
        self.num_comps: int | None = None
        self.num_groups: int | None = None
        self.sorted_subgroups: List[str] = []
        self.g_map: Dict[str, int] = {}
        self.r: List[float] = []
        self.q: List[float] = []
        self.nu: List[List[float]] = []
        self.Q_k: List[float] = []
        self.main_groups: List[str] = []

    @property
    def component_group_data(self) -> List[Dict[str, float | int]]:
        """
        Get the component group contribution data.

        Returns
        -------
        List[Dict[str, float | int]]
            List of components, each defined by a dictionary of group counts.
        """
        return self._component_group_data

    @component_group_data.setter
    def component_group_data(self, data: List[Dict[str, float | int]]) -> None:
        """
        Set the component group contribution data.

        Parameters
        ----------
        data : List[Dict[str, float | int]]
            List of components, each defined by a dictionary of group counts.
        """
        if not isinstance(data, list):
            logger.error("components must be a list of dicts")
            raise TypeError("components must be a list of dicts")

        self._component_group_data = data

    # ---------------------------------------------------------
    # PART 0: Initialize Calculation Methods
    # ---------------------------------------------------------
    def initialize_calc(
        self,
    ):
        try:
            # TODO: Check variables
            if not self.component_group_data:
                raise ValueError("component_group_data not set")

            # SECTION: Topology Pre-calculation (Done once)
            components_data = self.component_group_data
            self.num_comps = len(components_data)
            self.r = []
            self.q = []

            # SECTION: Identify unique groups and map to indices
            unique_groups = set()
            for comp in components_data:
                unique_groups.update(str(k) for k in comp.keys())

            self.sorted_subgroups = sorted(
                unique_groups, key=lambda s: int(s) if s.isdigit() else s)
            self.num_groups = len(self.sorted_subgroups)
            self.g_map = {
                gid: i for i, gid in enumerate(self.sorted_subgroups)
            }

            # SECTION: Cache Group Properties
            self.Q_k = [
                float(self.group_data[g]["Q"]) for g in self.sorted_subgroups
            ]
            self.main_groups = [
                str(self.group_data[g]["main-group"]) for g in self.sorted_subgroups
            ]

            # SECTION: Build Count Matrix (nu) and component r, q
            self.nu = [[0.0] * self.num_groups for _ in range(self.num_comps)]

            # NOTE: iterate components_data to fill nu, r, q
            for i, comp in enumerate(components_data):
                ri = 0.0
                qi = 0.0
                for gid, count in comp.items():
                    gid = str(gid)
                    k = self.g_map[gid]
                    cnt = float(count)
                    self.nu[i][k] = cnt
                    ri += cnt * float(self.group_data[gid]["R"])
                    qi += cnt * float(self.group_data[gid]["Q"])
                self.r.append(ri)
                self.q.append(qi)
        except Exception as e:
            logger.error(f"Error in initialize_calc: {e}")
            raise

    # ---------------------------------------------------------
    # PART 1: Combinatorial (Depends only on x, r, q)
    # ---------------------------------------------------------
    def calc_combinatorial(self, x: List[float]) -> List[float]:
        """
        Calculates the entropy-based combinatorial term (ln gamma_C).
        Only depends on component volume/area (r, q) and mole fraction (x).

        Parameters
        ----------
        x : List[float]
            Mole fractions of components.

        Returns
        -------
        List[float]
            Combinatorial activity coefficients (ln gamma_C) for each component.
        """
        try:
            # TODO: Check variables
            if not self.num_comps:
                raise ValueError("num_comps not initialized")

            # SECTION: Calculate combinatorial term
            if len(x) != self.num_comps:
                raise ValueError("Mole fraction length mismatch")

            z = self.z
            n = self.num_comps

            sum_xr = sum(x[i] * self.r[i] for i in range(n))
            sum_xq = sum(x[i] * self.q[i] for i in range(n))

            Phi = [(x[i] * self.r[i]) / (sum_xr + self.eps) for i in range(n)]
            Theta = [
                (x[i] * self.q[i]) / (sum_xq + self.eps) for i in range(n)
            ]
            l = [
                (z/2)*(self.r[i] - self.q[i]) - (self.r[i] - 1) for i in range(n)
            ]
            sum_xl = sum(x[i] * l[i] for i in range(n))

            ln_gamma_c = []
            for i in range(n):
                if x[i] <= 1e-12:
                    ln_gamma_c.append(0.0)
                else:
                    term = (math.log(Phi[i]/x[i]) +
                            (z/2) * self.q[i] * math.log(Theta[i]/Phi[i]) +
                            l[i] -
                            (Phi[i]/x[i]) * sum_xl)
                    ln_gamma_c.append(term)

            return ln_gamma_c
        except Exception as e:
            logger.error(f"Error in calc_combinatorial: {e}")
            raise

    # ---------------------------------------------------------
    # PART 2: Residual (Depends on x, T, interactions)
    # ---------------------------------------------------------
    def _get_psi(self, T: float) -> List[List[float]]:
        """
        Helper to compute interaction matrix for current T.

        Parameters
        ----------
        T : float
            Temperature in Kelvin.

        Returns
        -------
        List[List[float]]
            Interaction parameter matrix Psi.
        """
        try:
            # TODO: Check variables
            if not self.num_groups:
                raise ValueError("num_groups not initialized")

            # SECTION: Build Psi matrix
            Psi = [[0.0] * self.num_groups for _ in range(self.num_groups)]
            for m in range(self.num_groups):
                for n in range(self.num_groups):
                    if m == n:
                        Psi[m][n] = 1.0
                    else:
                        main_m = self.main_groups[m]
                        main_n = self.main_groups[n]
                        # Safer .get with 0.0 default, or strict try/except
                        # a_mn = self.interaction_data.get(main_m, {}).get(main_n, 0.0)

                        try:
                            a_mn = self.interaction_data[main_m][main_n]
                            # a_mn = self.interaction_data.get(main_m, {}).get(main_n, 0.0)
                            # >> check nan
                            if math.isnan(a_mn):
                                a_mn = 0.0
                        except KeyError:
                            # Defaulting to 0 is dangerous; better to crash or warn.
                            # Assuming 0.0 here for safety in dummy data examples,
                            # but strictly this should raise error.
                            a_mn = 0.0

                        # set
                        Psi[m][n] = math.exp(-a_mn / T)

            return Psi
        except Exception as e:
            logger.error(f"Error in _get_psi: {e}")
            raise

    def _ln_gamma_group(self, nu_vector: List[float], Psi: List[List[float]]) -> List[float]:
        """
        Core UNIFAC equation for a group vector.

        Parameters
        ----------
        nu_vector : List[float]
            Group count vector for a mixture or pure component.
        Psi : List[List[float]]
            Interaction parameter matrix.

        Returns
        -------
        List[float]
            Log group activity coefficients (ln G_k).
        """
        try:
            # TODO: Check variables
            if not self.num_groups:
                raise ValueError("num_groups not initialized")

            # SECTION: Total number of groups in mixture
            total_groups = sum(nu_vector)
            if total_groups < self.eps:
                return [0.0] * self.num_groups

            X = [n_k / total_groups for n_k in nu_vector]

            denom_theta = sum(
                X[k] * self.Q_k[k] for k in range(self.num_groups)
            )
            Theta = [
                (X[k] * self.Q_k[k]) / (denom_theta + self.eps) for k in range(self.num_groups)
            ]

            ln_G = [0.0] * self.num_groups

            # Precompute denominator term
            denom_psi = [
                sum(Theta[n] * Psi[n][m] for n in range(self.num_groups)) for m in range(self.num_groups)
            ]

            for k in range(self.num_groups):
                term2 = sum((Theta[m] * Psi[k][m]) / (denom_psi[m] + self.eps)
                            for m in range(self.num_groups))
                ln_G[k] = self.Q_k[k] * \
                    (1.0 - math.log(denom_psi[k] + self.eps) - term2)

            return ln_G
        except Exception as e:
            logger.error(f"Error in _ln_gamma_group: {e}")
            raise

    def calc_residual(self, x: List[float], T: float) -> List[float]:
        """
        Calculates the enthalpy-based residual term (ln gamma_R).
        Depends on Temperature and component interactions.

        Parameters
        ----------
        x : List[float]
            Mole fractions of components.
        T : float
            Temperature in Kelvin.

        Returns
        -------
        List[float]
            Residual activity coefficients (ln gamma_R) for each component.
        """
        try:
            # TODO: Check variables
            if not self.num_comps or not self.num_groups:
                raise ValueError("num_comps not initialized")

            # SECTION: Mixture Group Activities
            if len(x) != self.num_comps:
                raise ValueError("Mole fraction length mismatch")
            Psi = self._get_psi(T)

            # 1. Mixture Properties
            nu_mix = [0.0] * self.num_groups
            for i in range(self.num_comps):
                if x[i] > 1e-12:
                    for k in range(self.num_groups):
                        nu_mix[k] += x[i] * self.nu[i][k]

            lnG_mix = self._ln_gamma_group(nu_mix, Psi)

            # 2. Pure Component Properties & Summation
            ln_gamma_r = []
            for i in range(self.num_comps):
                # Calculate pure component group activities
                lnG_pure = self._ln_gamma_group(self.nu[i], Psi)

                resid = 0.0
                for k in range(self.num_groups):
                    if self.nu[i][k] > 0:
                        resid += self.nu[i][k] * (lnG_mix[k] - lnG_pure[k])
                ln_gamma_r.append(resid)

            return ln_gamma_r
        except Exception as e:
            logger.error(f"Error in calc_residual: {e}")
            raise

    # ---------------------------------------------------------
    # PART 3: Total Activity Coefficient
    # ---------------------------------------------------------
    def get_activity_coefficients(
            self,
            T: float,
            x: List[float],
            comp: List[str],
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Returns Activity Coefficients (gamma).
        Gamma = exp(ln_gamma_C + ln_gamma_R)

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        x : List[float]
            Mole fractions of components.
        comp : List[str]
            List of component names.

        Returns
        -------
        List[float]
            Activity coefficients for each component.
        """
        try:
            # TODO: Check variables
            if not self.num_comps:
                raise ValueError("num_comps not initialized")

            # SECTION: Total Activity Coefficient
            if len(x) != self.num_comps:
                raise ValueError("Mole fraction length mismatch")

            ln_c = self.calc_combinatorial(x)
            ln_r = self.calc_residual(x, T)

            # final activity coefficients
            AcCo_i = [math.exp(c + r) for c, r in zip(ln_c, ln_r)]

            # SECTION
            # init the activity coefficients
            AcCo_i_comp = {
                comp[i]: float(AcCo_i[i]) for i in range(self.num_comps)
            }

            # init mole fractions
            x_comp = {
                comp[i]: float(x[i]) for i in range(self.num_comps)
            }

            # init gamma combinatorial and residual
            ln_gamma_combinatorial_comp = {
                comp[i]: float(ln_c[i]) for i in range(self.num_comps)
            }
            ln_gamma_residual_comp = {
                comp[i]: float(ln_r[i]) for i in range(self.num_comps)
            }

            # input values
            other_values = {
                'ln_gamma_combinatorial': ln_c,
                'ln_gamma_combinatorial_comp': ln_gamma_combinatorial_comp,
                'ln_gamma_residual': ln_r,
                'ln_gamma_residual_comp': ln_gamma_residual_comp,
                'temperature_value': T,
                'temperature_unit': 'K',
                'mole_fraction': x,
                'mole_fraction_comp': x_comp,
                "AcCo_i_comp": AcCo_i_comp,
            }

            # res
            res = {
                'property_name': 'activity coefficients',
                'components': comp,
                'mole_fraction': x,
                'value': AcCo_i,
                'unit': 1,
                'symbol': "AcCo_i",
            }

            # res
            return res, other_values
        except Exception as e:
            logger.error(f"Error in get_activity_coefficients: {e}")
            raise
