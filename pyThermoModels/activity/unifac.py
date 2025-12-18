# import libs
import logging
from typing import Dict, Any, List, Tuple
import math

# NOTE: logger
logger = logging.getLogger(__name__)


class UNIFAC:
    """
    Classic UNIFAC (Fredenslund–Jones–Prausnitz) using:
    - subgroup parameters: R_k, Q_k, and subgroup -> main-group mapping
    - main-group interaction parameters a_mn (K) stored as nested dict

    Data formats
    ------------
    group_data:
        {
        "1": {"main-group": "1", "R": 0.9011, "Q": 0.8480},
        "2": {"main-group": "1", "R": 0.6744, "Q": 0.5400},
        ...
        }

    interaction_data (classic UNIFAC):
        interaction_data[main_m][main_n] = a_mn   # units: K
        Example:
        {
        "1": {"1": 0.0, "2": 200.0},
        "2": {"1": 100.0, "2": 0.0}
        }

    Notes
    -----
    - This implements *classic UNIFAC* equations.
    - If modified=True, it raises NotImplementedError (Dortmund needs different parameter sets & formulas).
    """

    def __init__(
        self,
        group_data: Dict[str, Dict[str, Any]],
        interaction_data: Dict[str, Dict[str, Any]],
        modified: bool = False,
        **kwargs
    ):
        """
        Initialize UNIFAC model.

        Parameters
        ----------
        group_data : dict
            Subgroup data with R, Q, and main-group mapping.
        interaction_data : dict
            Main-group interaction parameters a_mn (K).
        modified : bool, optional
            If True, use modified UNIFAC (Dortmund). Default is False.
        kwargs : dict
            Additional parameters (e.g., eps for numerical stability).
        """
        # NOTE: attributes
        self.group_data = group_data
        self.interaction_data = interaction_data
        self.modified = modified

        # NOTE: kwargs
        self.eps = kwargs.get("eps", 1e-30)

        # NOTE: Map subgroup -> main-group
        self.sub_to_main = {
            str(sub): str(info["main-group"]) for sub, info in group_data.items()
        }

        # NOTE: Cache subgroup R, Q as floats
        self.sub_R = {
            str(sub): float(info["R"]) for sub, info in group_data.items()
        }
        self.sub_Q = {
            str(sub): float(info["Q"]) for sub, info in group_data.items()
        }

        if self.modified:
            logger.error("Modified UNIFAC (Dortmund) not implemented.")
            raise NotImplementedError(
                "Modified UNIFAC (Dortmund) not implemented."
            )

    def calculate_component_RQ(
        self,
        comp_groups: Dict[str, float]
    ) -> Tuple[float, float]:
        """
        Calculate component R and Q from subgroup counts.

        Parameters
        ----------
        comp_groups : dict
            Dictionary of subgroup counts for the component.

        Returns
        -------
        R : float
            Total R value for the component.
        Q : float
            Total Q value for the component.
        """
        try:
            # NOTE: initialize R and Q
            R = 0.0
            Q = 0.0

            # NOTE: iterate over subgroup counts
            for gid, count in comp_groups.items():
                gid = str(gid)
                if gid not in self.group_data:
                    raise KeyError(f"Unknown subgroup id {gid}")
                R += self.sub_R[gid] * float(count)
                Q += self.sub_Q[gid] * float(count)

            # res
            return R, Q
        except KeyError as e:
            logger.error(f"Error in calculate_component_RQ: {e}")
            raise

    def combinatorial_term(
            self,
            x: List[float],
            r: List[float],
            q: List[float],
            z: float = 10,
            modified: bool = False,
            eps: float = 1e-30
    ) -> List[float]:
        """
        Uses the "Phi, theta" form of the combinatorial term.:

        Phi_i   = (r_i x_i) / Σ(r x)
        theta_i = (q_i x_i) / Σ(q x)

        ln(gamma_i^C) =
        1 - Phi_i/x_i + ln(Phi_i/x_i)
        - 5 q_i [ 1 - Phi_i/theta_i + ln(Phi_i/theta_i) ]

        z and modified are unused for this classic form (kept for signature compatibility).

        Parameters
        ----------
        x : list of float
            Mole fractions of components.
        r : list of float
            r parameters of components.
        q : list of float
            q parameters of components.
        z : float, optional
            Coordination number (unused in classic UNIFAC).
        modified : bool, optional
            If True, use modified UNIFAC (unused in classic UNIFAC).
        eps : float, optional
            Small value to prevent division by zero.
        """
        try:
            n = len(x)
            sx = sum(x)
            if abs(sx - 1.0) > 1e-8:
                raise ValueError(
                    f"x must be normalized (sum=1). Got sum(x)={sx}")

            sum_xr = sum(x[i] * r[i] for i in range(n))
            sum_xq = sum(x[i] * q[i] for i in range(n))

            Phi = [(x[i] * r[i]) / (sum_xr + eps) for i in range(n)]
            theta = [(x[i] * q[i]) / (sum_xq + eps) for i in range(n)]

            # l_i
            l = [(z/2.0) * (r[i] - q[i]) - (r[i] - 1.0) for i in range(n)]
            sum_xl = sum(x[i] * l[i] for i in range(n))

            # optional modified phi' (ONLY for the ln(phi/x) term)
            if modified:
                sum_xr075 = sum(x[i] * (r[i] ** 0.75) for i in range(n))
                phi_prime = [
                    (r[i] ** 0.75) / (sum_xr075 + eps) for i in range(n)
                ]
            else:
                phi_prime = Phi

            ln_gamma_c = [0.0] * n
            for i in range(n):
                if x[i] <= 0.0:
                    ln_gamma_c[i] = 0.0
                    continue

                term1 = math.log((phi_prime[i] + eps) / (x[i] + eps))
                term2 = (z/2.0) * q[i] * \
                    math.log((theta[i] + eps) / (Phi[i] + eps))
                term3 = l[i] - (Phi[i] / (x[i] + eps)) * sum_xl

                ln_gamma_c[i] = term1 + term2 + term3

            return ln_gamma_c
        except Exception as e:
            logger.error(f"Error in combinatorial_term: {e}")
            raise

    def residual_term(
            self,
            T: float,
            mole_fractions: List[float],
            composition: List[Dict[str, float]]
    ) -> List[float]:
        """
        Calculate residual part of activity coefficients using UNIFAC model.

        Steps:
        - Build mixture subgroup totals: N_k = Σ_i x_i ν_ik
        - Compute mixture X_k and Θ_k over *subgroups*
        - Use Ψ_mn built from *main-group* a_mn, mapped from subgroup->main
        - Compute lnΓ_k for mixture and for each pure component i
        - lnγ_i^R = Σ_k ν_ik (lnΓ_k(mix) - lnΓ_k(pure i))

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        mole_fractions : list of float
            Mole fractions of components.
        composition : list of dict
            List of components, each represented as a dict of subgroup counts.

        Returns
        -------
        ln_gamma_r : list of float
            Residual part of the natural logarithm of activity coefficients.
        """
        try:
            if T <= 0:
                raise ValueError("T must be > 0 K")

            x = mole_fractions
            ncomp = len(x)
            if len(composition) != ncomp:
                raise ValueError(
                    "composition length must match mole_fractions length")

            sx = sum(x)
            if abs(sx - 1.0) > 1e-8:
                raise ValueError(
                    f"x must be normalized (sum=1). Got sum(x)={sx}")

            # --- all subgroups present in the problem (keep deterministic order)
            subgroup_ids = sorted({str(g) for comp in composition for g in comp.keys(
            )}, key=lambda s: int(s) if s.isdigit() else s)
            ng = len(subgroup_ids)

            # arrays for Q_k, main-group ids per subgroup
            Qk = [self.sub_Q[g] for g in subgroup_ids]
            main_of = [self.sub_to_main[g] for g in subgroup_ids]

            # build ν_ik matrix aligned to subgroup_ids
            nu = [[0.0] * ng for _ in range(ncomp)]
            for i, comp in enumerate(composition):
                for g, cnt in comp.items():
                    g = str(g)
                    # ng is small typically; if you want, optimize with a dict map
                    j = subgroup_ids.index(g)
                    nu[i][j] = float(cnt)

            # helper: interaction lookup (main-group based)
            def a_main(m_main: str, n_main: str) -> float:
                try:
                    return float(self.interaction_data[m_main][n_main])
                except KeyError:
                    # missing => assume 0 (common in incomplete tables)
                    return 0.0

            # Psi between subgroups uses their main-groups
            Psi = [[math.exp(-a_main(main_of[m], main_of[n]) / T)
                    for n in range(ng)] for m in range(ng)]

            def lnGamma_from_group_counts(nu_counts: List[float]) -> List[float]:
                total_groups = sum(nu_counts)
                if total_groups <= 0:
                    return [0.0] * ng

                X = [
                    nu_counts[k] / (total_groups + self.eps) for k in range(ng)
                ]
                denom = sum(Qk[k] * X[k] for k in range(ng))
                Theta = [
                    (Qk[k] * X[k]) / (denom + self.eps) for k in range(ng)
                ]

                # D[n] = Σ_p Theta_p Psi[p][n]
                D = [
                    sum(Theta[p] * Psi[p][n] for p in range(ng)) + self.eps for n in range(ng)
                ]
                # S[m] = Σ_n Theta_n Psi[n][m]
                S = [
                    sum(Theta[n] * Psi[n][m] for n in range(ng)) + self.eps for m in range(ng)
                ]

                lnG = [0.0] * ng
                for m in range(ng):
                    second_sum = sum(
                        Theta[n] * Psi[m][n] / D[n] for n in range(ng)
                    )
                    lnG[m] = Qk[m] * (1.0 - math.log(S[m]) - second_sum)
                return lnG

            # mixture group counts: ν_mix[k] = Σ_i x_i ν_ik
            nu_mix = [0.0] * ng
            for i in range(ncomp):
                if x[i] <= 0:
                    continue
                for k in range(ng):
                    nu_mix[k] += x[i] * nu[i][k]

            lnG_mix = lnGamma_from_group_counts(nu_mix)

            # pure-i lnGamma
            lnG_pure = [lnGamma_from_group_counts(nu[i]) for i in range(ncomp)]

            # lnγ_i^R
            ln_gamma_r = [0.0] * ncomp
            for i in range(ncomp):
                s = 0.0
                for k in range(ng):
                    vik = nu[i][k]
                    if vik:
                        s += vik * (lnG_mix[k] - lnG_pure[i][k])
                ln_gamma_r[i] = s

            return ln_gamma_r
        except Exception as e:
            logger.error(f"Error in residual_term: {e}")
            raise

    def activity_coefficients(
            self,
            T: float,
            mole_fractions: List[float],
            composition: List[Dict[str, float]]
    ) -> List[float]:
        """
        Calculate full activity coefficients using UNIFAC model.

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        mole_fractions : list of float
            Mole fractions of components.
        composition : list of dict
            List of components, each represented as a dict of subgroup counts.

        Returns
        -------
        gamma : list of float
            Activity coefficients of components.
        """
        try:
            # component r,q
            r = []
            q = []

            # NOTE: calculate r and q for each component
            for comp_groups in composition:
                Ri, Qi = self.calculate_component_RQ(comp_groups)
                r.append(Ri)
                q.append(Qi)

            # SECTION: lnγ = lnγ^C + lnγ^R
            # NOTE: combinatorial part
            lnC = self.combinatorial_term(
                mole_fractions,
                r,
                q,
                z=10,
                modified=self.modified,
                eps=self.eps
            )

            # NOTE: residual part
            lnR = self.residual_term(T, mole_fractions, composition)

            # NOTE: total lnγ
            lnG = [lnC[i] + lnR[i] for i in range(len(mole_fractions))]

            # res
            return [math.exp(v) for v in lnG]
        except Exception as e:
            logger.error(f"Error in activity_coefficients: {e}")
            raise
