"""Excess-Gibbs diagnostics for implemented ENRTL contributions."""

from typing import Any, Dict, Literal

import numpy as np

from .local_composition import ENRTLLocalComposition


# SECTION: Local-composition excess-Gibbs evaluation
def local_composition_excess_gibbs_RT(
    local_composition: ENRTLLocalComposition,
    mole_fraction: np.ndarray,
    charges: np.ndarray,
    tau_ij: np.ndarray,
    alpha_ij: np.ndarray,
    mode: Literal["chen_evans_1986", "neutral_nrtl_limit"] = "chen_evans_1986",
) -> float:
    """Return the implemented local-composition contribution to ``G^E / RT``."""
    x = np.asarray(mole_fraction, dtype=float)
    z = np.asarray(charges, dtype=float)
    tau = np.asarray(tau_ij, dtype=float)
    alpha = np.asarray(alpha_ij, dtype=float)

    # NOTE: retain the local-composition class as the single equation owner.
    if np.any(z != 0):
        local_composition.cal_ln_gamma_lc(
            mole_fraction=x,
            charges=z,
            tau_ij=tau,
            alpha_ij=alpha,
            mode=mode,
        )
        gE_local_RT = local_composition.last_diagnostics.get(
            "gE_local_composition_RT"
        )
        # NOTE: diagnostics are heterogeneous, so narrow the kernel value here.
        if not isinstance(gE_local_RT, (int, float, np.floating)):
            raise TypeError(
                "Chen-Evans local-composition diagnostics must contain a numeric "
                "gE_local_composition_RT value"
            )
        return float(gE_local_RT)

    G_ij = local_composition.cal_G_ij(tau_ij=tau, alpha_ij=alpha)
    denominator = np.sum(x[:, np.newaxis] * G_ij, axis=0)
    if np.any(denominator <= 0.0) or not np.all(np.isfinite(denominator)):
        raise ValueError("NRTL excess-Gibbs denominator must be finite and positive")

    # ! This is the neutral NRTL limiting expression, not an ionic identity.
    return float(np.sum(x * np.sum(x[:, np.newaxis] * tau * G_ij, axis=0) / denominator))


# SECTION: Contribution and convention diagnostics
def build_excess_gibbs_diagnostics(
    *,
    local_composition_RT: float,
    charges: np.ndarray,
    composition_representation: str,
    local_composition_mode: str,
    long_range_model: str,
    ionic_strength_basis: str,
) -> Dict[str, Any]:
    """Describe which ENRTL excess-Gibbs contributions are currently valid."""
    ionic = bool(np.any(np.asarray(charges, dtype=float) != 0.0))

    # ? A long-range potential must be derived and validated before exposure.
    if ionic:
        total_RT = None
        total_status = "unavailable"
        total_reason = (
            "The current long-range model exposes ln(gamma) contributions but "
            "does not yet expose a validated excess-Gibbs potential."
        )
    else:
        # NOTE: all-neutral mixtures have zero ionic long-range contribution.
        total_RT = float(local_composition_RT)
        total_status = "complete"
        total_reason = None

    return {
        "local_composition": float(local_composition_RT),
        "long_range": None,
        "total": total_RT,
        "total_status": total_status,
        "total_reason": total_reason,
        "reference_state": "true_species_activity_model",
        "activity_convention": "natural_log_activity_coefficients",
        "species_basis": composition_representation,
        "local_composition_mode": local_composition_mode,
        "long_range_model": long_range_model,
        "ionic_strength_basis": ionic_strength_basis,
    }
