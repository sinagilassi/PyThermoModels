import numpy as np
import pytest

from pyThermoModels.activity.enrtl.local_composition import ENRTLLocalComposition


def test_neutral_nrtl_limit_returns_expected_binary_ln_gamma():
    model = ENRTLLocalComposition(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
    )
    x = np.asarray([0.4, 0.6])
    tau = np.asarray([[0.0, 0.2], [0.1, 0.0]])
    alpha = np.asarray([[0.0, 0.3], [0.3, 0.0]])

    result = model.cal_ln_gamma_lc(
        mole_fraction=x,
        charges=np.asarray([0, 0]),
        tau_ij=tau,
        alpha_ij=alpha,
        mode="chen_evans_1986",
    )

    assert result.shape == (2,)
    assert result.tolist() == pytest.approx([0.10621866, 0.04584418])


def test_ionic_chen_evans_path_returns_finite_local_contribution():
    model = ENRTLLocalComposition(
        components=["H2O-l", "Na{+}-aq", "Cl{-}-aq"],
        comp_idx={"H2O-l": 0, "Na{+}-aq": 1, "Cl{-}-aq": 2},
    )

    tau = np.asarray([
        [0.0, 0.10, 0.20],
        [0.15, 0.00, 0.30],
        [0.25, 0.40, 0.00],
    ])
    alpha = np.asarray([
        [0.0, 0.2, 0.2],
        [0.2, 0.0, 0.2],
        [0.2, 0.2, 0.0],
    ])

    result = model.cal_ln_gamma_lc(
        mole_fraction=np.asarray([0.98, 0.01, 0.01]),
        charges=np.asarray([0, 1, -1]),
        tau_ij=tau,
        alpha_ij=alpha,
    )

    assert result.shape == (3,)
    assert np.all(np.isfinite(result))
    assert "local_electroneutrality_residual" in model.last_diagnostics


def test_ionic_chen_evans_rejects_like_ion_parameters():
    model = ENRTLLocalComposition(
        components=["H2O-l", "Na{+}-aq", "K{+}-aq", "Cl{-}-aq"],
        comp_idx={"H2O-l": 0, "Na{+}-aq": 1, "K{+}-aq": 2, "Cl{-}-aq": 3},
    )
    tau = np.zeros((4, 4))
    tau[1, 2] = 0.1

    with pytest.raises(ValueError, match="like-ion"):
        model.cal_ln_gamma_lc(
            mole_fraction=np.asarray([0.96, 0.01, 0.01, 0.02]),
            charges=np.asarray([0, 1, 1, -1]),
            tau_ij=tau,
            alpha_ij=np.zeros((4, 4)),
        )
