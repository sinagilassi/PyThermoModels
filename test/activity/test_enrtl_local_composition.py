import numpy as np
import pytest

from pyThermoModels.activity.enrtl_local_composition import ENRTLLocalComposition


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
    assert result.tolist() == pytest.approx([0.08176653, 0.03630797])


def test_ionic_chen_evans_path_fails_instead_of_using_ordinary_nrtl():
    model = ENRTLLocalComposition(
        components=["H2O-l", "Na{+}-aq", "Cl{-}-aq"],
        comp_idx={"H2O-l": 0, "Na{+}-aq": 1, "Cl{-}-aq": 2},
    )

    with pytest.raises(NotImplementedError, match="ordinary NRTL"):
        model.cal_ln_gamma_lc(
            mole_fraction=np.asarray([0.98, 0.01, 0.01]),
            charges=np.asarray([0, 1, -1]),
            tau_ij=np.zeros((3, 3)),
            alpha_ij=np.zeros((3, 3)),
        )
