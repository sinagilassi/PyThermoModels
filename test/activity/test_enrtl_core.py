import numpy as np
import pytest

from pyThermoModels.activity.enrtl_core import ENRTLCore


def test_enrtl_electroneutrality_accepts_balanced_true_species():
    core = ENRTLCore()

    assert core.check_electroneutrality(
        composition=np.asarray([0.98, 0.01, 0.01]),
        charges=np.asarray([0, 1, -1]),
    )


def test_enrtl_electroneutrality_rejects_unbalanced_true_species():
    core = ENRTLCore()

    with pytest.raises(ValueError, match="not electrically neutral"):
        core.check_electroneutrality(
            composition=np.asarray([0.97, 0.02, 0.01]),
            charges=np.asarray([0, 1, -1]),
        )


def test_enrtl_ionic_strength_uses_squared_charge():
    core = ENRTLCore()

    value = core.ionic_strength(
        composition=np.asarray([0.1, 0.2]),
        charges=np.asarray([2, -1]),
        basis="molality",
    )

    assert value == pytest.approx(0.3)


def test_enrtl_combine_ln_gamma_sums_contributions():
    core = ENRTLCore()

    result = core.combine_ln_gamma(
        ln_gamma_long_range=np.asarray([0.1, -0.2]),
        ln_gamma_local_composition=np.asarray([0.3, 0.4]),
    )

    assert result.tolist() == pytest.approx([0.4, 0.2])
