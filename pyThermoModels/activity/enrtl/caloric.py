"""Guardrail for ENRTL caloric properties pending thermodynamic validation."""

from typing import NoReturn


def excess_enthalpy_not_available() -> NoReturn:
    """Raise until ENRTL derivative identities and reference states are validated."""
    # ! Do not infer H^E from activity coefficients before the convention is fixed.
    raise NotImplementedError(
        "ENRTL excess enthalpy is not available until derivative identities and "
        "reference-state conventions are independently validated."
    )
