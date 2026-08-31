# import libs
from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional

import numpy as np

from ..nrtl.parameter_core import NRTLParameterCore

ENRTLRole = Literal["solvent", "neutral_solute", "cation", "anion"]
ENRTLInteractionType = Literal[
    "self",
    "molecule_molecule",
    "molecule_cation",
    "molecule_anion",
    "cation_molecule",
    "anion_molecule",
    "cation_anion",
    "anion_cation",
    "like_cation",
    "like_anion",
]


@dataclass(frozen=True)
class ENRTLComponentInfo:
    """Thermodynamic role assigned to a true-species component."""

    key: str
    charge: int
    role: ENRTLRole


@dataclass(frozen=True)
class ENRTLInteractionParameters:
    """Chen-Evans local-composition parameter structure.

    `tau[i, j]` and `alpha[i, j]` use the same index convention as the
    existing NRTL implementation: row `i` is the neighboring/local species and
    column `j` is the central species for the NRTL-like weighting
    `G[i, j] = exp(-alpha[i, j] * tau[i, j])`.

    `interaction_mask[i, j]` marks interaction classes allowed by the
    Chen-Evans like-ion repulsion postulate. Like-ion and self terms are not
    used by the ionic local-composition kernel.
    """

    components: List[ENRTLComponentInfo]
    tau: np.ndarray
    alpha: np.ndarray
    interaction_mask: np.ndarray
    interaction_type: np.ndarray

    @classmethod
    def from_tau_alpha(
        cls,
        component_keys: List[str],
        charges: np.ndarray,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
    ) -> "ENRTLInteractionParameters":
        """Build and validate Chen-Evans interaction metadata."""
        comp_num = len(component_keys)
        if charges.shape != (comp_num,):
            raise ValueError(f"charges must have shape ({comp_num},)")
        if tau_ij.shape != (comp_num, comp_num):
            raise ValueError(f"tau_ij must have shape ({comp_num}, {comp_num})")
        if alpha_ij.shape != (comp_num, comp_num):
            raise ValueError(f"alpha_ij must have shape ({comp_num}, {comp_num})")
        if not np.all(np.isfinite(tau_ij)):
            raise ValueError("tau_ij contains non-finite values")
        if not np.all(np.isfinite(alpha_ij)):
            raise ValueError("alpha_ij contains non-finite values")

        components = [
            ENRTLComponentInfo(
                key=component_keys[i],
                charge=int(charges[i]),
                role=cls._role_from_charge(int(charges[i])),
            )
            for i in range(comp_num)
        ]
        interaction_type = np.empty((comp_num, comp_num), dtype=object)
        interaction_mask = np.ones((comp_num, comp_num), dtype=bool)

        for i in range(comp_num):
            for j in range(comp_num):
                kind = cls._interaction_type(components[i].role, components[j].role, i == j)
                interaction_type[i, j] = kind
                if kind in ("self", "like_cation", "like_anion"):
                    interaction_mask[i, j] = False

        prohibited = ~interaction_mask
        np.fill_diagonal(prohibited, False)
        if np.any(np.abs(tau_ij[prohibited]) > 0.0):
            bad_i, bad_j = np.argwhere(np.abs(tau_ij * prohibited) > 0.0)[0]
            raise ValueError(
                "Chen-Evans 1986 prohibits like-ion local-composition "
                "interaction parameters; nonzero tau_ij found for "
                f"{component_keys[bad_i]} -> {component_keys[bad_j]}."
            )

        return cls(
            components=components,
            tau=np.asarray(tau_ij, dtype=float),
            alpha=np.asarray(alpha_ij, dtype=float),
            interaction_mask=interaction_mask,
            interaction_type=interaction_type,
        )

    @staticmethod
    def _role_from_charge(charge: int) -> ENRTLRole:
        if charge > 0:
            return "cation"
        if charge < 0:
            return "anion"
        return "neutral_solute"

    @staticmethod
    def _interaction_type(
        neighbor_role: ENRTLRole,
        central_role: ENRTLRole,
        is_self: bool,
    ) -> ENRTLInteractionType:
        if is_self:
            return "self"
        if neighbor_role in ("solvent", "neutral_solute"):
            if central_role in ("solvent", "neutral_solute"):
                return "molecule_molecule"
            if central_role == "cation":
                return "molecule_cation"
            return "molecule_anion"
        if neighbor_role == "cation":
            if central_role == "cation":
                return "like_cation"
            if central_role == "anion":
                return "cation_anion"
            return "cation_molecule"
        if central_role == "anion":
            return "like_anion"
        if central_role == "cation":
            return "anion_cation"
        return "anion_molecule"


class ENRTLParameterCore(NRTLParameterCore):
    """Parameter-source utilities for ENRTL.

    ENRTL reuses the NRTL temperature-correlation and matrix conversion
    infrastructure, while electrolyte-only parameters stay separate from
    ordinary binary interaction matrices.
    """

    def __init__(
        self,
        components: List[str],
        comp_idx: Dict[str, int],
        datasource: Dict,
        equationsource: Dict,
        **kwargs: Any,
    ) -> None:
        """
        Initialize the ENRTL parameter core for a fixed component set.

        Parameters
        ----------
        components : List[str]
            Ordered list of true-species component keys.
        comp_idx : Dict[str, int]
            Mapping of component key to its index/position.
        datasource : Dict
            Data source used to resolve model parameters.
        equationsource : Dict
            Equation source used to resolve temperature-dependent correlations.
        **kwargs : Any
            Additional keyword arguments forwarded to `NRTLParameterCore`.
        """
        super().__init__(
            components=components,
            comp_idx=comp_idx,
            datasource=datasource,
            equationsource=equationsource,
            **kwargs,
        )

    def data_source_generator(
        self,
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Resolve the ENRTL-specific datasource block, merging in any override
        supplied through `model_input`.

        Parameters
        ----------
        mixture_ids : Optional[Dict[str, str]]
            Mapping of component key to mixture identifier, if applicable.
        **kwargs : Any
            May include `model_input`, a dict merged on top of the resolved
            datasource (used for run-time parameter overrides).

        Returns
        -------
        Dict[str, Any]
            Resolved, non-empty ENRTL datasource dictionary.

        Raises
        ------
        ValueError
            If the resolved datasource is not a dictionary or is empty.
        """
        # SECTION: resolve the ENRTL datasource block (case-insensitive key lookup)
        if "ENRTL" in self.datasource:
            datasource = self.datasource["ENRTL"]
        elif "enrtl" in self.datasource:
            datasource = self.datasource["enrtl"]
        else:
            datasource = {}

        # NOTE: model_input overrides take precedence over the static datasource
        datasource = {} if datasource is None else dict(datasource)
        if kwargs.get("model_input") is not None:
            datasource.update(kwargs["model_input"])

        if not isinstance(datasource, dict):
            raise ValueError("ENRTL datasource must be a dictionary")

        if len(datasource) == 0:
            raise ValueError("ENRTL datasource cannot be empty")

        return datasource
