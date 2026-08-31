# import libs
from typing import Any, Dict, List, Optional

from ..nrtl.parameter_core import NRTLParameterCore


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
