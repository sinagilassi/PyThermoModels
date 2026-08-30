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
        if "ENRTL" in self.datasource:
            datasource = self.datasource["ENRTL"]
        elif "enrtl" in self.datasource:
            datasource = self.datasource["enrtl"]
        else:
            datasource = {}

        datasource = {} if datasource is None else dict(datasource)
        if kwargs.get("model_input") is not None:
            datasource.update(kwargs["model_input"])

        if not isinstance(datasource, dict):
            raise ValueError("ENRTL datasource must be a dictionary")

        if len(datasource) == 0:
            raise ValueError("ENRTL datasource cannot be empty")

        return datasource

