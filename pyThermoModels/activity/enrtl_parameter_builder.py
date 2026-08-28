from typing import Any, Dict, List, Literal, Optional

from ..utils.utility import TauCorrelation
from .enrtl_parameter_core import ENRTLParameterCore
from .nrtl_parameter_builder import NRTLParameterBuilder


class ENRTLParameterBuilder(ENRTLParameterCore):
    """Resolve ENRTL parameters without treating charge as fitted data."""

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
        self._nrtl_parameter_builder = NRTLParameterBuilder(
            components=components,
            comp_idx=comp_idx,
            datasource=datasource,
            equationsource=equationsource,
            **kwargs,
        )

    def inputs_generator(
        self,
        temperature: Optional[List[float | str]] = None,
        tau_correlation: TauCorrelation = "gibbs_energy",
        symbol_delimiter: Literal["|", "_"] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        include_alpha: bool = True,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        return self._nrtl_parameter_builder.inputs_generator(
            temperature=temperature,
            tau_correlation=tau_correlation,
            symbol_delimiter=symbol_delimiter,
            mixture_ids=mixture_ids,
            include_alpha=include_alpha,
            **kwargs,
        )

    def resolve_long_range(self, model_input: Dict[str, Any]) -> Dict[str, Any]:
        long_range = model_input.get("long_range", {})
        if not isinstance(long_range, dict):
            raise TypeError("model_input['long_range'] must be a dictionary")

        return {
            "model": long_range.get("model", "pitzer_debye_huckel"),
            "basis": long_range.get("basis", model_input.get("activity_basis", "molality")),
            "A": long_range.get("A"),
            "B": long_range.get("B"),
            "A_phi": long_range.get("A_phi"),
            "b": long_range.get("b"),
            "ion_size": long_range.get("ion_size"),
        }
