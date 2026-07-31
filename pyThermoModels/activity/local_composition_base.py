# SECTION: imports
from typing import Any, Dict, List, Literal, TypeAlias, cast

import numpy as np
from pyThermoDB import TableMatrixData as _TableMatrixData

# SECTION: shared types
IJData: TypeAlias = np.ndarray | Dict[str, float] | Any
IJ_DATA_TYPES: tuple[type[Any], ...] = (
    np.ndarray,
    dict,
    cast(type[Any], _TableMatrixData),
)


class LocalCompositionBase:
    """Shared input and key handling for local-composition models."""

    # SECTION: initialization
    def __init__(
        self,
        components: List[str],
        component_idx: Dict[str, int],
        **kwargs
    ):
        # NOTE: component order and index mapping are owned by the parent model
        self.components = components
        self.comp_idx = component_idx
        self.comp_num = len(self.components)

    # SECTION: validators
    def _delimiter(self, symbol_delimiter: Literal["|", "_"]) -> str:
        if symbol_delimiter == "|":
            return " | "
        if symbol_delimiter == "_":
            return "_"
        raise ValueError("symbol_delimiter must be '|' or '_'")

    def _validate_temperature(self, temperature: float) -> None:
        # ! local-composition temperature correlations require Kelvin above zero
        if temperature <= 0:
            raise ValueError("temperature must be greater than 0 K")

    def _validate_same_source_type(self, *values: IJData) -> type:
        # ! all coefficient matrices must use the same source representation
        if not all(isinstance(value, IJ_DATA_TYPES) for value in values):
            raise TypeError(
                "all interaction parameters must be numpy array, dict or TableMatrixData"
            )

        value_type = type(values[0])
        if not all(isinstance(value, value_type) for value in values):
            raise TypeError(
                "all interaction parameters must be provided in the same format"
            )
        return value_type

    def _validate_ij_data(self, value: IJData, name: str) -> None:
        # ! single interaction matrix must be one of the supported ij formats
        if not isinstance(value, IJ_DATA_TYPES):
            raise TypeError(
                f"{name} must be numpy array, dict or TableMatrixData")

    # SECTION: pair helpers
    def _pair_key(
        self,
        i: int,
        j: int,
        symbol_delimiter: Literal["|", "_"]
    ) -> str:
        delimiter = self._delimiter(symbol_delimiter)
        return f"{self.components[i]}{delimiter}{self.components[j]}"

    def _component_position(self, i: int, j: int) -> tuple[int, int]:
        # NOTE: use parent model component ids for matrix placement
        row = self.comp_idx[self.components[i]]
        col = self.comp_idx[self.components[j]]
        return row, col

    # SECTION: data extraction
    def _table_value(
        self,
        data: Any,
        symbol: str,
        i: int,
        j: int
    ) -> float:
        component_i = self.components[i]
        component_j = self.components[j]
        # NOTE: support both underscore and dash table-matrix key conventions
        keys = (
            f"{symbol}_{component_i}_{component_j}",
            f"{symbol}_{component_i}-{component_j}",
        )

        value = None
        used_key = keys[0]
        for key in keys:
            value = data.ij(key)
            used_key = key
            if value is not None and value.get("value", None) is not None:
                return float(value["value"])

        raise ValueError(
            f"Invalid value for {symbol}: {value} for key: {used_key}"
        )

    def _ij_value(
        self,
        data: IJData,
        symbol: str,
        i: int,
        j: int,
        pair_key: str
    ) -> float:
        # NOTE: extract a numeric ij value independent of source format
        if isinstance(data, np.ndarray):
            return float(data[i, j])
        if isinstance(data, dict):
            return float(data[pair_key])
        return self._table_value(data, symbol, i, j)
