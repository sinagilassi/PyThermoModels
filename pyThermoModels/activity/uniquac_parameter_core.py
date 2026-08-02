# import libs
import logging
from typing import Any, Dict, List, Literal, Optional

import numpy as np
import pycuc

# locals
from .component_parameter_mixin import ComponentParameterMixin

# NOTE: setup logger
logger = logging.getLogger(__name__)


class UNIQUACParameterCore:
    """
    Shared parameter-source handling for UNIQUAC activity-model inputs.

    This class owns the datasource lookup, alias resolution, source extraction,
    and primitive type conversion used by `UNIQUACParameterBuilder`. The
    top-level `UNIQUAC` model should stay focused on calculation flow and
    result formatting.
    """

    def __init__(
        self,
        components: List[str],
        comp_idx: Dict[str, int],
        datasource: Dict,
        equationsource: Dict,
        **kwargs
    ):
        """
        Initialize the UNIQUAC parameter core.

        Parameters
        ----------
        components : List[str]
            Component identifiers in the model calculation order.
        comp_idx : Dict[str, int]
            Mapping from component identifier to matrix/vector position.
        datasource : Dict
            Data source containing UNIQUAC parameters or component records.
        equationsource : Dict
            Equation source reserved for equation-driven parameter extraction.
        **kwargs : dict
            Additional keyword arguments reserved for future extensions.
        """
        # SECTION: core component and source configuration
        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)
        self.datasource = datasource
        self.equationsource = equationsource

        # SECTION: component parameter converters
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx
        )
        # NOTE: expose only converters used by the parameter core
        self.to_i = self.component_parameter_mixin.to_i
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_matrix_ij_or = self.component_parameter_mixin.to_matrix_ij_or

    def data_source_generator(
        self,
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Generate a data source dictionary for UNIQUAC activity model.

        Parameters
        ----------
        mixture_ids : Optional[Dict[str, str]]
            Mixture identifiers, usually keyed by `Name` and `Formula`.
        **kwargs : dict
            Optional values, especially `model_input`, that override the
            datasource-selected values.

        Returns
        -------
        Dict[str, Any]
            Resolved datasource dictionary used for UNIQUAC input generation.
        """
        try:
            # SECTION: select datasource by model key or mixture id
            if "UNIQUAC" in self.datasource.keys():
                # NOTE: canonical uppercase model key
                datasource = self.datasource["UNIQUAC"]
            elif "uniquac" in self.datasource.keys():
                # NOTE: lowercase model key used by some sources
                datasource = self.datasource["uniquac"]
            elif mixture_ids is not None:
                datasource = {}
                # NOTE: prefer mixture Name over Formula, matching NRTL
                if 'Name' in mixture_ids.keys():
                    if mixture_ids.get('Name', None):
                        key_ = mixture_ids['Name']
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
                elif 'Formula' in mixture_ids.keys():
                    if mixture_ids.get('Formula', None):
                        key_ = mixture_ids['Formula']
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
            else:
                logger.warning(
                    "No UNIQUAC or uniquac key found in datasource, using model_input if provided."
                )
                datasource = {}

            # SECTION: validate and normalize selected datasource
            if datasource is not None and not isinstance(datasource, dict):
                raise ValueError("datasource must be a dictionary.")

            datasource = {} if datasource is None else dict(datasource)

            # NOTE: model_input values override persisted datasource values
            if kwargs.get('model_input') is not None:
                datasource.update(kwargs['model_input'])

            # SECTION: final datasource validation
            if datasource is None:
                raise ValueError("datasource cannot be None.")

            if not isinstance(datasource, dict):
                raise ValueError("datasource must be a dictionary.")

            if len(datasource) == 0:
                raise ValueError("datasource cannot be empty.")

            return datasource
        except Exception as e:
            raise Exception(
                f"Failed to generate UNIQUAC activity data source: {e}"
            ) from e

    def _extract_parameter_source(
        self,
        ids: List[str],
        datasource: Dict[str, Any]
    ) -> Optional[Any]:
        """
        Extract the first available parameter source from a list of aliases.

        Parameters
        ----------
        ids : List[str]
            Candidate parameter names ordered by precedence.
        datasource : Dict[str, Any]
            Resolved UNIQUAC datasource.

        Returns
        -------
        Optional[Any]
            First matching source value, or None when no alias exists.
        """
        try:
            # SECTION: search parameter aliases in precedence order
            src = None
            for id_ in ids:
                src = datasource.get(id_, None)
                if src is not None:
                    break
            return src
        except Exception as e:
            raise Exception(
                f"Failed to extract parameter source for UNIQUAC activity model: {e}"
            ) from e

    def _extract_component_parameter_source(
        self,
        ids: List[str],
        datasource: Dict[str, Any],
        components_ids: Optional[Dict[str, List[str]]] = None
    ) -> Optional[Any]:
        """
        Extract pure-component parameter sources with component-record fallback.

        Parameters
        ----------
        ids : List[str]
            Candidate direct parameter names, such as `r_i` then `r`.
        datasource : Dict[str, Any]
            Resolved UNIQUAC datasource.
        components_ids : Optional[Dict[str, List[str]]]
            Component identifier collections used to locate component records.

        Returns
        -------
        Optional[Any]
            Direct parameter source or a component-name dictionary.
        """
        # SECTION: direct datasource extraction
        src = self._extract_parameter_source(ids=ids, datasource=datasource)
        if src is not None:
            return src

        if components_ids is None:
            return None

        # SECTION: fallback extraction from component records
        property_name = ids[-1]
        for component_key in ("Name-State", "Formula-State"):
            if component_key not in components_ids:
                continue

            values = {}
            for component in components_ids[component_key]:
                # NOTE: strip phase/state suffix to match model component names
                component_name = component.rsplit('-', 1)[0]
                values[component_name] = self.datasource.get(
                    component, {}
                ).get(property_name, {}).get('value', None)

            if values:
                return values

        return None

    def extract_parameter_sources(
        self,
        datasource: Dict[str, Any],
        components_ids: Optional[Dict[str, List[str]]] = None,
        include_pure_component_parameters: bool = True
    ) -> Dict[str, Any]:
        """
        Extract raw UNIQUAC parameter sources from a resolved datasource.

        Parameters
        ----------
        datasource : Dict[str, Any]
            Resolved UNIQUAC datasource.
        components_ids : Optional[Dict[str, List[str]]]
            Component identifiers used for r/q component-record fallback.
        include_pure_component_parameters : bool
            If True, extract `r_i` and `q_i` sources.

        Returns
        -------
        Dict[str, Any]
            Raw source values for interaction, coefficient, and pure-component
            parameters.
        """
        # SECTION: binary interaction parameter sources
        dU_ij_src = self._extract_parameter_source(
            ids=['dU_ij', 'dU'],
            datasource=datasource
        )
        a_ij_src = self._extract_parameter_source(
            ids=['a_ij', 'a'],
            datasource=datasource
        )
        b_ij_src = self._extract_parameter_source(
            ids=['b_ij', 'b'],
            datasource=datasource
        )
        c_ij_src = self._extract_parameter_source(
            ids=['c_ij', 'c'],
            datasource=datasource
        )
        d_ij_src = self._extract_parameter_source(
            ids=['d_ij', 'd'],
            datasource=datasource
        )
        tau_ij_src = self._extract_parameter_source(
            ids=['tau_ij', 'tau'],
            datasource=datasource
        )

        # SECTION: pure-component parameter sources
        r_i_src = None
        q_i_src = None
        if include_pure_component_parameters:
            r_i_src = self._extract_component_parameter_source(
                ids=['r_i', 'r'],
                datasource=datasource,
                components_ids=components_ids
            )
            q_i_src = self._extract_component_parameter_source(
                ids=['q_i', 'q'],
                datasource=datasource,
                components_ids=components_ids
            )

        return {
            'dU_ij_src': dU_ij_src,
            'a_ij_src': a_ij_src,
            'b_ij_src': b_ij_src,
            'c_ij_src': c_ij_src,
            'd_ij_src': d_ij_src,
            'tau_ij_src': tau_ij_src,
            'r_i_src': r_i_src,
            'q_i_src': q_i_src,
        }

    def _has_source(self, src: Any) -> bool:
        """
        Check whether a source value should be treated as present.
        """
        return src is not None and src != 'None'

    def _to_i_or_none(self, data: Any, name: str) -> Optional[np.ndarray]:
        """
        Convert a pure-component parameter source to an ordered vector.

        Parameters
        ----------
        data : Any
            List, numpy array, or component-keyed dictionary source.
        name : str
            Parameter name used in validation errors.

        Returns
        -------
        Optional[np.ndarray]
            Ordered component vector or None when no source is provided.
        """
        # SECTION: optional source handling
        if not self._has_source(data):
            return None

        # SECTION: source conversion
        if isinstance(data, np.ndarray):
            return data
        if isinstance(data, list):
            return np.asarray(data, dtype=float)
        if isinstance(data, dict):
            return self.to_i(data)
        raise ValueError(f"{name} must be a list, numpy array, or dict.")

    def extract_parameter_values(
        self,
        symbol_delimiter: Literal["|", "_"] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        components_ids: Optional[Dict[str, List[str]]] = None,
        include_pure_component_parameters: bool = True,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Extract and normalize UNIQUAC parameter values.

        Parameters
        ----------
        symbol_delimiter : Literal["|", "_"]
            Delimiter used for component-pair dictionary keys.
        mixture_ids : Optional[Dict[str, str]]
            Mixture identifiers for datasource lookup.
        components_ids : Optional[Dict[str, List[str]]]
            Component identifiers for pure-component fallback extraction.
        include_pure_component_parameters : bool
            If True, extract `r_i` and `q_i`. Tau-only calculators set this to
            False because pure-component UNIQUAC parameters are not needed to
            calculate tau_ij.
        **kwargs : dict
            Optional values passed through to datasource generation.

        Returns
        -------
        Dict[str, Any]
            Normalized matrices/vectors plus `tau_ij_cal_method`.
        """
        # SECTION: resolve datasource and raw sources
        datasource = self.data_source_generator(
            mixture_ids=mixture_ids,
            **kwargs
        )

        sources = self.extract_parameter_sources(
            datasource=datasource,
            components_ids=components_ids,
            include_pure_component_parameters=include_pure_component_parameters
        )

        tau_ij = None
        dU_ij = None
        a_ij = None
        b_ij = None
        c_ij = None
        d_ij = None

        tau_ij_src = sources['tau_ij_src']
        dU_ij_src = sources['dU_ij_src']
        a_ij_src = sources['a_ij_src']
        b_ij_src = sources['b_ij_src']
        c_ij_src = sources['c_ij_src']
        d_ij_src = sources['d_ij_src']

        # SECTION: choose tau source strategy
        if self._has_source(tau_ij_src):
            # NOTE: tau is already provided, so no temperature correlation is needed
            tau_ij = self.to_matrix_ij_or(
                data=tau_ij_src,
                property_name='tau',
                symbol_delimiter=symbol_delimiter
            )
            tau_ij_cal_method = 0
        elif self._has_source(dU_ij_src):
            # NOTE: dU requires the UNIQUAC Gibbs-energy correlation
            dU_ij = self.to_matrix_ij_or(
                data=dU_ij_src,
                property_name='dU',
                symbol_delimiter=symbol_delimiter
            )
            tau_ij_cal_method = 1
        else:
            # NOTE: coefficient matrices are consumed by selected tau correlations
            if self._has_source(a_ij_src):
                a_ij = self.to_matrix_ij_or(
                    data=a_ij_src,
                    property_name='a',
                    symbol_delimiter=symbol_delimiter
                )
            if self._has_source(b_ij_src):
                b_ij = self.to_matrix_ij_or(
                    data=b_ij_src,
                    property_name='b',
                    symbol_delimiter=symbol_delimiter
                )
            if self._has_source(c_ij_src):
                c_ij = self.to_matrix_ij_or(
                    data=c_ij_src,
                    property_name='c',
                    symbol_delimiter=symbol_delimiter
                )
            if self._has_source(d_ij_src):
                d_ij = self.to_matrix_ij_or(
                    data=d_ij_src,
                    property_name='d',
                    symbol_delimiter=symbol_delimiter
                )
            tau_ij_cal_method = 2

        # SECTION: pure-component parameter vectors
        r_i = None
        q_i = None
        if include_pure_component_parameters:
            r_i = self._to_i_or_none(sources['r_i_src'], "r_i")
            q_i = self._to_i_or_none(sources['q_i_src'], "q_i")

        return {
            'tau_ij': tau_ij,
            'dU_ij': dU_ij,
            'a_ij': a_ij,
            'b_ij': b_ij,
            'c_ij': c_ij,
            'd_ij': d_ij,
            'r_i': r_i,
            'q_i': q_i,
            'tau_ij_cal_method': tau_ij_cal_method,
        }

    def validate_temperature(
        self,
        temperature: Optional[List[float | str]],
        unit: str = 'K'
    ) -> float:
        """
        Validate and convert a temperature input to the requested unit.

        Parameters
        ----------
        temperature : Optional[List[float | str]]
            Temperature as `[value, unit]`.
        unit : str
            Target temperature unit, by default Kelvin.

        Returns
        -------
        float
            Converted temperature value.
        """
        # SECTION: parse and convert temperature
        T = -1

        if temperature is not None:
            if not isinstance(temperature, list):
                raise ValueError(
                    "temperature must be a list of floats or strings."
                )
            if len(temperature) == 0:
                raise ValueError("temperature list is empty")

            T_value = float(temperature[0])
            T_unit = str(temperature[1])
            T = pycuc.convert_from_to(T_value, T_unit, unit)

        # SECTION: validate physical temperature range
        if T is None or T <= 0:
            raise ValueError(
                "temperature must be provided and greater than 0 K."
            )

        return T
