# import libs
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union, Literal, TypeAlias, cast
import pycuc
from pyThermoDB import (
    TableMatrixData,
    TableData,
    TableEquation,
    TableMatrixEquation
)
# locals
from .nrtl_parameter_core import NRTLParameterCore
from .nrtl_local_composition import NRTLLocalComposition
from .component_parameter_mixin import ComponentParameterMixin


# NOTE: set up logger for this module
logger = logging.getLogger(__name__)


class NRTLParameterBuilder(NRTLParameterCore):

    def __init__(
        self,
        components: List[str],
        comp_idx: Dict[str, int],
        datasource: Dict,
        equationsource: Dict,
        **kwargs
    ):
        # LINK: init NRTLParameterCore
        NRTLParameterCore.__init__(
            self,
            components=components,
            comp_idx=comp_idx,
            datasource=datasource,
            equationsource=equationsource,
            **kwargs
        )

        # SECTION: local composition model
        self.local_composition_model = NRTLLocalComposition(
            components=self.components,
            component_idx=self.comp_idx
        )
        # ! set local composition access methods
        # ? calculate dg_ij
        self.cal_dg_ij_M1 = self.local_composition_model.cal_dg_ij_M1
        # ? calculate tau_ij
        self.cal_tau_ij_M1 = self.local_composition_model.cal_tau_ij_M1
        self.cal_tau_ij_M2 = self.local_composition_model.cal_tau_ij_M2
        self.cal_tau_ij_M3 = self.local_composition_model.cal_tau_ij_M3
        self.cal_tau_ij_M4 = self.local_composition_model.cal_tau_ij_M4
        self.cal_tau_ij_M5 = self.local_composition_model.cal_tau_ij_M5

    # SECTION: inputs generator
    def inputs_generator(
        self,
        temperature: Optional[
            List[float | str]
        ] = None,
        tau_correlation: Literal[
            'M1', 'M2', 'M3', 'M4', 'M5'
        ] = 'M1',
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs
    ):
        '''
        Prepares inputs for the NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for component-pair dictionary keys. Default is "|".
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: data source
            # generate datasource
            datasource: Dict[str, Any] = self.data_source_generator(
                mixture_ids=mixture_ids,
                **kwargs
            )

            # ! set initial values
            a_ij = None
            b_ij = None
            c_ij = None
            d_ij = None
            dg_ij = None
            alpha_ij = None
            tau_ij = None

            def _has_source(src) -> bool:
                return src is not None and src != 'None'

            def _require_matrix(
                value: Optional[np.ndarray],
                parameter_name: str
            ) -> np.ndarray:
                if value is None:
                    raise ValueError(f"{parameter_name} cannot be None")
                return value

            # NOTE: check temperature
            # init temperature [K]
            T_K = self._validate_temperature(temperature, unit='K')

            # SECTION: extract parameter source
            parameter_sources = self.extract_parameter_sources(datasource)
            # >> unpack
            dg_ij_src = parameter_sources['dg_ij_src']
            a_ij_src = parameter_sources['a_ij_src']
            b_ij_src = parameter_sources['b_ij_src']
            c_ij_src = parameter_sources['c_ij_src']
            d_ij_src = parameter_sources['d_ij_src']
            alpha_ij_src = parameter_sources['alpha_ij_src']
            tau_ij_src = parameter_sources['tau_ij_src']

            # SECTION: check correlation methods

            allowed_tau_correlations = ('M1', 'M2', 'M3', 'M4', 'M5')
            allowed_dg_correlations = ('M1',)

            if tau_correlation not in allowed_tau_correlations:
                raise ValueError(
                    f"tau_correlation must be one of {allowed_tau_correlations}."
                )

            # SECTION: extract data

            # NOTE: check if dg_ij is provided
            if (
                tau_ij_src is not None
            ):
                # ! use tau_ij
                # set method
                tau_ij_cal_method = 0

            elif (
                tau_correlation == 'M1' and
                _has_source(dg_ij_src)
            ):  # >>> check if dg_ij is provided
                # ! use dg_ij
                # set method
                tau_ij_cal_method = 1

                # >> check
                if dg_ij_src is not None:
                    dg_ij = self.to_matrix_ij(
                        dg_ij_src,
                        property_name='dg',
                        symbol_delimiter=symbol_delimiter
                    )
                else:
                    raise ValueError(
                        "dg_ij must be provided for tau_correlation 'M1'."
                    )

            else:
                # ! use constants based on the selected correlation
                tau_ij_cal_method = 2

                tau_required_sources = {
                    'M2': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src, 'd_ij': d_ij_src},
                    'M3': {'a_ij': a_ij_src, 'b_ij': b_ij_src},
                    'M4': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src},
                    'M5': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src},
                }

                required_sources = tau_required_sources[tau_correlation]
                missing_sources = [
                    name for name, src in required_sources.items()
                    if not _has_source(src)
                ]

                if missing_sources:
                    raise ValueError(
                        f"Missing required source(s) for tau_correlation '{tau_correlation}': {', '.join(missing_sources)}."
                    )

                # ! extract a_ij values
                if _has_source(a_ij_src):
                    a_ij = self.to_matrix_ij_or(
                        data=a_ij_src,
                        property_name='a',
                        symbol_delimiter=symbol_delimiter
                    )

                # ! extract b_ij values
                if _has_source(b_ij_src):
                    b_ij = self.to_matrix_ij_or(
                        data=b_ij_src,
                        property_name='b',
                        symbol_delimiter=symbol_delimiter
                    )

                # ! extract c_ij values
                if _has_source(c_ij_src):
                    c_ij = self.to_matrix_ij_or(
                        data=c_ij_src,
                        property_name='c',
                        symbol_delimiter=symbol_delimiter
                    )

                # ! extract d_ij values
                if _has_source(d_ij_src):
                    d_ij = self.to_matrix_ij_or(
                        data=d_ij_src,
                        property_name='d',
                        symbol_delimiter=symbol_delimiter
                    )

            # SECTION: extract data
            # NOTE: α_ij, non-randomness parameter
            # check
            if alpha_ij_src is not None:
                # ! extract α_ij values
                alpha_ij = self.to_matrix_ij_or(
                    data=alpha_ij_src,
                    property_name='alpha',
                    symbol_delimiter=symbol_delimiter
                )
            else:
                # set default value
                alpha_ij = None

            # SECTION: calculate the binary interaction parameter matrix (tau_ij)
            # check
            if tau_ij_src is None:

                # NOTE: check method
                if tau_ij_cal_method == 1:
                    # ! >>> method 1: calculate tau_ij from dg_ij
                    # Check if dg_ij is None and convert values to float if needed
                    if dg_ij is None:
                        raise ValueError(
                            "dg_ij cannot be None for calculating tau_ij"
                        )

                    # >> calculate
                    tau_ij, tau_ij_comp = self.cal_tau_ij_M1(
                        temperature=T_K,
                        dg_ij=dg_ij,
                        symbol_delimiter=symbol_delimiter
                    )
                elif tau_ij_cal_method == 2:
                    # ! >>> method 2: calculate tau_ij from constants a, b, c, d based on the selected correlation
                    required_arrays = {
                        'M2': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij, 'd_ij': d_ij},
                        'M3': {'a_ij': a_ij, 'b_ij': b_ij},
                        'M4': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij},
                        'M5': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij},
                    }[tau_correlation]

                    missing_arrays = [
                        name for name, value in required_arrays.items()
                        if value is None
                    ]

                    if missing_arrays:
                        raise ValueError(
                            f"Missing required coefficient matrix/matrices for tau_correlation '{tau_correlation}': {', '.join(missing_arrays)}."
                        )

                    # NOTE: calculate tau_ij based on the selected correlation
                    if tau_correlation == 'M2':
                        a_ij_m2 = _require_matrix(a_ij, "a_ij")
                        b_ij_m2 = _require_matrix(b_ij, "b_ij")
                        c_ij_m2 = _require_matrix(c_ij, "c_ij")
                        d_ij_m2 = _require_matrix(d_ij, "d_ij")

                        # >> calculate tau_ij
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M2(
                            temperature=T_K,
                            a_ij=a_ij_m2,
                            b_ij=b_ij_m2,
                            c_ij=c_ij_m2,
                            d_ij=d_ij_m2,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M3':
                        a_ij_m3 = _require_matrix(a_ij, "a_ij")
                        b_ij_m3 = _require_matrix(b_ij, "b_ij")

                        # >> calculate tau_ij
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M3(
                            temperature=T_K,
                            a_ij=a_ij_m3,
                            b_ij=b_ij_m3,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M4':
                        a_ij_m4 = _require_matrix(a_ij, "a_ij")
                        b_ij_m4 = _require_matrix(b_ij, "b_ij")
                        c_ij_m4 = _require_matrix(c_ij, "c_ij")

                        # >> calculate tau_ij
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M4(
                            temperature=T_K,
                            a_ij=a_ij_m4,
                            b_ij=b_ij_m4,
                            c_ij=c_ij_m4,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M5':
                        a_ij_m5 = _require_matrix(a_ij, "a_ij")
                        b_ij_m5 = _require_matrix(b_ij, "b_ij")
                        c_ij_m5 = _require_matrix(c_ij, "c_ij")

                        # >> calculate tau_ij
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M5(
                            temperature=T_K,
                            a_ij=a_ij_m5,
                            b_ij=b_ij_m5,
                            c_ij=c_ij_m5,
                            symbol_delimiter=symbol_delimiter
                        )
                    else:
                        raise ValueError(
                            f"tau_correlation '{tau_correlation}' not supported!"
                        )
                else:
                    raise ValueError(
                        "tau_ij_cal_method not supported!"
                    )
            else:
                # ! check if tau_ij is provided
                # >>> method 0
                # check types
                tau_ij = self.to_matrix_ij_or(
                    data=tau_ij_src,
                    property_name='tau',
                    symbol_delimiter=symbol_delimiter
                )

            # NOTE: nrtl inputs
            inputs = {
                "alpha_ij": alpha_ij,
                "tau_ij": tau_ij,
                "dg_ij": dg_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
            }

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate NRTL activity inputs: {e}") from e
