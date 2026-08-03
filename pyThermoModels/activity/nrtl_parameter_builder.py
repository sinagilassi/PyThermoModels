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
from ..utils.utility import TauCorrelation, map_tau_correlation_to_method
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

    def __repr__(self) -> str:
        return f"NRTLParameterBuilder(components={self.components}, comp_idx={self.comp_idx})"

    # SECTION: Manage parameters

    def manage_parameters(
            self,
    ):
        pass

    # SECTION: inputs generator

    def inputs_generator(
        self,
        temperature: Optional[
            List[float | str]
        ] = None,
        tau_correlation: TauCorrelation = "gibbs_energy",
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        include_alpha: bool = True,
        **kwargs
    ):
        '''
        Prepares inputs for the NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        tau_correlation : TauCorrelation, optional
            Descriptive tau-correlation name. Default is `gibbs_energy`,
            which maps to NRTL method M1.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for component-pair dictionary keys. Default is "|".
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            tau_method = map_tau_correlation_to_method(tau_correlation)

            # NOTE: utils
            def _require_matrix(matrix, name):
                if matrix is None:
                    raise ValueError(
                        f"{name} is required for tau_correlation '{tau_correlation}' but is None."
                    )
                return self.to_matrix_ij_or(
                    data=matrix,
                    property_name=name,
                    symbol_delimiter=symbol_delimiter
                )

            # SECTION: temperature validation
            T_K = self.validate_temperature(
                temperature=temperature,
                unit='K',
            )

            # SECTION: extract parameters
            parameters = self.extract_parameter_values(
                mixture_ids=mixture_ids,
                symbol_delimiter=symbol_delimiter,
                include_alpha=include_alpha,
                **kwargs
            )
            # >>> unpack parameters
            tau_ij = parameters.get('tau_ij')
            dg_ij = parameters.get('dg_ij')
            a_ij = parameters.get('a_ij')
            b_ij = parameters.get('b_ij')
            c_ij = parameters.get('c_ij')
            d_ij = parameters.get('d_ij')
            alpha_ij = parameters.get('alpha_ij')
            # calculation method for tau_ij
            tau_ij_cal_method: int = parameters.get('tau_ij_cal_method', 0)

            # SECTION: calculate the binary interaction parameter matrix (tau_ij)
            # check
            if tau_ij is None:
                # NOTE: check method
                if (
                    tau_ij_cal_method == 1 and
                    tau_method == 'M1'
                ):
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
                    }[tau_method]

                    missing_arrays = [
                        name for name, value in required_arrays.items()
                        if value is None
                    ]

                    if missing_arrays:
                        raise ValueError(
                            f"Missing required coefficient matrix/matrices for tau_correlation '{tau_correlation}': {', '.join(missing_arrays)}."
                        )

                    # NOTE: calculate tau_ij based on the selected correlation
                    if tau_method == 'M3':
                        # ! M3
                        a_ij_m3 = _require_matrix(a_ij, "a_ij")
                        b_ij_m3 = _require_matrix(b_ij, "b_ij")

                        # >> calculate tau_ij
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M3(
                            temperature=T_K,
                            a_ij=a_ij_m3,
                            b_ij=b_ij_m3,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_method == 'M4':
                        # ! M4
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
                    elif tau_method == 'M5':
                        # ! M5
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
                        # ! default case: M2
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
                else:
                    raise ValueError(
                        "tau_ij_cal_method not supported!"
                    )
            else:
                # NOTE: tau_ij is provided, validate it
                tau_ij = _require_matrix(tau_ij, "tau_ij")

            # NOTE: nrtl inputs
            inputs = {
                "tau_ij": tau_ij,
            }

            if include_alpha:
                # NOTE: validate alpha_ij only for full NRTL model inputs
                alpha_ij = _require_matrix(alpha_ij, "alpha_ij")
                inputs["alpha_ij"] = alpha_ij

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate NRTL activity inputs: {e}") from e
