# import libs
import logging
from typing import Dict, List, Literal, Optional

# locals
from ..utils.utility import (
    TauCorrelation,
    map_uniquac_tau_correlation_to_method,
)
from .uniquac_local_composition import UNIQUACLocalComposition
from .uniquac_parameter_core import UNIQUACParameterCore

# NOTE: setup logger
logger = logging.getLogger(__name__)


class UNIQUACParameterBuilder(UNIQUACParameterCore):
    """
    Build UNIQUAC model inputs from datasource and runtime model input.

    The builder owns the mapping from public tau-correlation names to
    UNIQUAC-specific local-composition methods. It keeps `UNIQUAC` focused on
    calculation orchestration and result formatting.
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
        Initialize the UNIQUAC parameter builder.

        Parameters
        ----------
        components : List[str]
            Component identifiers in model calculation order.
        comp_idx : Dict[str, int]
            Mapping from component identifier to matrix/vector index.
        datasource : Dict
            Data source for UNIQUAC model parameters.
        equationsource : Dict
            Equation source reserved for equation-driven parameter extraction.
        **kwargs : dict
            Additional keyword arguments passed through from `UNIQUAC`.
        """
        # LINK: init UNIQUACParameterCore
        UNIQUACParameterCore.__init__(
            self,
            components=components,
            comp_idx=comp_idx,
            datasource=datasource,
            equationsource=equationsource,
            **kwargs
        )

        # SECTION: local composition model
        self.local_composition_model = UNIQUACLocalComposition(
            components=self.components,
            component_idx=self.comp_idx
        )
        # ! set local composition access methods
        # ? calculate dU_ij
        self.cal_dU_ij_M1 = self.local_composition_model.cal_dU_ij_M1
        # ? calculate tau_ij
        self.cal_tau_ij_M1 = self.local_composition_model.cal_tau_ij_M1
        self.cal_tau_ij_M2 = self.local_composition_model.cal_tau_ij_M2
        self.cal_tau_ij_M3 = self.local_composition_model.cal_tau_ij_M3
        self.cal_tau_ij_M4 = self.local_composition_model.cal_tau_ij_M4
        self.cal_tau_ij_M5 = self.local_composition_model.cal_tau_ij_M5
        self.cal_tau_ij_M6 = self.local_composition_model.cal_tau_ij_M6

    def __repr__(self) -> str:
        return (
            f"UNIQUACParameterBuilder("
            f"components={self.components}, comp_idx={self.comp_idx})"
        )

    def inputs_generator(
        self,
        temperature: Optional[List[float | str]] = None,
        tau_correlation: TauCorrelation = "extended_temperature",
        symbol_delimiter: Literal["|", "_"] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        components_ids: Optional[Dict[str, List[str]]] = None,
        include_pure_component_parameters: bool = True,
        **kwargs
    ):
        """
        Prepare inputs for the UNIQUAC activity model.

        Parameters
        ----------
        temperature : Optional[List[float | str]]
            Temperature as `[value, unit]`; required when tau must be
            calculated from dU or coefficient matrices.
        tau_correlation : TauCorrelation
            Descriptive tau-correlation name. The UNIQUAC default is
            `extended_temperature`, which maps to method M2.
        symbol_delimiter : Literal["|", "_"]
            Delimiter used for component-pair dictionary keys.
        mixture_ids : Optional[Dict[str, str]]
            Mixture identifiers used for datasource lookup.
        components_ids : Optional[Dict[str, List[str]]]
            Component identifiers used for r/q component-record fallback.
        include_pure_component_parameters : bool
            If True, extract `r_i` and `q_i`. Tau-only calculators set this to
            False because pure-component parameters are not needed for tau_ij.
        **kwargs : dict
            Additional runtime values, especially `model_input`.

        Returns
        -------
        Dict[str, Any]
            UNIQUAC input dictionary containing `tau_ij`, `r_i`, `q_i`, and
            any intermediate coefficient matrices available from the source.
        """
        try:
            # NOTE: convert a required ij matrix and attach a parameter-specific error
            def _require_matrix(matrix, name):
                if matrix is None:
                    raise ValueError(
                        f"{name} is required for tau_correlation "
                        f"{tau_correlation!r} but is None."
                    )
                return self.to_matrix_ij_or(
                    data=matrix,
                    property_name=name,
                    symbol_delimiter=symbol_delimiter
                )

            # SECTION: extract source parameter values
            parameters = self.extract_parameter_values(
                mixture_ids=mixture_ids,
                components_ids=components_ids,
                symbol_delimiter=symbol_delimiter,
                include_pure_component_parameters=include_pure_component_parameters,
                **kwargs
            )

            tau_ij = parameters.get('tau_ij')
            dU_ij = parameters.get('dU_ij')
            a_ij = parameters.get('a_ij')
            b_ij = parameters.get('b_ij')
            c_ij = parameters.get('c_ij')
            d_ij = parameters.get('d_ij')
            r_i = parameters.get('r_i')
            q_i = parameters.get('q_i')
            tau_ij_cal_method = parameters.get('tau_ij_cal_method', 0)

            # SECTION: build tau_ij when it is not supplied directly
            if tau_ij is None:
                T_K = self.validate_temperature(temperature=temperature, unit='K')
                tau_method = map_uniquac_tau_correlation_to_method(
                    tau_correlation
                )

                if tau_ij_cal_method == 1 and tau_method == 'M1':
                    # NOTE: gibbs_energy maps to dU-based UNIQUAC M1
                    dU_ij = _require_matrix(dU_ij, "dU_ij")
                    tau_ij, _ = self.cal_tau_ij_M1(
                        temperature=T_K,
                        dU_ij=dU_ij,
                        symbol_delimiter=symbol_delimiter
                    )
                elif tau_ij_cal_method == 2:
                    # NOTE: coefficient-based UNIQUAC correlations use ln(tau)
                    if tau_method == 'M2':
                        tau_ij, _ = self.cal_tau_ij_M2(
                            temperature=T_K,
                            a_ij=_require_matrix(a_ij, "a_ij"),
                            b_ij=_require_matrix(b_ij, "b_ij"),
                            c_ij=_require_matrix(c_ij, "c_ij"),
                            d_ij=_require_matrix(d_ij, "d_ij"),
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_method == 'M4':
                        tau_ij, _ = self.cal_tau_ij_M4(
                            temperature=T_K,
                            a_ij=_require_matrix(a_ij, "a_ij"),
                            b_ij=_require_matrix(b_ij, "b_ij"),
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_method == 'M5':
                        tau_ij, _ = self.cal_tau_ij_M5(
                            temperature=T_K,
                            a_ij=_require_matrix(a_ij, "a_ij"),
                            b_ij=_require_matrix(b_ij, "b_ij"),
                            c_ij=_require_matrix(c_ij, "c_ij"),
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_method == 'M6':
                        tau_ij, _ = self.cal_tau_ij_M6(
                            temperature=T_K,
                            a_ij=_require_matrix(a_ij, "a_ij"),
                            b_ij=_require_matrix(b_ij, "b_ij"),
                            c_ij=_require_matrix(c_ij, "c_ij"),
                            symbol_delimiter=symbol_delimiter
                        )
                    else:
                        raise ValueError(
                            f"tau_correlation {tau_correlation!r} requires "
                            "dU_ij for UNIQUAC."
                        )
                else:
                    raise ValueError(
                        "tau_ij_cal_method not supported for selected "
                        "UNIQUAC tau_correlation."
                    )
            else:
                # NOTE: normalize directly supplied tau to an ordered matrix
                tau_ij = _require_matrix(tau_ij, "tau_ij")

            # SECTION: package model inputs
            inputs = {
                "dU_ij": dU_ij,
                "tau_ij": tau_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
            }

            if include_pure_component_parameters:
                inputs["r_i"] = r_i
                inputs["q_i"] = q_i

            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate UNIQUAC activity inputs: {e}"
            ) from e
