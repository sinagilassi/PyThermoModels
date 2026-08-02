# import libs
import logging
import numpy as np
from typing import Any, Dict, Literal, Optional, Tuple, cast, List
import pycuc
# locals
from .component_parameter_mixin import ComponentParameterMixin

# NOTE: setup logger
logger = logging.getLogger(__name__)


class NRTLParameterCore:

    def __init__(
        self,
        components: List[str],
        comp_idx: Dict[str, int],
        datasource: Dict,
        equationsource: Dict,
        **kwargs
    ):
        # set attributes
        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)

        self.datasource = datasource
        self.equationsource = equationsource

        # SECTION: component parameter mixin
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx
        )
        # ! set component access methods
        self.to_ij = self.component_parameter_mixin.to_ij
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_matrix_ij_or = self.component_parameter_mixin.to_matrix_ij_or

    # SECTION: inputs generator
    # NOTE: data source generator
    def data_source_generator(
        self,
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Generate a data source dictionary for NRTL activity model.

        Parameters
        ----------
        mixture_ids : Optional[Dict[str, str]], optional
            Dictionary containing mixture identifiers. Default is None.

        Returns
        -------
        Dict[str, Any]
            A dictionary containing the data source for the NRTL activity model.
        """
        try:
            # SECTION: check src
            # check NRTL & nrtl keys in datasource
            if "NRTL" in self.datasource.keys():
                # ! NRTL provided
                datasource = self.datasource["NRTL"]
            elif "nrtl" in self.datasource.keys():
                # ! nrtl provided
                datasource = self.datasource["nrtl"]
            elif mixture_ids is not None:
                # ! mixture_ids provided
                # init datasource
                datasource = {}
                # set datasource by mixture ids
                if 'Name' in mixture_ids.keys():
                    # check not empty
                    if mixture_ids.get('Name', None):
                        key_ = mixture_ids['Name']
                        # check key in datasource
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
                elif 'Formula' in mixture_ids.keys():
                    # check not empty
                    if mixture_ids.get('Formula', None):
                        key_ = mixture_ids['Formula']
                        # check key in datasource
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
            else:
                # ! no keys found, use model_input if provided
                # log
                logger.warning(
                    "No NRTL or nrtl key found in datasource, using model_input if provided."
                )
                datasource = {}

            if (
                datasource is not None and
                not isinstance(datasource, dict)
            ):
                raise ValueError(
                    "datasource must be a dictionary."
                )

            # NOTE: set datasource to empty dict if None
            datasource = {} if datasource is None else dict(datasource)

            # NOTE: check model inputs
            # ! when constants are provided in model_input, they override the datasource
            if kwargs.get('model_input') is not None:
                # update the datasource
                datasource.update(kwargs['model_input'])

            # NOTE: final datasource validation
            if datasource is None:
                raise ValueError(
                    "datasource cannot be None."
                )

            if not isinstance(datasource, dict):
                raise ValueError(
                    "datasource must be a dictionary."
                )

            if len(datasource) == 0:
                raise ValueError(
                    "datasource cannot be empty."
                )

            # return
            return datasource
        except Exception as e:
            raise Exception(
                f"Failed to generate NRTL activity data source: {e}"
            ) from e

    # NOTE: extract parameter source
    def _extract_parameter_source(
        self,
        ids: List[str],
        datasource: Dict[str, Any]
    ) -> Optional[Any]:
        try:
            # source
            src: Optional[Dict[str, Any]] = None

            # iterate over ids and extract parameter sources
            for id_ in ids:
                src = datasource.get(id_, None)
                # >> check
                if src is not None:
                    # >> found, break
                    break

            return src
        except Exception as e:
            raise Exception(
                f"Failed to extract parameter source for NRTL activity model: {e}"
            ) from e

    def extract_parameter_sources(
            self,
            datasource: Dict[str, Any],
            include_alpha: bool = True
    ):
        # NOTE: method 1
        # ! Δg_ij, interaction energy parameter
        dg_ij_src = self._extract_parameter_source(
            ids=['dg_ij', 'dg'],
            datasource=datasource
        )

        # NOTE: method 2
        # ! constants a
        a_ij_src = self._extract_parameter_source(
            ids=['a_ij', 'a'],
            datasource=datasource
        )

        # ! constants b
        b_ij_src = self._extract_parameter_source(
            ids=['b_ij', 'b'],
            datasource=datasource
        )

        # ! constants c
        c_ij_src = self._extract_parameter_source(
            ids=['c_ij', 'c'],
            datasource=datasource
        )

        # ! constants d
        d_ij_src = self._extract_parameter_source(
            ids=['d_ij', 'd'],
            datasource=datasource
        )

        # NOTE: α_ij, non-randomness parameter
        # ! check if alpha_ij is provided

        alpha_ij_src = None
        if include_alpha:
            alpha_ij_src = self._extract_parameter_source(
                ids=['alpha_ij', 'alpha'],
                datasource=datasource
            )

        # NOTE: tau_ij, binary interaction parameter
        # ! check if tau_ij is provided
        tau_ij_src = self._extract_parameter_source(
            ids=['tau_ij', 'tau'],
            datasource=datasource
        )

        # res
        return {
            'dg_ij_src': dg_ij_src,
            'a_ij_src': a_ij_src,
            'b_ij_src': b_ij_src,
            'c_ij_src': c_ij_src,
            'd_ij_src': d_ij_src,
            'alpha_ij_src': alpha_ij_src,
            'tau_ij_src': tau_ij_src
        }

    # NOTE: extract parameter values

    def extract_parameter_values(
        self,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        include_alpha: bool = True,
        **kwargs

    ):
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

        # SECTION: extract parameter source
        parameter_sources = self.extract_parameter_sources(
            datasource=datasource,
            include_alpha=include_alpha
        )
        # >> unpack
        dg_ij_src = parameter_sources['dg_ij_src']
        a_ij_src = parameter_sources['a_ij_src']
        b_ij_src = parameter_sources['b_ij_src']
        c_ij_src = parameter_sources['c_ij_src']
        d_ij_src = parameter_sources['d_ij_src']
        alpha_ij_src = parameter_sources['alpha_ij_src']
        tau_ij_src = parameter_sources['tau_ij_src']

        # SECTION: extract data

        # NOTE: >>> check if tau_ij is provided
        if tau_ij_src is not None:
            # ! use tau_ij

            # extract tau_ij values
            tau_ij = self.to_matrix_ij_or(
                data=tau_ij_src,
                property_name='tau',
                symbol_delimiter=symbol_delimiter
            )

            # set tau_ij_cal_method to 0 (no calculation needed)
            tau_ij_cal_method = 0

        elif dg_ij_src is not None:  # >>> check if dg_ij is provided
            # ! use dg_ij

            # extract dg_ij values
            dg_ij = self.to_matrix_ij_or(
                data=dg_ij_src,
                property_name='dg',
                symbol_delimiter=symbol_delimiter
            )

            # set tau_ij_cal_method to 1
            tau_ij_cal_method = 1

        else:
            # >>> check if a_ij, b_ij, c_ij, d_ij are provided

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

            # set tau_ij_cal_method to 2
            tau_ij_cal_method = 2

        # SECTION: extract data
        # NOTE: α_ij, non-randomness parameter
        # check
        if include_alpha and alpha_ij_src is not None:
            # ! extract α_ij values
            alpha_ij = self.to_matrix_ij_or(
                data=alpha_ij_src,
                property_name='alpha',
                symbol_delimiter=symbol_delimiter
            )
        elif include_alpha:
            # ! set alpha_ij to default value 0.3
            alpha_ij_cte = 0.3 * \
                (np.ones((self.comp_num, self.comp_num)) - np.eye(self.comp_num))

            # ! convert alpha_ij_cte to matrix_ij
            alpha_ij = self.to_matrix_ij_or(
                data=alpha_ij_cte,
                property_name='alpha',
                symbol_delimiter=symbol_delimiter
            )

        # NOTE: results
        return {
            'tau_ij': tau_ij,
            'dg_ij': dg_ij,
            'a_ij': a_ij,
            'b_ij': b_ij,
            'c_ij': c_ij,
            'd_ij': d_ij,
            'alpha_ij': alpha_ij,
            'tau_ij_cal_method': tau_ij_cal_method
        }

    # SECTION: temperature validation

    def validate_temperature(
            self,
            temperature: Optional[List[float | str]],
            unit: str = 'K'
    ) -> float:

        # init temperature [K]
        T = -1  # invalid temperature

        # >> check
        if temperature is not None:
            # check if temperature is a list
            if not isinstance(temperature, list):
                raise ValueError(
                    "temperature must be a list of floats or strings."
                )

            # temperature
            T_value = float(temperature[0])
            T_unit = str(temperature[1])

            # convert temperature to Kelvin
            T = pycuc.convert_from_to(
                T_value,
                T_unit,
                unit
            )

        # >> check if temperature is valid
        if T is None or T <= 0:
            raise ValueError(
                "temperature must be provided and greater than 0 K."
            )

        return T
