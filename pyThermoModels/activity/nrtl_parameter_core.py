# import libs
import logging
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
    ) -> Optional[Dict[str, Any]]:
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
            datasource: Dict[str, Any]
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

    # SECTION: temperature validation
    def _validate_temperature(
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
