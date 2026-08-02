# import libs
import time
import logging
from typing import Dict, Literal, List, Optional
from pythermodb_settings.models import Component, Temperature, ComponentKey
from pythermodb_settings.utils import set_component_id, create_mixture_id, measure_time
from pyThermoLinkDB.models import ModelSource
# local
from ..docs import ThermoModelCore
from ..activity import NRTL, UNIQUAC
from ..utils import set_feed_specification
from ..utils.utility import TauCorrelation, map_tau_correlation_to_method

# NOTE: setup logger
logger = logging.getLogger(__name__)


# SECTION: Input configuration
# NOTE: Validate inputs
def _validate_inputs(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
):
    # ! components
    if (
        not isinstance(components, list) or
        not all(isinstance(c, Component) for c in components)
    ):
        raise ValueError(
            "Invalid components input. Must be a list of Component objects."
        )
    if len(components) == 0:
        raise ValueError("Components list is empty.")

    # ! temperature
    if not isinstance(temperature, Temperature):
        raise ValueError(
            "Invalid temperature input. Must be a Temperature object."
        )

    # ! model source
    if not isinstance(model_source, ModelSource):
        raise ValueError(
            "Invalid model_source input. Must be a ModelSource object."
        )


# NOTE: Component input configuration


def _set_component_input_configuration(
    components: List[Component],
    component_key: Literal[
        "Name-State", "Formula-State"
    ],
    mixture_key: Literal[
        "Name", "Formula"
    ],
    separator_symbol: str,
    delimiter: Literal[
        "|",
    ],
    verbose: bool = False,
):
    try:
        # NOTE: component configuration
        # >> component ids
        # ! default using component key
        components_ids = [
            set_component_id(c, component_key, separator_symbol) for c in components
        ]

        # >> name-state
        components_names_state = [
            set_component_id(c, 'Name-State', separator_symbol) for c in components
        ]

        # >> formula-state
        components_formulas_state = [
            set_component_id(c, 'Formula-State', separator_symbol) for c in components
        ]

        # >> component names
        components_names = [
            set_component_id(c, 'Name', separator_symbol) for c in components
        ]

        # >> component formulas
        components_formulas = [
            set_component_id(c, 'Formula', separator_symbol) for c in components
        ]

        # set component ids
        components_ids_dict = {
            "ids": components_ids,
            "Name-State": components_names_state,
            "Formula-State": components_formulas_state,
            "Name": components_names,
            "Formula": components_formulas
        }

        # NOTE: mixture id
        # ! default using mixture key (sorted alphabetically)
        mixture_id = create_mixture_id(
            components=components,
            mixture_key=mixture_key,
            delimiter=delimiter
        )
        # ! by name
        mixture_name = create_mixture_id(
            components=components,
            mixture_key='Name',
            delimiter=delimiter
        )
        # ! by formula
        mixture_formula = create_mixture_id(
            components=components,
            mixture_key='Formula',
            delimiter=delimiter
        )

        # set mixture id
        mixture_ids_dict = {
            "Name": mixture_name,
            "Formula": mixture_formula
        }

        # NOTE: mole fraction
        feed_spec: Dict[str, float] = set_feed_specification(
            components=components,
            component_key="Name"
        )

        # >>> log
        if verbose:
            logger.info("Input preparation successful")

        # res
        return {
            "mole_fraction": feed_spec,
            "mixture_id": mixture_id,
            "components_ids": components_ids,
            "components_ids_dict": components_ids_dict,
            "mixture_ids_dict": mixture_ids_dict
        }

    except Exception as e:
        logger.error(f"Input preparation failed!, {e}")
        raise


# SECTION: Binary interaction parameter (tau) calculators
@measure_time
def calc_tau_ij_using_nrtl_model(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
    tau_correlation: TauCorrelation = "gibbs_energy",
    component_key: Literal[
        "Name-State", "Formula-State"
    ] = "Name-State",
    mixture_key: Literal[
        "Name", "Formula"
    ] = "Name",
    separator_symbol: str = '-',
    delimiter: Literal[
        "|",
    ] = "|",
    message: Optional[str] = None,
    verbose: bool = False,
    output_format: ComponentKey = "Name",
    **kwargs
):
    '''
        Calculate binary interaction parameters (tau_ij) used with NRTL model for a given mixture of components at a specified temperature using the provided model source.

        Parameters
        ----------
        components : List[Component]
            List of Component objects.
        pressure : Pressure
            Pressure object containing pressure value and unit.
        temperature : Temperature
            Temperature object containing temperature value and unit.
        model_source : ModelSource
            datasource and equationsource needed for activity model calculation.
                - datasource: dict, datasource for the component (`generated by PyThermoDB`)
                - equationsource: dict, equationsource for the component (`generated by PyThermoDB`)
        model_name : Literal['NRTL', 'UNIQUAC']
            Name of the activity model to use. Options are 'NRTL' or 'UNIQUAC'.
        tau_correlation : Literal['gibbs_energy', 'extended_temperature', 'inverse_temperature', 'inverse_temperature_squared', 'inverse_log_temperature'], optional
            Correlation method for calculating tau_ij, by default 'M1'.
                - 'gibbs_energy': Calculate tau_ij using dg_ij.
                - 'extended_temperature': Calculate tau_ij from constants a, b, c, d based on the selected correlation.
                - 'inverse_temperature': Calculate tau_ij from constants a, b based on the selected correlation.
                - 'inverse_temperature_squared': Calculate tau_ij from constants a, b, c based on the selected correlation.
                - 'inverse_log_temperature': Calculate tau_ij from constants a, b, c based on the selected correlation.
        component_key : Literal['Name-State', 'Formula-State'], optional
            Key to identify components, by default 'Name-State'.
                - 'Name-State': Component name with state (e.g., 'ethanol-l') is used.
                - 'Formula-State': Component formula with state (e.g., 'C2H5OH-l') is used.
        mixture_key : Literal['Name', 'Formula'], optional
            Key to identify mixture, by default 'Name'.
                - 'Name': Component names (e.g., 'ethanol | butyl-methyl-ether') are used.
                - 'Formula': Component formulas (e.g., 'C2H5OH | C5H12O') are used.
        separator_symbol : str, optional
            Symbol to separate component name and state, by default '-'.
                - e.g., 'ethanol-l' or 'C2H5OH-l'
        delimiter : Literal['|'], optional
            Delimiter to separate components in mixture, by default '|'.
                - e.g., 'ethanol | butyl-methyl-ether' or 'C2H5OH | C5H12O'
        message : Optional[str], optional
            Message for logging, by default None.
                - If None, a default message will be used.
        verbose : bool, optional
            If True, detailed logs will be printed, by default False.
        output_format : Literal['Name', 'Formula'], optional
            Format for output keys, by default 'Name'.
        **kwargs : dict
            Additional keyword arguments.
            - mode : Literal['silent', 'log', 'attach'], optional
                Mode for time measurement logging. Default is 'silent'.
        '''
    try:
        # >> log
        if verbose:
            logger.info(
                f"Calculating tau_ij for components: {[c.name for c in components]} at temperature: {temperature.value} {temperature.unit}"
            )

        # SECTION: validate inputs
        _validate_inputs(
            components=components,
            temperature=temperature,
            model_source=model_source
        )

        # NOTE: message config
        if message is None:
            message = f"Calculating activity coefficient using NRTL model"

        # SECTION: set input configuration
        component_input_config = _set_component_input_configuration(
            components=components,
            component_key=component_key,
            mixture_key=mixture_key,
            separator_symbol=separator_symbol,
            delimiter=delimiter,
            verbose=verbose
        )
        # >> unpack
        mole_fraction = component_input_config["mole_fraction"]
        mixture_id = component_input_config["mixture_id"]
        components_ids = component_input_config["components_ids"]
        components_ids_dict = component_input_config["components_ids_dict"]
        mixture_ids_dict = component_input_config["mixture_ids_dict"]

        # NOTE: model inputs
        model_input = {
            "mole_fraction": mole_fraction,
            "temperature": [temperature.value, temperature.unit],
            "mixture_id": mixture_id,
            "components_ids": components_ids,
            "components_ids_dict": components_ids_dict,
            "mixture_ids_dict": mixture_ids_dict
        }

        # NOTE: model source
        model_source_dict = {
            "datasource": model_source.data_source,
            "equationsource": model_source.equation_source
        }

        # SECTION: initialize activity model
        try:
            # NOTE: thermo manager
            ThermoModelCore_ = ThermoModelCore()
            # NOTE: initialize activity model
            activity_models = ThermoModelCore_.select_activities(
                components=components,
                model_name='NRTL',
                model_source=model_source_dict,
                mixture_ids=mixture_ids_dict,
                components_ids=components_ids_dict,
                **kwargs
            )

            # >>> log
            if verbose:
                logger.info(f"NRTL model initialization successful")
                logger.info(f"Activity model: {activity_models}")
        except Exception as e:
            logger.error(f"Initialization failed!, {e}")
            raise

        # SECTION: calculate activity coefficient
        # NOTE: map tau_correlation to method
        tau_correlation_selected = map_tau_correlation_to_method(
            tau_correlation
        )

        try:
            # NOTE: check nrtl
            if isinstance(activity_models, NRTL):
                # NOTE: set ids
                # ! mixture
                activity_models.mixture_ids = mixture_ids_dict
                # ! components
                activity_models.components_ids = components_ids_dict

                # NOTE: calculate activity coefficient
                res = activity_models.check_and_build_inputs(
                    model_input=model_input,
                    required_keys=['tau_ij', 'alpha_ij'],
                    tau_correlation=tau_correlation_selected,
                    symbol_delimiter=delimiter,
                    return_all=False,
                    **kwargs
                )

                # ! tau_ij
                # ? numpy array
                tau_ij_array = res['tau_ij']
                # ? dict with component keys
                tau_ij_dict = activity_models.to_dict_ij_ext(
                    data=tau_ij_array,
                    components=components,
                    component_key=output_format,  # ? result format
                    symbol_delimiter=delimiter,
                )

                # ! alpha_ij
                # alpha_ij = activity_models.to_dict_ij(
                #     data=res['alpha_ij'],
                #     symbol_delimiter=delimiter,
                # )

                # >>> log
                if verbose:
                    logger.info(f"Tau calculation successful")
                    logger.info(f"Result: {tau_ij_dict}")

                # res
                return tau_ij_array, tau_ij_dict

            else:
                raise ValueError(
                    f"Activity model is not an instance of NRTL. Found: {type(activity_models)}"
                )

        except Exception as e:
            logger.error(f"Tau calculation failed!, {e}")
            raise

    except Exception as e:
        logger.error(f"Error in logging: {e}")
