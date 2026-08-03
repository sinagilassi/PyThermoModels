# import libs
import logging
from typing import Dict, Literal, List, Optional

from pythermodb_settings.models import Component, Temperature, ComponentKey
from pythermodb_settings.utils import (
    create_mixture_id,
    measure_time,
    set_component_id,
)
from pyThermoLinkDB.models import ModelSource

# local
from ..activity import NRTL, UNIQUAC
from ..docs import ThermoModelCore
from ..utils import set_feed_specification
from ..utils.utility import TauCorrelation, map_tau_correlation_to_method, map_uniquac_tau_correlation_to_method

# NOTE: setup logger
logger = logging.getLogger(__name__)


# SECTION: Input configuration
# NOTE: Validate inputs
def _validate_inputs(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
) -> None:
    """
    Validate shared activity-parameter calculator inputs.

    Parameters
    ----------
    components : List[Component]
        Components included in the liquid-mixture activity model.
    temperature : Temperature
        Operating temperature used for temperature-dependent tau correlations.
    model_source : ModelSource
        Model source containing datasource and equationsource dictionaries.

    Raises
    ------
    ValueError
        If any required input has an unsupported type or empty value.
    """
    # SECTION: components
    if (
        not isinstance(components, list) or
        not all(isinstance(c, Component) for c in components)
    ):
        raise ValueError(
            "Invalid components input. Must be a list of Component objects."
        )
    if len(components) == 0:
        raise ValueError("Components list is empty.")

    # SECTION: temperature
    if not isinstance(temperature, Temperature):
        raise ValueError(
            "Invalid temperature input. Must be a Temperature object."
        )

    # SECTION: model source
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
) -> Dict:
    """
    Build component, mixture, and feed identifiers for activity calculations.

    Parameters
    ----------
    components : List[Component]
        Components included in the activity model.
    component_key : Literal["Name-State", "Formula-State"]
        Component identifier style used for datasource lookup.
    mixture_key : Literal["Name", "Formula"]
        Mixture identifier style used for datasource lookup.
    separator_symbol : str
        Separator between component identity and state.
    delimiter : Literal["|"]
        Delimiter used between components in mixture identifiers.
    verbose : bool
        If True, log input-preparation progress.

    Returns
    -------
    Dict
        Component ids, mixture ids, and default feed specification.
    """
    try:
        # SECTION: component ids
        # ! default using component key
        components_ids = [
            set_component_id(c, component_key, separator_symbol)
            for c in components
        ]

        # NOTE: keep alternate component ids for model-specific fallbacks
        components_names_state = [
            set_component_id(c, 'Name-State', separator_symbol)
            for c in components
        ]
        components_formulas_state = [
            set_component_id(c, 'Formula-State', separator_symbol)
            for c in components
        ]
        components_names = [
            set_component_id(c, 'Name', separator_symbol)
            for c in components
        ]
        components_formulas = [
            set_component_id(c, 'Formula', separator_symbol)
            for c in components
        ]

        components_ids_dict = {
            "ids": components_ids,
            "Name-State": components_names_state,
            "Formula-State": components_formulas_state,
            "Name": components_names,
            "Formula": components_formulas,
        }

        # SECTION: mixture ids
        # ! default using mixture key (sorted alphabetically)
        mixture_id = create_mixture_id(
            components=components,
            mixture_key=mixture_key,
            delimiter=delimiter,
        )
        mixture_name = create_mixture_id(
            components=components,
            mixture_key='Name',
            delimiter=delimiter,
        )
        mixture_formula = create_mixture_id(
            components=components,
            mixture_key='Formula',
            delimiter=delimiter,
        )

        mixture_ids_dict = {
            "Name": mixture_name,
            "Formula": mixture_formula,
        }

        # SECTION: mole fraction placeholder
        feed_spec: Dict[str, float] = set_feed_specification(
            components=components,
            component_key="Name",
        )

        if verbose:
            logger.info("Input preparation successful")

        return {
            "mole_fraction": feed_spec,
            "mixture_id": mixture_id,
            "components_ids": components_ids,
            "components_ids_dict": components_ids_dict,
            "mixture_ids_dict": mixture_ids_dict,
        }

    except Exception as e:
        logger.error(f"Input preparation failed!, {e}")
        raise


# SECTION: Model-specific tau configuration
def _activity_tau_configuration(
    model_name: Literal['NRTL', 'UNIQUAC'],
    tau_correlation: Optional[TauCorrelation],
) -> Dict:
    """
    Return required input keys and internal tau correlation for a model.

    Parameters
    ----------
    model_name : Literal["NRTL", "UNIQUAC"]
        Activity model used to generate tau parameters.
    tau_correlation : Optional[TauCorrelation]
        Public descriptive tau-correlation name. If None, the selected model
        chooses its own default.

    Returns
    -------
    Dict
        Required input keys, model class, and correlation value expected by the
        selected model instance.
    """
    # SECTION: NRTL configuration
    if model_name == 'NRTL':
        # NOTE: default NRTL tau calculation uses Gibbs-energy data
        tau_correlation = (
            "gibbs_energy" if tau_correlation is None else tau_correlation
        )
        return {
            "model_type": NRTL,
            # NOTE: NRTL still expects raw M1-M5 method ids internally
            "tau_correlation": map_tau_correlation_to_method(tau_correlation),
            "generator_kwargs": {
                "include_alpha": False,
            },
        }

    # SECTION: UNIQUAC configuration
    if model_name == 'UNIQUAC':
        # NOTE: default UNIQUAC tau calculation uses extended-temperature data
        tau_correlation = (
            "extended_temperature"
            if tau_correlation is None
            else tau_correlation
        )
        return {
            "model_type": UNIQUAC,
            # NOTE: UNIQUAC maps descriptive names inside its parameter builder
            "tau_correlation": map_uniquac_tau_correlation_to_method(tau_correlation),
            "generator_kwargs": {
                "include_pure_component_parameters": False,
            },
        }

    raise ValueError(
        f"Unsupported activity model '{model_name}'. "
        "Supported models are 'NRTL' and 'UNIQUAC'."
    )


# SECTION: Binary interaction parameter (tau) calculators
@measure_time
def calc_tau_ij(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
    model_name: Literal['NRTL', 'UNIQUAC'],
    tau_correlation: Optional[TauCorrelation] = None,
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
    """
    Calculate tau_ij for NRTL or UNIQUAC using a model source.

    Parameters
    ----------
    components : List[Component]
        Components included in the activity model.
    temperature : Temperature
        Operating temperature used for tau generation.
    model_source : ModelSource
        Source containing activity-model datasource and equationsource.
    model_name : Literal["NRTL", "UNIQUAC"]
        Activity model used to calculate tau_ij.
    tau_correlation : Optional[TauCorrelation]
        Public descriptive tau-correlation name. If None, NRTL defaults to
        `gibbs_energy` and UNIQUAC defaults to `extended_temperature`.
    component_key : Literal["Name-State", "Formula-State"]
        Component identifier style used for datasource lookup.
    mixture_key : Literal["Name", "Formula"]
        Mixture identifier style used for datasource lookup.
    separator_symbol : str
        Separator between component identity and state.
    delimiter : Literal["|"]
        Delimiter used for component-pair keys.
    message : Optional[str]
        Optional logging message.
    verbose : bool
        If True, log detailed progress.
    output_format : ComponentKey
        Component identifier style used for output dictionary keys.
    **kwargs : dict
        Additional keyword arguments forwarded to model initialization and
        parameter generation.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Tuple[np.ndarray, Dict[str, float]]
        Tau matrix ordered by model components and output dictionary keyed by
        requested component format.

    Notes
    -----
    This function generates only `tau_ij`. It intentionally does not extract
    full activity-model parameters such as NRTL `alpha_ij` or UNIQUAC
    `r_i`/`q_i`.
    """
    try:
        # SECTION: validate shared inputs
        _validate_inputs(
            components=components,
            temperature=temperature,
            model_source=model_source,
        )

        # NOTE: model-specific required keys and correlation mapping
        tau_config = _activity_tau_configuration(
            model_name=model_name,
            tau_correlation=tau_correlation,
        )

        if message is None:
            message = f"Calculating tau_ij using {model_name} model"

        if verbose:
            logger.info(message)
            logger.info(
                f"Components: {[c.name for c in components]}, "
                f"temperature: {temperature.value} {temperature.unit}"
            )

        # SECTION: set input configuration
        component_input_config = _set_component_input_configuration(
            components=components,
            component_key=component_key,
            mixture_key=mixture_key,
            separator_symbol=separator_symbol,
            delimiter=delimiter,
            verbose=verbose,
        )
        mole_fraction = component_input_config["mole_fraction"]
        mixture_id = component_input_config["mixture_id"]
        components_ids = component_input_config["components_ids"]
        components_ids_dict = component_input_config["components_ids_dict"]
        mixture_ids_dict = component_input_config["mixture_ids_dict"]

        # SECTION: model input and model source
        model_input = {
            "mole_fraction": mole_fraction,
            "temperature": [temperature.value, temperature.unit],
            "mixture_id": mixture_id,
            "components_ids": components_ids,
            "components_ids_dict": components_ids_dict,
            "mixture_ids_dict": mixture_ids_dict,
        }

        model_source_dict = {
            "datasource": model_source.data_source,
            "equationsource": model_source.equation_source,
        }

        # SECTION: initialize selected activity model
        try:
            ThermoModelCore_ = ThermoModelCore()
            activity_model = ThermoModelCore_.select_activities(
                components=components,
                model_name=model_name,
                model_source=model_source_dict,
                mixture_ids=mixture_ids_dict,
                components_ids=components_ids_dict,
                **kwargs
            )

            if verbose:
                logger.info(f"{model_name} model initialization successful")
                logger.info(f"Activity model: {activity_model}")
        except Exception as e:
            logger.error(f"Initialization failed!, {e}")
            raise

        # SECTION: calculate tau_ij
        # model_type = tau_config["model_type"]
        if not isinstance(activity_model, (NRTL, UNIQUAC)):
            raise TypeError(
                f"Activity model is not an instance of {model_name}. "
                f"Found: {type(activity_model)}"
            )

        # NOTE: set ids explicitly for datasource fallback behavior
        activity_model.mixture_ids = mixture_ids_dict
        activity_model.components_ids = components_ids_dict

        res = activity_model.inputs_generator(
            temperature=[temperature.value, temperature.unit],
            tau_correlation=tau_config["tau_correlation"],
            symbol_delimiter=delimiter,
            mixture_ids=mixture_ids_dict,
            components_ids=components_ids_dict,
            model_input=model_input,
            **tau_config["generator_kwargs"],
            **kwargs
        )

        tau_ij_array = res['tau_ij']
        tau_ij_dict = activity_model.to_dict_ij_ext(
            data=tau_ij_array,
            components=components,
            component_key=output_format,
            symbol_delimiter=delimiter,
        )

        if verbose:
            logger.info("Tau calculation successful")
            logger.info(f"Result: {tau_ij_dict}")

        return tau_ij_array, tau_ij_dict

    except Exception as e:
        logger.error(f"Tau calculation failed!, {e}")
        raise

# ! NRTL tau_ij calculator


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
    """
    Calculate NRTL binary interaction parameters (`tau_ij`).

    Parameters
    ----------
    components : List[Component]
        Components included in the liquid-mixture activity model.
    temperature : Temperature
        Operating temperature used to evaluate the NRTL tau correlation.
    model_source : ModelSource
        Model source containing the activity-model datasource and
        equationsource dictionaries.
    tau_correlation : TauCorrelation
        Correlation used to generate `tau_ij`. Defaults to `gibbs_energy`. The tau correlation can be one of the following:
            - `gibbs_energy`
            - `extended_temperature`
            - `inverse_temperature`
            - `inverse_temperature_squared`
            - `inverse_log_temperature`
    component_key : Literal["Name-State", "Formula-State"]
        Component identifier style used for datasource lookup.
    mixture_key : Literal["Name", "Formula"]
        Mixture identifier style used for datasource lookup.
    separator_symbol : str
        Separator between component identity and state.
    delimiter : Literal["|"]
        Delimiter used between component ids in binary interaction keys.
    message : Optional[str]
        Optional logging message.
    verbose : bool
        If True, log detailed calculation progress.
    output_format : ComponentKey
        Component identifier style used for output dictionary keys.
    **kwargs : dict
        Additional keyword arguments forwarded to model initialization and
        parameter generation.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Tuple[np.ndarray, Dict[str, float]]
        Tau matrix ordered by model components and output dictionary keyed by
        the requested component format.

    Notes
    -----
    This tau-only helper does not extract `alpha_ij`. It requires either
    direct `tau_ij`, `dg_ij`, or coefficient matrices compatible with
    `tau_correlation`.
    """
    # NOTE: thin compatibility wrapper around the shared dispatcher
    return calc_tau_ij(
        components=components,
        temperature=temperature,
        model_source=model_source,
        model_name='NRTL',
        tau_correlation=tau_correlation,
        component_key=component_key,
        mixture_key=mixture_key,
        separator_symbol=separator_symbol,
        delimiter=delimiter,
        message=message,
        verbose=verbose,
        output_format=output_format,
        **kwargs
    )

# ! UNIQUAC tau_ij calculator


def calc_tau_ij_using_uniquac_model(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
    tau_correlation: TauCorrelation = "extended_temperature",
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
    """
    Calculate UNIQUAC binary interaction parameters (`tau_ij`).

    Parameters
    ----------
    components : List[Component]
        Components included in the liquid-mixture activity model.
    temperature : Temperature
        Operating temperature used to evaluate the UNIQUAC tau correlation.
    model_source : ModelSource
        Model source containing the activity-model datasource and
        equationsource dictionaries.
    tau_correlation : TauCorrelation
        Correlation used to generate `tau_ij`. Defaults to `extended_temperature`. The tau correlation can be one of the following:
            - `gibbs_energy`
            - `extended_temperature`
            - `inverse_temperature`
            - `inverse_temperature_squared`
            - `inverse_log_temperature`
    component_key : Literal["Name-State", "Formula-State"]
        Component identifier style used for datasource lookup.
    mixture_key : Literal["Name", "Formula"]
        Mixture identifier style used for datasource lookup.
    separator_symbol : str
        Separator between component identity and state.
    delimiter : Literal["|"]
        Delimiter used between component ids in binary interaction keys.
    message : Optional[str]
        Optional logging message.
    verbose : bool
        If True, log detailed calculation progress.
    output_format : ComponentKey
        Component identifier style used for output dictionary keys.
    **kwargs : dict
        Additional keyword arguments forwarded to model initialization and
        parameter generation.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Tuple[np.ndarray, Dict[str, float]]
        Tau matrix ordered by model components and output dictionary keyed by
        the requested component format.

    Notes
    -----
    This tau-only helper does not extract `r_i` or `q_i`. It requires either
    direct `tau_ij`, `dU_ij`, or coefficient matrices compatible with
    `tau_correlation`. The default UNIQUAC correlation is
    `extended_temperature`.
    """
    # NOTE: thin model-specific wrapper around the shared dispatcher
    return calc_tau_ij(
        components=components,
        temperature=temperature,
        model_source=model_source,
        model_name='UNIQUAC',
        tau_correlation=tau_correlation,
        component_key=component_key,
        mixture_key=mixture_key,
        separator_symbol=separator_symbol,
        delimiter=delimiter,
        message=message,
        verbose=verbose,
        output_format=output_format,
        **kwargs
    )
