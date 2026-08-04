# import libs
import time
import logging
from typing import Dict, Literal, List, Optional, Sequence, TypeAlias, cast
from pythermodb_settings.models import Component, Temperature, Pressure
from pythermodb_settings.utils import set_component_id, create_mixture_id, measure_time
from pyThermoLinkDB.models import ModelSource
# local
from ..docs import ThermoModelCore
from ..activity import NRTL, UNIQUAC
from ..utils import set_feed_specification
from ..utils.utility import TauCorrelation

# NOTE: logger
logger = logging.getLogger(__name__)

# SECTION: activity coefficient calculation

TauSourcePriorityItem: TypeAlias = Literal[
    "tau",
    "energy",
    "dg",
    "dU",
    "coefficients",
]


def _activity_tau_correlation_default(
    model_name: str,
    datasource: Dict,
    mixture_ids: Optional[Dict[str, str]] = None,
    tau_source_priority: Sequence[TauSourcePriorityItem] = (
        "tau",
        "energy",
        "coefficients",
    ),
) -> TauCorrelation:
    """
    Select the default tau correlation for high-level activity calculations.
    """
    energy_keys = {
        "NRTL": ("dg_ij", "dg"),
        "UNIQUAC": ("dU_ij", "dU"),
    }
    direct_tau_keys = ("tau_ij", "tau")
    default_by_model: Dict[str, TauCorrelation] = {
        "NRTL": "gibbs_energy",
        "UNIQUAC": "extended_temperature",
    }

    if model_name not in default_by_model:
        return "gibbs_energy"

    source_aliases = {
        "tau": "tau",
        "energy": "energy",
        "dg": "energy",
        "dU": "energy",
        "coefficients": "coefficients",
    }

    if (
        not isinstance(tau_source_priority, Sequence) or
        isinstance(tau_source_priority, (str, bytes))
    ):
        raise TypeError("tau_source_priority must be a list or tuple.")

    normalized_priority = []
    for source in tau_source_priority:
        if source not in source_aliases:
            raise ValueError(
                "Unsupported tau source priority item "
                f"{source!r}. Expected one of: {', '.join(source_aliases)}."
            )
        normalized_source = cast(
            Literal["tau", "energy", "coefficients"],
            source_aliases[source],
        )
        if normalized_source in normalized_priority:
            raise ValueError(
                f"Duplicate tau source priority item: {source!r}."
            )
        normalized_priority.append(normalized_source)

    if len(normalized_priority) == 0:
        raise ValueError("tau_source_priority cannot be empty.")

    datasource_ = {}
    if isinstance(datasource, dict):
        if model_name in datasource:
            datasource_ = datasource.get(model_name) or {}
        elif model_name.lower() in datasource:
            datasource_ = datasource.get(model_name.lower()) or {}
        elif mixture_ids:
            for key in ("Name", "Formula"):
                mixture_id = mixture_ids.get(key)
                if mixture_id in datasource:
                    datasource_ = datasource.get(mixture_id) or {}
                    break

    if not isinstance(datasource_, dict):
        return default_by_model[model_name]

    has_direct_tau = any(
        datasource_.get(key) is not None and datasource_.get(key) != 'None'
        for key in direct_tau_keys
    )
    has_energy_source = any(
        datasource_.get(key) is not None and datasource_.get(key) != 'None'
        for key in energy_keys[model_name]
    )
    has_coefficients = any(
        datasource_.get(key) is not None and datasource_.get(key) != 'None'
        for key in ("a_ij", "a", "b_ij", "b", "c_ij", "c", "d_ij", "d")
    )

    source_to_correlation: Dict[
        Literal["tau", "energy", "coefficients"],
        tuple[TauCorrelation, bool],
    ] = {
        "tau": ("direct_tau", has_direct_tau),
        "energy": ("gibbs_energy", has_energy_source),
        "coefficients": ("extended_temperature", has_coefficients),
    }
    for source in normalized_priority:
        correlation, is_available = source_to_correlation[source]
        if is_available:
            return correlation

    return default_by_model[model_name]


@measure_time
def calc_activity_coefficient(
    components: List[Component],
    pressure: Pressure,
    temperature: Temperature,
    model_source: ModelSource,
    model_name: Literal['NRTL', 'UNIQUAC'],
    tau_correlation: Optional[TauCorrelation] = None,
    tau_source_priority: Sequence[TauSourcePriorityItem] = (
        "tau",
        "energy",
        "coefficients",
    ),
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
    **kwargs
):
    '''
    Calculate activity coefficient using `NRTL` or `UNIQUAC` models.

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
    tau_correlation : Literal['direct_tau', 'gibbs_energy', 'extended_temperature', 'inverse_temperature', 'inverse_temperature_squared', 'inverse_log_temperature'], optional
        Correlation method for calculating tau_ij. If omitted, the source is
        inspected and direct tau tables are preferred over generated tau.
            - 'direct_tau': Use tau_ij directly from the source table.
            - 'gibbs_energy': Calculate tau_ij using dg_ij or dU_ij.
            - 'extended_temperature': Calculate tau_ij from constants a, b, c, d based on the selected correlation.
            - 'inverse_temperature': Calculate tau_ij from constants a, b based on the selected correlation.
            - 'inverse_temperature_squared': Calculate tau_ij from constants a, b, c based on the selected correlation.
            - 'inverse_log_temperature': Calculate tau_ij from constants a, b, c based on the selected correlation.
    tau_source_priority : Sequence[Literal['tau', 'energy', 'dg', 'dU', 'coefficients']], optional
        Ordered source preference used only when tau_correlation is omitted.
        Default is ('tau', 'energy', 'coefficients'). The 'dg' and 'dU'
        entries are accepted as aliases for 'energy'.
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
    **kwargs : dict
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    res: Dict[str, float | Dict]
        Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients as:
        - property_name: str
            Name of the property calculated.
        - components: List[str]
            List of component names.
        - mole_fraction: Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        - value: Dict[str, float]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
        - unit: float
            Unit of the property calculated.
        - symbol: str
            Symbol of the property calculated.
        - message: str
            Message to be displayed.
    other_values: Dict[str, float | Dict]
        Dictionary of other values used for the calculation.
    G_ex: Dict[str, float | Dict]
        Dictionary of excess Gibbs free energy.

    Notes
    -----
    For `NRTL` model, the other_values will contain:
        - AcCo_i_comp: Dict[str, float]
        - tau_ij: np.ndarray
        - tau_ij_comp: Dict[str, float]
        - alpha_ij: np.ndarray
        - alpha_ij_comp: Dict[str, float]
        - G_ij: np.ndarray
        - G_ij_comp: Dict[str, float]
        - calculation_mode: str

    For `UNIQUAC` model, the other_values will contain:
        - AcCo_i_comp: Dict[str, float]
        - tau_ij: np.ndarray
        - tau_ij_comp: Dict[str, float]
        - r_i: np.ndarray
        - r_i_comp: Dict[str, float]
        - q_i: np.ndarray
        - q_i_comp: Dict[str, float]
        - calculation_mode: str

    The G_ex will contain:
        - property_name: str
        - components: List[str]
        - mole_fraction: List[float]
        - value: float
        - unit: str
        - symbol: str
        - message: str

    Supported correlations
    ----------------------
    For `NRTL` model, the following correlations are supported for calculating tau_ij:

    direct_tau : M0
        Use a supplied tau_ij table directly.

    gibbs_energy : M1
        Gibbs-energy-based correlation:

        ``tau_ij = dg_ij / (R * T)``

    extended_temperature : M2
        Extended temperature-dependent correlation:

        ``tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T``

    inverse_temperature : M3
        Linear inverse-temperature correlation:

        ``tau_ij = a_ij + b_ij / T``

    inverse_temperature_squared : M4
        Second-order inverse-temperature correlation:

        ``tau_ij = a_ij + b_ij / T + c_ij / T**2``

    inverse_log_temperature : M5
        Inverse- and logarithmic-temperature correlation:

        ``tau_ij = a_ij + b_ij / T + c_ij * ln(T)``
    '''
    try:
        # LINK: start time
        start_time = time.time()

        # >>> log
        if verbose:
            logger.info(
                f"Calculating activity coefficient using {model_name} model"
            )

        # SECTION: validate inputs
        # ! components
        if (
            not isinstance(components, list) or
            not all(isinstance(c, Component) for c in components)
        ):
            raise ValueError(
                "Invalid components input. Must be a list of Component objects.")
        if len(components) == 0:
            raise ValueError("Components list is empty.")
        # ! pressure
        if not isinstance(pressure, Pressure):
            raise ValueError(
                "Invalid pressure input. Must be a Pressure object.")
        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError(
                "Invalid temperature input. Must be a Temperature object.")
        # ! model source
        if not isinstance(model_source, ModelSource):
            raise ValueError(
                "Invalid model_source input. Must be a ModelSource object.")

        # NOTE: message config
        if message is None:
            message = f"Calculating activity coefficient using {model_name} model"

        # SECTION: input preparation
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

            # NOTE: model source
            model_source_dict = {
                "datasource": model_source.data_source,
                "equationsource": model_source.equation_source
            }

            # NOTE: model input
            model_input = {
                "mole_fraction": feed_spec,
                "pressure": [pressure.value, pressure.unit],
                "temperature": [temperature.value, temperature.unit],
                "mixture_id": mixture_id,
                "components_ids": components_ids,
                "components_ids_dict": components_ids_dict,
                "mixture_ids_dict": mixture_ids_dict
            }

            # >>> log
            if verbose:
                logger.info("Input preparation successful")
                logger.info(f"Model input: {model_input}")

        except Exception as e:
            logger.error(f"Input preparation failed!, {e}")
            raise

        # SECTION: initialize activity model
        try:
            # NOTE: thermo manager
            ThermoModelCore_ = ThermoModelCore()
            # NOTE: initialize activity model
            activity_models = ThermoModelCore_.select_activities(
                components=components,
                model_name=model_name,
                model_source=model_source_dict,
                mixture_ids=mixture_ids_dict,
                components_ids=components_ids_dict,
                **kwargs
            )

            # >>> log
            if verbose:
                logger.info(f"{model_name} model initialization successful")
                logger.info(f"Activity model: {activity_models}")
        except Exception as e:
            logger.error(f"Initialization failed!, {e}")
            raise

        if tau_correlation is None:
            tau_correlation = _activity_tau_correlation_default(
                model_name=model_name,
                datasource=model_source.data_source,
                mixture_ids=mixture_ids_dict,
                tau_source_priority=tau_source_priority,
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
                res, others = activity_models.cal(
                    model_input=model_input,
                    tau_correlation=tau_correlation,
                    message=message,
                    **kwargs
                )

                # NOTE: calculate excess Gibbs energy
                G_ex = activity_models.excess_gibbs_free_energy()
            elif isinstance(activity_models, UNIQUAC):
                # NOTE: set ids
                # ! mixture
                activity_models.mixture_ids = mixture_ids_dict
                # ! components
                activity_models.components_ids = components_ids_dict

                # NOTE: calculate activity coefficient
                res, others = activity_models.cal(
                    model_input=model_input,
                    tau_correlation=tau_correlation,
                    message=message,
                    **kwargs
                )

                # NOTE: calculate excess Gibbs energy
                G_ex = activity_models.excess_gibbs_free_energy()
            else:
                raise TypeError(
                    f"activity_models is not `NRTL` or `UNIQUAC`, but {type(activity_models)}")

            # LINK: end time
            end_time = time.time()
            elapsed_time = end_time - start_time

            # >>> log
            if verbose:
                logger.info(
                    f"Activity coefficient calculation successful, elapsed time: {elapsed_time:.2f} seconds")

            # return
            return res, others, G_ex
        except Exception as e:
            logger.error(f"calculation failed!, {e}")
            raise
    except Exception as e:
        raise Exception("calculation failed!, ", e)


def calc_activity_coefficient_using_nrtl_model(
    components: List[Component],
    pressure: Pressure,
    temperature: Temperature,
    tau_ij: Dict[str, float | int],
    alpha_ij: Dict[str, float | int],
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
    **kwargs
):
    '''
    Calculate activity coefficient using `NRTL` model when tau_ij and alpha_ij are provided.

    Parameters
    ----------
    components : List[Component]
        List of Component objects.
    pressure : Pressure
        Pressure object containing pressure value and unit.
    temperature : Temperature
        Temperature object containing temperature value and unit.
    tau_ij : Dict[str, float | int]
        Interaction parameter matrix (tau_ij) as a dictionary where keys are in the format 'component1|component2' and values are the interaction parameters.
    alpha_ij : Dict[str, float | int]
        Non-randomness parameter matrix (alpha_ij) as a dictionary where keys are in the format 'component1|component2' and values are the non-randomness parameters.
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
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    res: Dict[str, float | Dict]
        Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients as:
        - property_name: str
            Name of the property calculated.
        - components: List[str]
            List of component names.
        - mole_fraction: Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        - value: Dict[str, float]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
        - unit: float
            Unit of the property calculated.
        - symbol: str
            Symbol of the property calculated.
        - message: str
            Message to be displayed.
    other_values: Dict[str, float | Dict]
        Dictionary of other values used for the calculation.
    G_ex: Dict[str, float | Dict]
        Dictionary of excess Gibbs free energy.

    Notes
    -----
    For `NRTL` model, the other_values will contain:
        - AcCo_i_comp: Dict[str, float]
        - tau_ij: np.ndarray
        - tau_ij_comp: Dict[str, float]
        - alpha_ij: np.ndarray
        - alpha_ij_comp: Dict[str, float]
        - G_ij: np.ndarray
        - G_ij_comp: Dict[str, float]
        - calculation_mode: str

    The G_ex will contain:
        - property_name: str
        - components: List[str]
        - mole_fraction: List[float]
        - value: float
        - unit: str
        - symbol: str
        - message: str
    '''
    try:
        # LINK: start time
        start_time = time.time()

        # NOTE: model name
        model_name = 'NRTL'

        # >>> log
        if verbose:
            logger.info(
                f"Calculating activity coefficient using {model_name} model"
            )

        # SECTION: validate inputs
        # ! components
        if (
            not isinstance(components, list) or
            not all(isinstance(c, Component) for c in components)
        ):
            raise ValueError(
                "Invalid components input. Must be a list of Component objects.")
        if len(components) == 0:
            raise ValueError("Components list is empty.")
        # ! pressure
        if not isinstance(pressure, Pressure):
            raise ValueError(
                "Invalid pressure input. Must be a Pressure object.")
        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError(
                "Invalid temperature input. Must be a Temperature object.")

        # ! tau_ij
        if (
            not isinstance(tau_ij, dict) or
            not all(isinstance(k, str) and isinstance(v, (int, float))
                    for k, v in tau_ij.items())
        ):
            raise ValueError(
                "Invalid tau_ij input. Must be a dictionary with string keys and numeric values.")

        if len(tau_ij) == 0:
            raise ValueError("tau_ij dictionary is empty.")

        # ! alpha_ij
        if (
            not isinstance(alpha_ij, dict) or
            not all(isinstance(k, str) and isinstance(v, (int, float)) for
                    k, v in alpha_ij.items())
        ):
            raise ValueError(
                "Invalid alpha_ij input. Must be a dictionary with string keys and numeric values.")

        if len(alpha_ij) == 0:
            raise ValueError("alpha_ij dictionary is empty.")

        # NOTE: message config
        if message is None:
            message = f"Calculating activity coefficient using {model_name} model"

        # SECTION: input preparation
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
            # ! default using mixture key
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

            # NOTE: model input
            model_input = {
                "mole_fraction": feed_spec,
                "pressure": [pressure.value, pressure.unit],
                "temperature": [temperature.value, temperature.unit],
                "mixture_id": mixture_id,
                "components_ids": components_ids,
                "components_ids_dict": components_ids_dict,
                "mixture_ids_dict": mixture_ids_dict,
                "tau_ij": tau_ij,
                "alpha_ij": alpha_ij
            }

            # >>> log
            if verbose:
                logger.info("Input preparation successful")
                logger.info(f"Model input: {model_input}")

        except Exception as e:
            logger.error(f"Input preparation failed!, {e}")
            raise

        # SECTION: initialize activity model
        try:
            # NOTE: thermo manager
            ThermoModelCore_ = ThermoModelCore()
            # NOTE: initialize activity model
            activity_models = ThermoModelCore_.select_activities(
                components=components,
                model_name=model_name,
                model_source=None,
                mixture_ids=mixture_ids_dict,
                components_ids=components_ids_dict,
                **kwargs
            )

            # >>> log
            if verbose:
                logger.info(f"{model_name} model initialization successful")
                logger.info(f"Activity model: {activity_models}")
        except Exception as e:
            logger.error(f"Initialization failed!, {e}")
            raise

        # SECTION: calculate activity coefficient
        try:
            # NOTE: check nrtl
            if isinstance(activity_models, NRTL):
                # NOTE: set ids
                # ! mixture
                activity_models.mixture_ids = mixture_ids_dict
                # ! components
                activity_models.components_ids = components_ids_dict

                # NOTE: calculate activity coefficient
                res, others = activity_models.cal(
                    model_input=model_input,
                    message=message,
                    **kwargs
                )

                # NOTE: calculate excess Gibbs energy
                G_ex = activity_models.excess_gibbs_free_energy()
            else:
                raise TypeError(
                    f"activity_models is not `NRTL`, but {type(activity_models)}")

            # LINK: end time
            end_time = time.time()
            elapsed_time = end_time - start_time

            # >>> log
            if verbose:
                logger.info(
                    f"Activity coefficient calculation successful, elapsed time: {elapsed_time:.2f} seconds")

            # return
            return res, others, G_ex
        except Exception as e:
            logger.error(f"calculation failed!, {e}")
            raise
    except Exception as e:
        raise Exception("calculation failed!, ", e)


def calc_activity_coefficient_using_uniquac_model(
    components: List[Component],
    pressure: Pressure,
    temperature: Temperature,
    tau_ij: Dict[str, float | int],
    r_i: Dict[str, float | int],
    q_i: Dict[str, float | int],
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
    **kwargs
):
    '''
    Calculate activity coefficient using `UNIQUAC` model when tau_ij, r_i, and q_i are provided.

    Parameters
    ----------
    components : List[Component]
        List of Component objects.
    pressure : Pressure
        Pressure object containing pressure value and unit.
    temperature : Temperature
        Temperature object containing temperature value and unit.
    tau_ij : Dict[str, float | int]
        Interaction parameter matrix (tau_ij) as a dictionary where keys are in the format 'component1|component2' and values are the interaction parameters.
    r_i : Dict[str, float | int]
        Volume parameter (r_i) as a dictionary where keys are component names and values are the volume parameters.
    q_i : Dict[str, float | int]
        Surface area parameter (q_i) as a dictionary where keys are component names and values are the surface area parameters.
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
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    res: Dict[str, float | Dict]
        Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients as:
        - property_name: str
            Name of the property calculated.
        - components: List[str]
            List of component names.
        - mole_fraction: Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        - value: Dict[str, float]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
        - unit: float
            Unit of the property calculated.
        - symbol: str
            Symbol of the property calculated.
        - message: str
            Message to be displayed.
    other_values: Dict[str, float | Dict]
        Dictionary of other values used for the calculation.
    G_ex: Dict[str, float | Dict]
        Dictionary of excess Gibbs free energy.

    Notes
    -----
    For `UNIQUAC` model, the other_values will contain:
        - AcCo_i_comp: Dict[str, float]
        - tau_ij: np.ndarray
        - tau_ij_comp: Dict[str, float]
        - r_i: np.ndarray
        - r_i_comp: Dict[str, float]
        - q_i: np.ndarray
        - q_i_comp: Dict[str, float]
        - calculation_mode: str

    The G_ex will contain:
        - property_name: str
        - components: List[str]
        - mole_fraction: List[float]
        - value: float
        - unit: str
        - symbol: str
        - message: str
    '''
    try:
        # LINK: start time
        start_time = time.time()

        # NOTE: model name
        model_name = 'UNIQUAC'

        # >>> log
        if verbose:
            logger.info(
                f"Calculating activity coefficient using {model_name} model"
            )

        # SECTION: validate inputs
        # ! components
        if (
            not isinstance(components, list) or
            not all(isinstance(c, Component) for c in components)
        ):
            raise ValueError(
                "Invalid components input. Must be a list of Component objects.")
        if len(components) == 0:
            raise ValueError("Components list is empty.")
        # ! pressure
        if not isinstance(pressure, Pressure):
            raise ValueError(
                "Invalid pressure input. Must be a Pressure object.")
        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError(
                "Invalid temperature input. Must be a Temperature object.")

        # ! tau_ij
        if (
            not isinstance(tau_ij, dict) or
            not all(isinstance(k, str) and isinstance(v, (int, float))
                    for k, v in tau_ij.items())
        ):
            raise ValueError(
                "Invalid tau_ij input. Must be a dictionary with string keys and numeric values.")

        if len(tau_ij) == 0:
            raise ValueError("tau_ij dictionary is empty.")

        # ! r_i
        if (
            not isinstance(r_i, dict) or
            not all(isinstance(k, str) and isinstance(v, (int, float))
                    for k, v in r_i.items())
        ):
            raise ValueError(
                "Invalid r_i input. Must be a dictionary with string keys and numeric values.")

        if len(r_i) == 0:
            raise ValueError("r_i dictionary is empty.")

        # ! q_i
        if (
            not isinstance(q_i, dict) or
            not all(isinstance(k, str) and isinstance(v, (int, float))
                    for k, v in q_i.items())
        ):
            raise ValueError(
                "Invalid q_i input. Must be a dictionary with string keys and numeric values.")

        if len(q_i) == 0:
            raise ValueError("q_i dictionary is empty.")

        # NOTE: message config
        if message is None:
            message = f"Calculating activity coefficient using {model_name} model"

        # SECTION: input preparation
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
            # ! default using mixture key
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

            # NOTE: model input
            model_input = {
                "mole_fraction": feed_spec,
                "pressure": [pressure.value, pressure.unit],
                "temperature": [temperature.value, temperature.unit],
                "mixture_id": mixture_id,
                "components_ids": components_ids,
                "components_ids_dict": components_ids_dict,
                "mixture_ids_dict": mixture_ids_dict,
                "tau_ij": tau_ij,
                "r_i": r_i,
                "q_i": q_i
            }

            # >>> log
            if verbose:
                logger.info("Input preparation successful")
                logger.info(f"Model input: {model_input}")

        except Exception as e:
            logger.error(f"Input preparation failed!, {e}")
            raise

        # SECTION: initialize activity model
        try:
            # NOTE: thermo manager
            ThermoModelCore_ = ThermoModelCore()
            # NOTE: initialize activity model
            activity_models = ThermoModelCore_.select_activities(
                components=components,
                model_name=model_name,
                model_source=None,
                mixture_ids=mixture_ids_dict,
                components_ids=components_ids_dict,
                **kwargs
            )

            # >>> log
            if verbose:
                logger.info(f"{model_name} model initialization successful")
                logger.info(f"Activity model: {activity_models}")
        except Exception as e:
            logger.error(f"Initialization failed!, {e}")
            raise

        # SECTION: calculate activity coefficient
        try:
            # NOTE: check uniquac
            if isinstance(activity_models, UNIQUAC):
                # NOTE: set ids
                # ! mixture
                activity_models.mixture_ids = mixture_ids_dict
                # ! components
                activity_models.components_ids = components_ids_dict

                # NOTE: calculate activity coefficient
                res, others = activity_models.cal(
                    model_input=model_input,
                    message=message,
                    **kwargs
                )

                # NOTE: calculate excess Gibbs energy
                G_ex = activity_models.excess_gibbs_free_energy()
            else:
                raise TypeError(
                    f"activity_models is not `UNIQUAC`, but {type(activity_models)}")

            # LINK: end time
            end_time = time.time()
            elapsed_time = end_time - start_time

            # >>> log
            if verbose:
                logger.info(
                    f"Activity coefficient calculation successful, elapsed time: {elapsed_time:.2f} seconds")

            # return
            return res, others, G_ex
        except Exception as e:
            logger.error(f"calculation failed!, {e}")
            raise
    except Exception as e:
        raise Exception("calculation failed!, ", e)
