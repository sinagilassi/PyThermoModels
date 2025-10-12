# import libs
import logging
from typing import Dict, List, Literal, Tuple, Optional
from pythermodb_settings.models import Component, Temperature
import numpy as np
from pyThermoDB.core import TableMatrixData
import pycuc
from pythermodb_settings.utils import set_component_id
# local
from .nrtl import NRTL
from .uniquac import UNIQUAC

# NOTE: logger
logger = logging.getLogger(__name__)


def calc_dg_ij_using_nrtl_model(
    components: List[Component],
    temperature: Temperature,
    a_ij: Dict[str, float] | TableMatrixData,
    b_ij: Dict[str, float] | TableMatrixData,
    c_ij: Dict[str, float] | TableMatrixData,
    mixture_delimiter: Literal[
        "|", "_"
    ] = "|",
    component_key: Literal[
        'Name',
        'Formula',
        'Name-State',
        'Formula-State',
        "Name-Formula-State",
        "Formula-Name-State"
    ] = 'Name',
    component_delimiter: str = '-',
    verbose: bool = False,
    message: Optional[str] = None
) -> Tuple[
    np.ndarray,
    Dict[str, float],
    Dict[str, float]
]:
    """
    Calculate interaction energy parameter `dg_ij` matrix dependent of temperature using `NRTL` model.

    Parameters
    ----------
    components : List[Component]
        List of Component instances.
    temperature : Temperature
        Temperature value and unit.
    a_ij : Dict[str, float] | TableMatrixData
        Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
    b_ij : Dict[str, float] | TableMatrixData
        Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
    c_ij : Dict[str, float] | TableMatrixData
        Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
    mixture_delimiter : Literal["|", "_"]
        Delimiter for the mixture id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    component_delimiter : str
        Separator symbol used in component id when combining multiple keys. Default is '-' such as CO2-g.
    verbose : bool
        If True, print detailed logs. Default is False.
    message : Optional[str]
        Optional message to include in logs. Default is None.

    Returns
    -------
    dg_ij : np.ndarray
        Interaction energy parameter `dg_ij` matrix for NRTL model with respect to component ids sorted alphabetically.
    dg_ij_comp : dict
        Dictionary of interaction energy parameters where keys are component pairs and values are their respective dg_ij values.
    dg_ij_comp_upd : dict
        Dictionary of interaction energy parameters with updated component ids based on `component_key` and `component_delimiter`.

    Notes
    -----
    1. The interaction energy parameter matrix is calculated using the formula:

        dg_ij = a_ij + b_ij * T + c_ij * T^2

    where T is the temperature [K].

    2. All parameters including a_ij, b_ij, c_ij must be in the same format (numpy array, dict or TableMatrixData).
    """
    try:
        # SECTION: inputs config
        # ! components
        if not isinstance(components, list) or not all(isinstance(comp, Component) for comp in components):
            raise ValueError(
                "components must be a list of Component instances.")

        if len(components) < 2:
            raise ValueError(
                "At least two components are required to calculate interaction parameters.")

        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError("temperature must be an instance of Temperature.")

        # ! a_ij, b_ij, c_ij
        if not (
            isinstance(a_ij, (dict, TableMatrixData)) and
            isinstance(b_ij, (dict, TableMatrixData)) and
            isinstance(c_ij, (dict, TableMatrixData))
        ):
            raise ValueError(
                "a_ij, b_ij, and c_ij must all be of the same type: either dict or TableMatrixData."
            )

        # SECTION: model config
        # NOTE: components ids
        # ! by default: component name
        components_ids = [comp.name.strip() for comp in components]

        # NOTE: nrtl model instance
        nrtl = NRTL(
            components=components_ids
        )

        # SECTION: calculation
        # NOTE: temperature conversion
        # temperature [K]
        temperature_value_K = pycuc.convert_from_to(
            value=temperature.value,
            from_unit=temperature.unit,
            to_unit="K"
        )

        # NOTE: calculate dg_ij
        dg_ij, dg_ij_comp = nrtl.cal_dg_ij_M1(
            temperature=temperature_value_K,
            a_ij=a_ij,
            b_ij=b_ij,
            c_ij=c_ij,
            symbol_delimiter=mixture_delimiter,
        )

        # NOTE: sort dg_ij matrix based on component ids
        # components
        comp_names = nrtl.components
        # component id
        comp_idx = nrtl.comp_idx

        # init
        dg_ij_comp_upd = {}

        # iterate over dg_ij_comp keys
        for k, v in dg_ij_comp.items():
            # split key into components
            comp_i, comp_j = k.split(mixture_delimiter)

            # generate new component ids based on component_key
            comp_i = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_i.strip().lower()
            ][0]

            comp_j = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_j.strip().lower()
            ][0]

            # new key
            new_key = f"{comp_i}{mixture_delimiter}{comp_j}"

            # >> update dictionary
            dg_ij_comp_upd[new_key] = float(v)

        # NOTE: optional message
        if message is None:
            message = f"Calculated dg_ij matrix for {'-'.join(components_ids)} using NRTL model at {temperature_value_K} K."

        if verbose:
            logger.info(message)

        return dg_ij, dg_ij_comp, dg_ij_comp_upd
    except Exception as e:
        raise Exception(f"Error in calculating using NRTL model: {str(e)}")


def calc_tau_ij_with_dg_ij_using_nrtl_model(
    components: List[Component],
    temperature: Temperature,
    dg_ij: Dict[str, float] | TableMatrixData,
    mixture_delimiter: Literal[
        "|", "_"
    ] = "|",
    component_key: Literal[
        'Name',
        'Formula',
        'Name-State',
        'Formula-State',
        "Name-Formula-State",
        "Formula-Name-State"
    ] = 'Name',
    component_delimiter: str = '-',
    verbose: bool = False,
    message: Optional[str] = None
) -> Tuple[
    np.ndarray,
    Dict[str, float],
    Dict[str, float]
]:
    """
    Calculate interaction parameter `tau_ij` matrix dependent of temperature with `dg_ij` using `NRTL` model.

    Parameters
    ----------
    components : List[Component]
        List of Component instances.
    temperature : Temperature
        Temperature value and unit.
    dg_ij : Dict[str, float] | TableMatrixData
        Interaction energy parameter dg_ij matrix where dg_ij[i][j] between component i and j.
    mixture_delimiter : Literal["|", "_"]
        Delimiter for the mixture id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    component_delimiter : str
        Separator symbol used in component id when combining multiple keys. Default is '-' such as CO2-g.
    verbose : bool
        If True, print detailed logs. Default is False.
    message : Optional[str]
        Optional message to include in logs. Default is None.

    Returns
    -------
    tau_ij : np.ndarray
        Interaction parameter `tau_ij` matrix for NRTL model with respect to component ids sorted alphabetically.
    tau_ij_comp : dict
        Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.
    tau_ij_comp_upd : dict
        Dictionary of interaction parameters with updated component ids based on `component_key` and `component_delimiter`.
    """
    try:
        # SECTION: inputs config
        # ! components
        if not isinstance(components, list) or not all(isinstance(comp, Component) for comp in components):
            raise ValueError(
                "components must be a list of Component instances.")

        if len(components) < 2:
            raise ValueError(
                "At least two components are required to calculate interaction parameters.")

        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError("temperature must be an instance of Temperature.")

        # ! dg_ij
        if not (isinstance(dg_ij, (dict, TableMatrixData))):
            raise ValueError(
                "dg_ij must be either dict or TableMatrixData."
            )

        # SECTION: model config
        # NOTE: components ids
        # ! by default: component name
        components_ids = [comp.name.strip() for comp in components]

        # NOTE: nrtl model instance
        nrtl = NRTL(
            components=components_ids
        )

        # SECTION: calculation
        # NOTE: temperature conversion
        # temperature [K]
        temperature_value_K = pycuc.convert_from_to(
            value=temperature.value,
            from_unit=temperature.unit,
            to_unit="K"
        )

        # NOTE: calculate tau_ij
        tau_ij, tau_ij_comp = nrtl.cal_tau_ij_M1(
            temperature=temperature_value_K,
            dg_ij=dg_ij,
            symbol_delimiter=mixture_delimiter,
        )

        # NOTE: sort tau_ij matrix based on component ids
        # components
        comp_names = nrtl.components
        # component id
        comp_idx = nrtl.comp_idx

        # init
        tau_ij_comp_upd = {}

        # iterate over tau_ij_comp keys
        for k, v in tau_ij_comp.items():
            # split key into components
            comp_i, comp_j = k.split(mixture_delimiter)

            # generate new component ids based on component_key
            comp_i = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_i.strip().lower()
            ][0]

            comp_j = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_j.strip().lower()
            ][0]

            # new key
            new_key = f"{comp_i}{mixture_delimiter}{comp_j}"

            # >> update dictionary
            tau_ij_comp_upd[new_key] = float(v)

        # NOTE: optional message
        if message is None:
            message = f"Calculated tau_ij matrix for {'-'.join(components_ids)} using NRTL model at {temperature_value_K} K."

        if verbose:
            logger.info(message)

        return tau_ij, tau_ij_comp, tau_ij_comp_upd
    except Exception as e:
        raise Exception(f"Error in calculating using NRTL model: {str(e)}")


def calc_dU_ij_using_uniquac_model(
    components: List[Component],
    temperature: Temperature,
    a_ij: Dict[str, float] | TableMatrixData,
    b_ij: Dict[str, float] | TableMatrixData,
    c_ij: Dict[str, float] | TableMatrixData,
    mixture_delimiter: Literal[
        "|", "_"
    ] = "|",
    component_key: Literal[
        'Name',
        'Formula',
        'Name-State',
        'Formula-State',
        "Name-Formula-State",
        "Formula-Name-State"
    ] = 'Name',
    component_delimiter: str = '-',
    verbose: bool = False,
    message: Optional[str] = None
) -> Tuple[
    np.ndarray,
    Dict[str, float],
    Dict[str, float]
]:
    """
    Calculate interaction energy parameter `dU_ij` matrix dependent of temperature. using `UNIQUAC` model.

    Parameters
    ----------
    components : List[Component]
        List of Component instances.
    temperature : Temperature
        Temperature value and unit.
    a_ij : Dict[str, float] | TableMatrixData
        Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
    b_ij : Dict[str, float] | TableMatrixData
        Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
    c_ij : Dict[str, float] | TableMatrixData
        Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
    mixture_delimiter : Literal["|", "_"]
        Delimiter for the mixture id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    component_delimiter : str
        Separator symbol used in component id when combining multiple keys. Default is '-' such as CO2-g.
    verbose : bool
        If True, print detailed logs. Default is False.
    message : Optional[str]
        Optional message to include in logs. Default is None.

    Returns
    -------
    dU_ij : np.ndarray
        Interaction energy parameter `dU_ij` matrix for UNIQUAC model.
    dU_ij_comp : dict
        Dictionary of interaction energy parameters where keys are component pairs and values are their respective dU_ij values.
    dU_ij_comp_upd : dict
        Dictionary of interaction energy parameters with updated component ids based on `component_key` and `component_delimiter`.

    Notes
    -----
    1. The interaction energy parameter matrix is calculated using the formula:

        dU_ij = a_ij + b_ij * T + c_ij * T^2

        where T is the temperature [K].

    2. All parameters including a_ij, b_ij, c_ij must be in the same format (dict or TableMatrixData).
    """
    try:
        # SECTION: inputs config
        # ! components
        if not isinstance(components, list) or not all(isinstance(comp, Component) for comp in components):
            raise ValueError(
                "components must be a list of Component instances.")

        if len(components) < 2:
            raise ValueError(
                "At least two components are required to calculate interaction parameters.")

        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError("temperature must be an instance of Temperature.")

        # ! a_ij, b_ij, c_ij
        if not (
            isinstance(a_ij, (dict, TableMatrixData)) and
            isinstance(b_ij, (dict, TableMatrixData)) and
            isinstance(c_ij, (dict, TableMatrixData))
        ):
            raise ValueError(
                "a_ij, b_ij, and c_ij must all be of the same type: either dict or TableMatrixData."
            )

        # SECTION: model config
        # NOTE: components ids
        # ! by default: component name
        components_ids = [comp.name.strip() for comp in components]

        # NOTE: uniquac model instance
        uniquac = UNIQUAC(
            components=components_ids
        )

        # SECTION: calculation
        # NOTE: temperature conversion
        # temperature [K]
        temperature_value_K = pycuc.convert_from_to(
            value=temperature.value,
            from_unit=temperature.unit,
            to_unit="K"
        )

        # NOTE: calculate dU_ij
        dU_ij, dU_ij_comp = uniquac.cal_dU_ij_M1(
            temperature=temperature_value_K,
            a_ij=a_ij,
            b_ij=b_ij,
            c_ij=c_ij,
            symbol_delimiter=mixture_delimiter,
        )

        # NOTE: sort dU_ij matrix based on component ids
        # components
        comp_names = uniquac.components
        # component id
        comp_idx = uniquac.comp_idx

        # init
        dU_ij_comp_upd = {}

        # iterate over dU_ij_comp keys
        for k, v in dU_ij_comp.items():
            # split key into components
            comp_i, comp_j = k.split(mixture_delimiter)

            # generate new component ids based on component_key
            comp_i = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_i.strip().lower()
            ][0]

            comp_j = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_j.strip().lower()
            ][0]

            # new key
            new_key = f"{comp_i}{mixture_delimiter}{comp_j}"

            # >> update dictionary
            dU_ij_comp_upd[new_key] = float(v)

        # NOTE: optional message
        if message is None:
            message = f"Calculated dU_ij matrix for {'-'.join(components_ids)} using UNIQUAC model at {temperature_value_K} K."

        if verbose:
            logger.info(message)

        return dU_ij, dU_ij_comp, dU_ij_comp_upd
    except Exception as e:
        raise Exception(f"Error in calculating using UNIQUAC model: {str(e)}")


def calc_tau_ij_with_dU_ij_using_uniquac_model(
    components: List[Component],
    temperature: Temperature,
    dU_ij: Dict[str, float] | TableMatrixData,
    mixture_delimiter: Literal[
        "|", "_"
    ] = "|",
    component_key: Literal[
        'Name',
        'Formula',
        'Name-State',
        'Formula-State',
        "Name-Formula-State",
        "Formula-Name-State"
    ] = 'Name',
    component_delimiter: str = '-',
    verbose: bool = False,
    message: Optional[str] = None
) -> Tuple[
    np.ndarray,
    Dict[str, float],
    Dict[str, float]
]:
    """
    Calculate interaction parameter `tau_ij` matrix dependent of temperature using `UNIQUAC` model.

    Parameters
    ----------
    components : List[Component]
        List of Component instances.
    temperature : Temperature
        Temperature value and unit.
    dU_ij : Dict[str, float] | TableMatrixData
        Interaction energy parameter dU_ij matrix where dU_ij[i][j] between component i and j.
    mixture_delimiter : Literal["|", "_"]
        Delimiter for the mixture id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    component_delimiter : str
        Separator symbol used in component id when combining multiple keys. Default is '-' such as CO2
    verbose : bool
        If True, print detailed logs. Default is False.
    message : Optional[str]
        Optional message to include in logs. Default is None.

    Returns
    -------
    tau_ij : np.ndarray
        Interaction parameter `tau_ij` matrix for UNIQUAC model with respect to component ids sorted alphabetically.
    tau_ij_comp : dict
        Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.
    tau_ij_comp_upd : dict
        Dictionary of interaction parameters with updated component ids based on `component_key` and `component_delimiter`.

    Notes
    -----
    1. The tau_ij matrix is calculated using the formula:

    tau_ij = exp(-dU_ij / (R * T))

    where R is the universal gas constant [J/mol/K] and T is the temperature [K].

    2. dU_ij must be in the format of either dict or TableMatrixData.
    """
    try:
        # SECTION: inputs config
        # ! components
        if not isinstance(components, list) or not all(isinstance(comp, Component) for comp in components):
            raise ValueError(
                "components must be a list of Component instances.")

        if len(components) < 2:
            raise ValueError(
                "At least two components are required to calculate interaction parameters.")

        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError("temperature must be an instance of Temperature.")

        # ! dU_ij
        if not (isinstance(dU_ij, (dict, TableMatrixData))):
            raise ValueError(
                "dU_ij must be either dict or TableMatrixData."
            )

        # SECTION: model config
        # NOTE: components ids
        # ! by default: component name
        components_ids = [comp.name.strip() for comp in components]

        # NOTE: uniquac model instance
        uniquac = UNIQUAC(
            components=components_ids
        )

        # SECTION: calculation
        # NOTE: temperature conversion
        # temperature [K]
        temperature_value_K = pycuc.convert_from_to(
            value=temperature.value,
            from_unit=temperature.unit,
            to_unit="K"
        )

        # NOTE: calculate tau_ij
        tau_ij, tau_ij_comp = uniquac.cal_tau_ij_M1(
            temperature=temperature_value_K,
            dU_ij=dU_ij,
            symbol_delimiter=mixture_delimiter,
        )

        # NOTE: sort tau_ij matrix based on component ids
        # components
        comp_names = uniquac.components
        # component id
        comp_idx = uniquac.comp_idx

        # init
        tau_ij_comp_upd = {}

        # iterate over tau_ij_comp keys
        for k, v in tau_ij_comp.items():
            # split key into components
            comp_i, comp_j = k.split(mixture_delimiter)

            # generate new component ids based on component_key
            comp_i = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_i.strip().lower()
            ][0]

            comp_j = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_j.strip().lower()
            ][0]

            # new key
            new_key = f"{comp_i}{mixture_delimiter}{comp_j}"

            # >> update dictionary
            tau_ij_comp_upd[new_key] = float(v)

        # NOTE: optional message
        if message is None:
            message = f"Calculated tau_ij matrix for {'-'.join(components_ids)} using UNIQUAC model at {temperature_value_K} K."

        if verbose:
            logger.info(message)

        return tau_ij, tau_ij_comp, tau_ij_comp_upd
    except Exception as e:
        raise Exception(f"Error in calculating using UNIQUAC model: {str(e)}")


def calc_tau_ij(
    components: List[Component],
    temperature: Temperature,
    a_ij: Dict[str, float] | TableMatrixData,
    b_ij: Dict[str, float] | TableMatrixData,
    c_ij: Dict[str, float] | TableMatrixData,
    d_ij: Dict[str, float] | TableMatrixData,
    model: Literal['NRTL', 'UNIQUAC'] = 'NRTL',
    mixture_delimiter: Literal[
        "|", "_"
    ] = "|",
    component_key: Literal[
        'Name',
        'Formula',
        'Name-State',
        'Formula-State',
        "Name-Formula-State",
        "Formula-Name-State"
    ] = 'Name',
    component_delimiter: str = '-',
    verbose: bool = False,
    message: Optional[str] = None
):
    """
    Calculate interaction parameter `tau_ij` matrix dependent of temperature using specified activity model.

    Parameters
    ----------
    components : List[Component]
        List of Component instances.
    temperature : Temperature
        Temperature value and unit.
    a_ij : Dict[str, float] | TableMatrixData
        Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
    b_ij : Dict[str, float] | TableMatrixData
        Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
    c_ij : Dict[str, float] | TableMatrixData
        Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
    d_ij : Dict[str, float] | TableMatrixData
        Interaction parameter d_ij matrix where d_ij[i][j] between component i and j.
    model : Literal['NRTL', 'UNIQUAC']
        Activity model to use for calculation. Default is 'NRTL'.
    mixture_delimiter : Literal["|", "_"]
        Delimiter for the mixture id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    component_delimiter : str
        Separator symbol used in component id when combining multiple keys. Default is '-' such as CO2-g.
    verbose : bool
        If True, print detailed logs. Default is False.
    message : Optional[str]
        Optional message to include in logs. Default is None.

    Returns
    -------
    tau_ij : np.ndarray
        Interaction parameter `tau_ij` matrix for specified activity model with respect to component ids sorted alphabetically.
    tau_ij_comp : dict
        Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.
    tau_ij_comp_upd : dict
        Dictionary of interaction parameters with updated component ids based on `component_key` and `component_delimiter`.

    Raises
    ------
    ValueError
        If the specified model is not supported.

    Notes
    -----
    1. For both 'NRTL' and 'UNIQUAC' models, the interaction parameter matrix is calculated using the formula:

        tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T

    where T is the temperature [K].

    2. All parameters including a_ij, b_ij, c_ij, and d_ij must be in the same format (dict or TableMatrixData).
    3. The universal gas constant R = 8.314 J/mol/K is used in the UNIQUAC model for tau_ij calculation.
    """
    try:
        # SECTION: inputs config
        # ! components
        if not isinstance(components, list) or not all(isinstance(comp, Component) for comp in components):
            raise ValueError(
                "components must be a list of Component instances.")

        if len(components) < 2:
            raise ValueError(
                "At least two components are required to calculate interaction parameters.")

        # ! temperature
        if not isinstance(temperature, Temperature):
            raise ValueError("temperature must be an instance of Temperature.")

        # ! a_ij, b_ij, c_ij, d_ij
        if not (
            isinstance(a_ij, (dict, TableMatrixData)) and
            isinstance(b_ij, (dict, TableMatrixData)) and
            isinstance(c_ij, (dict, TableMatrixData)) and
            isinstance(d_ij, (dict, TableMatrixData))
        ):
            raise ValueError(
                "a_ij, b_ij, and c_ij must all be of the same type: either dict or TableMatrixData."
            )

        # SECTION: model config
        # NOTE: temperature conversion
        # temperature [K]
        temperature_value_K = pycuc.convert_from_to(
            value=temperature.value,
            from_unit=temperature.unit,
            to_unit="K"
        )

        # NOTE: components ids
        # ! by default: component name
        components_ids = [comp.name.strip() for comp in components]

        # SECTION: calculation
        if model == 'NRTL':
            # NOTE: init NRTL model instance
            nrtl = NRTL(
                components=components_ids
            )

            # NOTE: calculate tau_ij
            tau_ij, tau_ij_comp = nrtl.cal_tau_ij_M2(
                temperature=temperature_value_K,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                d_ij=d_ij,
                symbol_delimiter=mixture_delimiter,
            )

        elif model == 'UNIQUAC':
            # NOTE: init UNIQUAC model instance
            uniquac = UNIQUAC(
                components=components_ids
            )

            # NOTE: calculate tau_ij
            tau_ij, tau_ij_comp = uniquac.cal_tau_ij_M2(
                temperature=temperature_value_K,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                d_ij=d_ij,
                symbol_delimiter=mixture_delimiter,
            )

        else:
            raise ValueError(
                f"Unsupported model '{model}'. Supported models are 'NRTL' and 'UNIQUAC'."
            )

        # NOTE: sort tau_ij matrix based on component ids
        # init
        tau_ij_comp_upd = {}

        # iterate over tau_ij_comp keys
        for k, v in tau_ij_comp.items():
            # split key into components
            comp_i, comp_j = k.split(mixture_delimiter)

            # generate new component ids based on component_key
            comp_i = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_i.strip().lower()
            ][0]

            comp_j = [
                set_component_id(
                    component=comp,
                    component_key=component_key,
                    separator_symbol=component_delimiter
                ) for comp in components if comp.name.strip().lower() == comp_j.strip().lower()
            ][0]

            # new key
            new_key = f"{comp_i}{mixture_delimiter}{comp_j}"

            # >> update dictionary
            tau_ij_comp_upd[new_key] = float(v)

        # NOTE: optional message
        if message is None:
            message = f"Calculated tau_ij matrix for {'-'.join(components_ids)} using {model} model at {temperature_value_K} K."

        if verbose:
            logger.info(message)

        return tau_ij, tau_ij_comp, tau_ij_comp_upd
    except Exception as e:
        raise Exception(f"Error in calculating using {model} model: {str(e)}")
