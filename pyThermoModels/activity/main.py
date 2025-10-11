# import libs
import logging
from typing import Dict, List, Literal
from pythermodb_settings.models import Component, Temperature
import numpy as np
from pyThermoDB.core import TableMatrixData
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
    symbol_delimiter: Literal[
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
    verbose: bool = False
):
    """
    Calculate interaction energy parameter `dg_ij` matrix dependent of temperature using `NRTL` model.

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin [K].
    a_ij : Dict[str, float] | TableMatrixData
        Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
    b_ij : Dict[str, float] | TableMatrixData
        Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
    c_ij : Dict[str, float] | TableMatrixData
        Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
    symbol_delimiter : Literal["|", "_"]
        Delimiter for the component id. Default is "|".
    component_key : Literal['Name', 'Formula', 'Name-State', 'Formula-State', "Name-Formula-State", "Formula-Name-State"]
        Key to identify components. Default is 'Name'.
    verbose : bool
        If True, print detailed logs. Default is False.

    Returns
    -------
    dg_ij : np.ndarray
        Interaction energy parameter `dg_ij` matrix for NRTL model with respect to component ids sorted alphabetically.
    dg_ij_comp : dict
        Dictionary of interaction energy parameters where keys are component pairs and values are their respective dg_ij values.

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
        if not (isinstance(a_ij, (dict, TableMatrixData)) and isinstance(b_ij, (dict, TableMatrixData)) and isinstance(c_ij, (dict, TableMatrixData))):
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

    except Exception as e:
        raise Exception(f"Error in calculating using NRTL model: {str(e)}")
