# import libs
import logging
from typing import Literal
from pythermodb_settings.models import Component
# local
from ..models import ComponentIdentity

# NOTE: logger
logger = logging.getLogger(__name__)


def create_component_id(
    component: Component,
    separator_symbol: str = '-'
) -> ComponentIdentity:
    '''
    Create component name-state and formula-state identifiers.

    Parameters
    ----------
    component : Component
        The component for which to create the identifiers.
    separator_symbol : str, optional
        The symbol to use as a separator between the name/formula and

    Returns
    -------
    ComponentIdentity
        The component identity containing name-state and formula-state
        identifiers.
    '''
    try:
        # NOTE: extract component name
        component_name = component.name.strip()
        component_formula = component.formula.strip()
        component_state = component.state.strip().lower()

        # >> separator
        separator_symbol = separator_symbol.strip()

        # SECTION: create component identifiers
        name_state = f"{component_name}{separator_symbol}{component_state}"
        formula_state = f"{component_formula}{separator_symbol}{component_state}"

        return ComponentIdentity(
            name_state=name_state,
            formula_state=formula_state
        )
    except Exception as e:
        logger.error(
            f"Failed to create component identifiers for "
            f"'{component}': {e}"
        )
        raise e


def set_component_id(
    component: Component,
    component_key: Literal[
        'Name-State', 'Formula-State'
    ],
    separator_symbol: str = '-',
) -> str:
    '''
    Set component identifier based on the specified key.

    Parameters
    ----------
    component : Component
        The component for which to set the identifier.
    component_key : str
        The key to determine which identifier to use.
        Options are:
            - 'Name-State': Use the name-state identifier.
            - 'Formula-State': Use the formula-state identifier.
    separator_symbol : str, optional
        The symbol to use as a separator between the name/formula and state.
        Default is '-'.

    Returns
    -------
    str
        The component identifier based on the specified key.
    '''
    try:
        # NOTE: create component id
        component_idx: ComponentIdentity = create_component_id(
            component=component,
            separator_symbol=separator_symbol
        )

        # set component id
        if component_key == "Name-State":
            return component_idx.name_state
        elif component_key == "Formula-State":
            return component_idx.formula_state
        else:
            raise ValueError(
                f"Invalid component_key '{component_key}'. "
                f"Must be 'Name-State' or 'Formula-State'."
            )
    except Exception as e:
        logger.error(
            f"Failed to set component identifier for "
            f"'{component}': {e}"
        )
        raise e
