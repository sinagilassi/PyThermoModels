# import libs
import logging
from typing import List, Literal, cast
# locals
from ..models import (
    PropertyValue,
    ComponentGasFugacityPhaseResult,
    ComponentGasFugacityResult,
    ComponentLiquidFugacityPhaseResult,
    ComponentLiquidFugacityResult,
    MixtureFugacityResult
)

# NOTE: logger
logger = logging.getLogger(__name__)


def parse_gas_fugacity_calc_result(
    res: dict,
    phase_names: List[
        Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']
    ] = ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID'],
) -> ComponentGasFugacityResult:
    '''
    Parse gas fugacity result dictionary into GasFugacityResult model

    Parameters
    ----------
    res : dict
        Result dictionary from fugacity calculation
    phase_names : List[Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']], optional
        List of phase names to consider, by default ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']

    Returns
    -------
    ComponentGasFugacityResult
        Parsed gas fugacity result model
    '''
    try:
        # NOTE: extract phase results
        # parse phase results
        phase_results = {}

        # normalize phase names
        phase_names_normalized = [phase.upper() for phase in phase_names]

        # NOTE: res normalization
        res_normalized = {k.upper(): v for k, v in res.items()}

        # NOTE: phase
        res_phase = res_normalized.get('PHASE', None)
        # >> check
        if res_phase is None or isinstance(res_phase, list) is False:
            logger.warning("Phase information missing or invalid in result")
            res_phase = []

        # NOTE: component
        res_component = res_normalized.get('COMPONENT', None)
        # >> check
        if res_component is None or isinstance(res_component, list) is False:
            logger.warning(
                "Component information missing or invalid in result")
            res_component = []

        # iterate phase names
        for phase_name in phase_names_normalized:
            if phase_name in res_normalized.keys():
                # extract phase data
                phase_data = res_normalized[phase_name]

                # update
                phase_result = ComponentGasFugacityPhaseResult(
                    mole_fraction=phase_data['mole_fraction'],
                    temperature=PropertyValue(
                        value=phase_data['temperature']['value'],
                        unit=phase_data['temperature']['unit'],
                        symbol=phase_data['temperature']['symbol']
                    ),
                    pressure=PropertyValue(
                        value=phase_data['pressure']['value'],
                        unit=phase_data['pressure']['unit'],
                        symbol=phase_data['pressure']['symbol']
                    ),
                    molar_volume=PropertyValue(
                        value=phase_data['molar_volume']['value'],
                        unit=phase_data['molar_volume']['unit'],
                        symbol=phase_data['molar_volume']['symbol']
                    ),
                    compressibility_coefficient=PropertyValue(
                        value=phase_data['compressibility_coefficient']['value'],
                        unit=phase_data['compressibility_coefficient']['unit'],
                        symbol=phase_data['compressibility_coefficient']['symbol']
                    ),
                    fugacity_coefficient=PropertyValue(
                        value=phase_data['fugacity_coefficient']['value'],
                        unit=phase_data['fugacity_coefficient']['unit'],
                        symbol=phase_data['fugacity_coefficient']['symbol']
                    ),
                    fugacity=PropertyValue(
                        value=phase_data['fugacity']['value'],
                        unit=phase_data['fugacity']['unit'],
                        symbol=phase_data['fugacity']['symbol']
                    ),
                    roots=PropertyValue(
                        value=phase_data['roots']['value'],
                        unit=phase_data['roots']['unit'],
                        symbol='Z_roots'  # assuming symbol for roots
                    ),
                    mode=phase_data['mode'],
                    phase=cast(
                        Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID'], phase_name),
                    eos_model=phase_data['eos_model']
                )

                # add to phase results
                phase_results[phase_name] = phase_result

        # create GasFugacityResult
        gas_fugacity_result = ComponentGasFugacityResult(
            phase=res_phase,
            component=res_component,
            results=phase_results
        )

        return gas_fugacity_result
    except Exception as e:
        logger.error(f"Failed to parse gas fugacity result: {e}")
        raise e


def parse_liquid_fugacity_calc_result(
    res: dict,
    phase_names: List[
        Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']
    ] = ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID'],
) -> ComponentLiquidFugacityResult:
    '''
    Parse liquid fugacity result dictionary into LiquidFugacityResult model

    Parameters
    ----------
    res : dict
        Result dictionary from fugacity calculation
    phase_names : List[Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']], optional
        List of phase names to consider, by default ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']

    Returns
    -------
    ComponentLiquidFugacityResult
        Parsed liquid fugacity result model
    '''
    try:
        # NOTE: extract phase results
        # parse phase results
        phase_results = {}

        # normalize phase names
        phase_names_normalized = [phase.upper() for phase in phase_names]

        # NOTE: res normalization
        res_normalized = {k.upper(): v for k, v in res.items()}

        # NOTE: phase
        res_phase = res_normalized.get('PHASE', None)
        # >> check
        if res_phase is None or isinstance(res_phase, list) is False:
            logger.warning("Phase information missing or invalid in result")
            res_phase = []

        # NOTE: component
        res_component = res_normalized.get('COMPONENT', None)
        # >> check
        if res_component is None or isinstance(res_component, list) is False:
            logger.warning(
                "Component information missing or invalid in result")
            res_component = []

        # iterate phase names
        for phase_name in phase_names_normalized:
            if phase_name in res_normalized.keys():
                # extract phase data
                phase_data = res_normalized[phase_name]

                # update
                phase_result = ComponentLiquidFugacityPhaseResult(
                    mole_fraction=phase_data['mole_fraction'],
                    temperature=PropertyValue(
                        value=phase_data['temperature']['value'],
                        unit=phase_data['temperature']['unit'],
                        symbol=phase_data['temperature']['symbol']
                    ),
                    pressure=PropertyValue(
                        value=phase_data['pressure']['value'],
                        unit=phase_data['pressure']['unit'],
                        symbol=phase_data['pressure']['symbol']
                    ),
                    vapor_pressure=PropertyValue(
                        value=phase_data['vapor_pressure']['value'],
                        unit=phase_data['vapor_pressure']['unit'],
                        symbol=phase_data['vapor_pressure']['symbol']
                    ),
                    molar_volume=PropertyValue(
                        value=phase_data['molar_volume']['value'],
                        unit=phase_data['molar_volume']['unit'],
                        symbol=phase_data['molar_volume']['symbol']
                    ),
                    compressibility_coefficient=PropertyValue(
                        value=phase_data['compressibility_coefficient']['value'],
                        unit=phase_data['compressibility_coefficient']['unit'],
                        symbol=phase_data['compressibility_coefficient']['symbol']
                    ),
                    fugacity_coefficient_sat=PropertyValue(
                        value=phase_data['fugacity_coefficient_sat']['value'],
                        unit=phase_data['fugacity_coefficient_sat']['unit'],
                        symbol=phase_data['fugacity_coefficient_sat']['symbol']
                    ),
                    fugacity_coefficient=PropertyValue(
                        value=phase_data['fugacity_coefficient']['value'],
                        unit=phase_data['fugacity_coefficient']['unit'],
                        symbol=phase_data['fugacity_coefficient']['symbol']
                    ),
                    Poynting_term=PropertyValue(
                        value=phase_data['Poynting_term']['value'],
                        unit=phase_data['Poynting_term']['unit'],
                        symbol=phase_data['Poynting_term']['symbol']
                    ),
                    fugacity_sat=PropertyValue(
                        value=phase_data['fugacity_sat']['value'],
                        unit=phase_data['fugacity_sat']['unit'],
                        symbol=phase_data['fugacity_sat']['symbol']
                    ),
                    fugacity=PropertyValue(
                        value=phase_data['fugacity']['value'],
                        unit=phase_data['fugacity']['unit'],
                        symbol=phase_data['fugacity']['symbol']
                    ),
                    roots=PropertyValue(
                        value=phase_data['roots']['value'],
                        unit=phase_data['roots']['unit'],
                        symbol='Z_roots'  # assuming symbol for roots
                    ),
                    mode=phase_data['mode'],
                    phase=cast(
                        Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID'], phase_name),
                    eos_model=phase_data['eos_model']
                )

                # add to phase results
                phase_results[phase_name] = phase_result

        # create LiquidFugacityResult
        liquid_fugacity_result = ComponentLiquidFugacityResult(
            phase=res_phase,
            component=res_component,
            results=phase_results
        )

        return liquid_fugacity_result
    except Exception as e:
        logger.error(f"Failed to parse liquid fugacity result: {e}")
        raise e


def parse_mixture_fugacity_calc_result(
    res: dict,
    phase_names: List[
        Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']
    ] = ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']
) -> MixtureFugacityResult:
    '''
    Parse mixture fugacity result dictionary into MixtureFugacityResult model

    Parameters
    ----------
    res : dict
        Result dictionary from fugacity calculation
    phase_names : List[Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']], optional
        List of phase names to consider, by default ['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID']

    Returns
    -------
    MixtureFugacityResult
        Parsed mixture fugacity result model
    '''
    try:
        # NOTE: extract phase results
        # parse phase results
        phase_results = {}

        # normalize phase names
        phase_names_normalized = [phase.upper() for phase in phase_names]

        # NOTE: res normalization
        res_normalized = {k.upper(): v for k, v in res.items()}

        # NOTE: phase
        res_phase = res_normalized.get('PHASE', None)
        # >> check
        if res_phase is None or isinstance(res_phase, list) is False:
            logger.warning("Phase information missing or invalid in result")
            res_phase = []

        # NOTE: components
        res_components = res_normalized.get('COMPONENT', None)
        # >> check
        if res_components is None or isinstance(res_components, list) is False:
            logger.warning(
                "Components information missing or invalid in result")
            res_components = []

        # iterate phase names
        for phase_name in phase_names_normalized:
            if phase_name in res_normalized.keys():
                # extract phase data
                phase_data = res_normalized[phase_name]

                # update
                phase_result = {}
                for component, comp_data in phase_data.items():
                    phase_result[component] = ComponentGasFugacityPhaseResult(
                        mole_fraction=comp_data['mole_fraction'],
                        temperature=PropertyValue(
                            value=comp_data['temperature']['value'],
                            unit=comp_data['temperature']['unit'],
                            symbol=comp_data['temperature']['symbol']
                        ),
                        pressure=PropertyValue(
                            value=comp_data['pressure']['value'],
                            unit=comp_data['pressure']['unit'],
                            symbol=comp_data['pressure']['symbol']
                        ),
                        molar_volume=PropertyValue(
                            value=comp_data['molar_volume']['value'],
                            unit=comp_data['molar_volume']['unit'],
                            symbol=comp_data['molar_volume']['symbol']
                        ),
                        compressibility_coefficient=PropertyValue(
                            value=comp_data['compressibility_coefficient']['value'],
                            unit=comp_data['compressibility_coefficient']['unit'],
                            symbol=comp_data['compressibility_coefficient']['symbol']
                        ),
                        fugacity_coefficient=PropertyValue(
                            value=comp_data['fugacity_coefficient']['value'],
                            unit=comp_data['fugacity_coefficient']['unit'],
                            symbol=comp_data['fugacity_coefficient']['symbol']
                        ),
                        fugacity=PropertyValue(
                            value=comp_data['fugacity']['value'],
                            unit=comp_data['fugacity']['unit'],
                            symbol=comp_data['fugacity']['symbol']
                        ),
                        roots=PropertyValue(
                            value=comp_data['roots']['value'],
                            unit=comp_data['roots']['unit'],
                            symbol='Z_roots'  # assuming symbol for roots
                        ),
                        mode=comp_data['mode'],
                        phase=cast(
                            Literal['VAPOR', 'LIQUID', 'SUPERCRITICAL', 'VAPOR-LIQUID'], phase_name),
                        eos_model=comp_data['eos_model']
                    )

                # add to phase results
                phase_results[phase_name] = phase_result

        # create MixtureFugacityResult
        mixture_fugacity_result = MixtureFugacityResult(
            components=res_components,
            phase=res_phase,
            results=phase_results
        )

        return mixture_fugacity_result
    except Exception as e:
        logger.error(f"Failed to parse mixture fugacity result: {e}")
        raise e
