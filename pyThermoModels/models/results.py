# import libs
from pydantic import BaseModel, Field
from typing import Dict, List, Union, Literal

# SECTION: gas fugacity results model for single component
# NOTE: The component might be in different phases as vapor, liquid, vapor-liquid, supercritical, etc.
# {
#     'mole_fraction': 1.0,
#     'temperature': {'value': 300.1, 'unit': 'K', 'symbol': 'T'},
#     'pressure': {'value': 999000.0, 'unit': 'Pa', 'symbol': 'P'},
#     'molar_volume': {'value': 0.0020494767707729513, 'unit': 'm3/mol', 'symbol': 'MoVo'},
#     'compressibility_coefficient': {'value': 0.8205552301471566, 'unit': 'dimensionless', 'symbol': 'Z'},
#     'fugacity_coefficient': {'value': 0.8493394106824814, 'unit': 'dimensionless', 'symbol': 'phi'},
#     'fugacity': {'value': 848490.0712717989, 'unit': 'Pa', 'symbol': 'Fug_PURE'},
#     'roots': {'value': [0.8205552301471566], 'unit': 'dimensionless'},
#     'mode': 'SINGLE',
#     'phase': 'VAPOR',
#     'eos_model': 'PR'
# }

# NOTE: fugacity result


class PropertyValue(BaseModel):
    '''Model for a property with value, unit, and symbol'''
    value: Union[float, int, List[float]]
    unit: str
    symbol: str


class ComponentGasFugacityPhaseResult(BaseModel):
    '''
    Phase fugacity result model for fugacity calculations
    '''
    mole_fraction: float = Field(
        ...,
        description="Mole fraction of the component"
    )
    temperature: PropertyValue = Field(
        ..., description="Temperature with value, unit, and symbol"
    )
    pressure: PropertyValue = Field(
        ...,
        description="Pressure with value, unit, and symbol"
    )
    molar_volume: PropertyValue = Field(
        ...,
        description="Molar volume with value, unit, and symbol"
    )
    compressibility_coefficient: PropertyValue = Field(
        ...,
        description="Compressibility coefficient (Z-factor)"
    )
    fugacity_coefficient: PropertyValue = Field(
        ..., description="Fugacity coefficient (phi)")
    fugacity: PropertyValue = Field(
        ...,
        description="Fugacity value with unit and symbol"
    )
    roots: PropertyValue = Field(
        ...,
        description="Roots of the equation of state")
    mode: Literal[
        "SINGLE",
        "MIXTURE"
    ] = Field(..., description="Calculation mode")
    phase: Literal[
        "VAPOR",
        "LIQUID",
        "SUPERCRITICAL",
        "VAPOR-LIQUID"
    ] = Field(..., description="Phase of the substance")
    eos_model: str = Field(
        ...,
        description="Equation of state model used (e.g., 'PR', 'SRK')"
    )

    class Config:
        '''Pydantic model configuration'''
        json_schema_extra = {
            "example": {
                'mole_fraction': 1.0,
                'temperature': {'value': 300.1, 'unit': 'K', 'symbol': 'T'},
                'pressure': {'value': 999000.0, 'unit': 'Pa', 'symbol': 'P'},
                'molar_volume': {'value': 0.0020494767707729513, 'unit': 'm3/mol', 'symbol': 'MoVo'},
                'compressibility_coefficient': {'value': 0.8205552301471566, 'unit': 'dimensionless', 'symbol': 'Z'},
                'fugacity_coefficient': {'value': 0.8493394106824814, 'unit': 'dimensionless', 'symbol': 'phi'},
                'fugacity': {'value': 848490.0712717989, 'unit': 'Pa', 'symbol': 'Fug_PURE'},
                'roots': {'value': [0.8205552301471566], 'unit': 'dimensionless'},
                'mode': 'SINGLE',
                'phase': 'VAPOR',
                'eos_model': 'PR'
            }
        }


class ComponentGasFugacityResult(BaseModel):
    '''
    '''
    phase: list[str] = Field(
        ...,
        description="List of phases considered in the calculation, e.g., ['vapor', 'liquid']"
    )
    component: list[str] = Field(
        ...,
        description="Component name for which fugacity is calculated"
    )
    results: Dict[str, ComponentGasFugacityPhaseResult] = Field(
        ...,
        description="Dictionary mapping phase names to their respective GasFugacityPhaseResult models"
    )

# SECTION: liquid fugacity results model for single component
# NOTE: The component is assumed to be in liquid phase only for this model
# {
#     'liquid': {
#         'mole_fraction': 1.0,
#         'temperature': {'value': 300.1, 'unit': 'K', 'symbol': 'T'},
#         'pressure': {'value': 1000000.0, 'unit': 'Pa', 'symbol': 'P'},
#         'vapor_pressure': {'value': 999723.1044, 'unit': 'Pa', 'symbol': 'VaPr'},
#         'molar_volume': {'value': 8.847861260184554e-05, 'unit': 'm3/mol', 'symbol': 'MoVo'},
#         'compressibility_coefficient': {'value': 0.03545009167303335, 'unit': 'dimensionless', 'symbol': 'Z'},
#         'fugacity_coefficient_sat': {'value': 0.9638260398645724, 'unit': 'dimensionless', 'symbol': 'phi_SAT'},
#         'fugacity_coefficient': {'value': 0.9635686216131526, 'unit': 'dimensionless', 'symbol': 'phi'},
#         'Poynting_term': {'value': 1.0000098187413604, 'unit': 'dimensionless', 'symbol': 'Poynting'},
#         'fugacity_sat': {'value': 963559.1606749684, 'unit': 'Pa', 'symbol': 'Fug_SAT'},
#         'fugacity': {'value': 963568.6216131526, 'unit': 'Pa', 'symbol': 'Fug_PURE'},
#         'roots': {'value': [0.03545009167303335, 0.8204021644184413], 'unit': 'dimensionless', 'symbol': 'Z_i'},
#         'mode': 'SINGLE',
#         'phase': 'LIQUID',
#         'eos_model': 'PR'
#     },
#     'phase': ['liquid'],
#     'component': ['propane']
# }


class ComponentLiquidFugacityPhaseResult(BaseModel):
    '''
    Phase fugacity result model for liquid fugacity calculations
    '''
    mole_fraction: float = Field(
        ...,
        description="Mole fraction of the component"
    )
    temperature: PropertyValue = Field(
        ..., description="Temperature with value, unit, and symbol"
    )
    pressure: PropertyValue = Field(
        ...,
        description="Pressure with value, unit, and symbol"
    )
    vapor_pressure: PropertyValue = Field(
        ...,
        description="Vapor pressure with value, unit, and symbol"
    )
    molar_volume: PropertyValue = Field(
        ...,
        description="Molar volume with value, unit, and symbol"
    )
    compressibility_coefficient: PropertyValue = Field(
        ...,
        description="Compressibility coefficient (Z-factor)"
    )
    fugacity_coefficient_sat: PropertyValue = Field(
        ..., description="Fugacity coefficient at saturation (phi_SAT)")
    fugacity_coefficient: PropertyValue = Field(
        ..., description="Fugacity coefficient (phi)")
    Poynting_term: PropertyValue = Field(
        ...,
        description="Poynting correction term"
    )
    fugacity_sat: PropertyValue = Field(
        ...,
        description="Fugacity at saturation with unit and symbol"
    )
    fugacity: PropertyValue = Field(
        ...,
        description="Fugacity value with unit and symbol"
    )
    roots: PropertyValue = Field(
        ...,
        description="Roots of the equation of state")
    mode: Literal[
        "SINGLE",
        "MIXTURE"
    ] = Field(..., description="Calculation mode")
    phase: Literal[
        "LIQUID",
        "VAPOR",
        "SUPERCRITICAL",
        "VAPOR-LIQUID"
    ] = Field(..., description="Phase of the substance")
    eos_model: str = Field(
        ...,
        description="Equation of state model used (e.g., 'PR', 'SRK')"
    )

    class Config:
        '''Pydantic model configuration'''
        json_schema_extra = {
            "example": {
                'mole_fraction': 1.0,
                'temperature': {'value': 300.1, 'unit': 'K', 'symbol': 'T'},
                'pressure': {'value': 1000000.0, 'unit': 'Pa', 'symbol': 'P'},
                'vapor_pressure': {'value': 999723.1044, 'unit': 'Pa', 'symbol': 'VaPr'},
                'molar_volume': {'value': 8.847861260184554e-05, 'unit': 'm3/mol', 'symbol': 'MoVo'},
                'compressibility_coefficient': {'value': 0.03545009167303335, 'unit': 'dimensionless', 'symbol': 'Z'},
                'fugacity_coefficient_sat': {'value': 0.9638260398645724, 'unit': 'dimensionless', 'symbol': 'phi_SAT'},
                'fugacity_coefficient': {'value': 0.9635686216131526, 'unit': 'dimensionless', 'symbol': 'phi'},
                'Poynting_term': {'value': 1.0000098187413604, 'unit': 'dimensionless', 'symbol': 'Poynting'},
                'fugacity_sat': {'value': 963559.1606749684, 'unit': 'Pa', 'symbol': 'Fug_SAT'},
                'fugacity': {'value': 963568.6216131526, 'unit': 'Pa', 'symbol': 'Fug_PURE'},
                'roots': {'value': [0.03545009167303335, 0.8204021644184413], 'unit': 'dimensionless', 'symbol': 'Z_i'},
                'mode': 'SINGLE',
                'phase': 'LIQUID',
                'eos_model': 'PR'
            }
        }


class ComponentLiquidFugacityResult(BaseModel):
    '''
    '''
    phase: list[str] = Field(
        ...,
        description="List of phases considered in the calculation, e.g., ['liquid']"
    )
    component: list[str] = Field(
        ...,
        description="Component name for which fugacity is calculated"
    )
    results: Dict[str, ComponentLiquidFugacityPhaseResult] = Field(
        ...,
        description="Dictionary mapping phase names to their respective LiquidFugacityPhaseResult models"
    )


# SECTION: mixture fugacity results model for multiple components
# NOTE: The mixture might be in different phases as vapor and liquid
# {
#     'vapor': {
#         'CO2': {
#             'mole_fraction': 0.15,
#             'pressure': {'value': 1000000.0, 'unit': 'Pa', 'symbol': 'P'},
#             'temperature': {'value': 444.0, 'unit': 'K', 'symbol': 'T'},
#             'molar_volume': {'value': 0.00343938395902379, 'unit': 'm3/mol', 'symbol': 'MoVo'},
#             'compressibility_coefficient': {'value': 0.931671941173366, 'unit': 'dimensionless', 'symbol': 'Z'},
#             'fugacity_coefficient': {'value': 1.0637414054916856, 'unit': 'dimensionless', 'symbol': 'phi'},
#             'fugacity': {'value': 159561.21082375283, 'unit': 'Pa', 'symbol': 'Fug_MIX'},
#             'mode': 'MIXTURE',
#             'phase': 'VAPOR',
#             'eos_model': 'RK'
#         },
#         'n-butane': {
#             'mole_fraction': 0.85,
#             'pressure': {'value': 1000000.0, 'unit': 'Pa', 'symbol': 'P'},
#             'temperature': {'value': 444.0, 'unit': 'K', 'symbol': 'T'},
#             'molar_volume': {'value': 0.00343938395902379, 'unit': 'm3/mol', 'symbol': 'MoVo'},
#             'compressibility_coefficient': {'value': 0.931671941173366, 'unit': 'dimensionless', 'symbol': 'Z'},
#             'fugacity_coefficient': {'value': 1.0122167294730493, 'unit': 'dimensionless', 'symbol': 'phi'},
#             'fugacity': {'value': 860384.2200520919, 'unit': 'Pa', 'symbol': 'Fug_MIX'},
#             'mode': 'MIXTURE',
#             'phase': 'VAPOR',
#             'eos_model': 'RK'
#         }
#     },
#     'phase': ['vapor'],
#     'component': ['CO2', 'n-butane']
# }

class MixtureFugacityResult(BaseModel):
    '''
    '''
    phase: list[str] = Field(
        ...,
        description="List of phases considered in the calculation, e.g., ['vapor', 'liquid']"
    )
    components: list[str] = Field(
        ...,
        description="List of component names for which fugacity is calculated"
    )
    results: Dict[str, Dict[str, ComponentGasFugacityPhaseResult]] = Field(
        ...,
        description="Dictionary mapping phase names to another dictionary that maps component names to their respective GasFugacityPhaseResult models"
    )
