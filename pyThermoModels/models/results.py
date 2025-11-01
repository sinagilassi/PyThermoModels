# import libs
from pydantic import BaseModel, Field
from typing import Dict, List, Union, Literal, Optional

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
    symbol: Optional[str] = None


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
    results: Dict[str, ComponentLiquidFugacityPhaseResult | ComponentGasFugacityPhaseResult] = Field(
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

# SECTION: component eos root results model
# {
#     'component_name': 'propane-g',
#     'pressure': 1000000.0,
#     'pressure_unit': 'Pa',
#     'temperature': 300.1,
#     'temperature_unit': 'K',
#     'root': 2,
#     'root-no': '1 real root (liquid)',
#     'phase': 'LIQUID',
#     'vapor_pressure': 999723.1044,
#     'vapor_pressure_unit': 'Pa',
#     'critical_temperature': 369.83,
#     'critical_temperature_unit': 'K',
#     'critical_pressure': 4248000.0,
#     'critical_pressure_unit': 'Pa',
#     'tolerance': 0.1,
#     'vapor_pressure_check': -276.8956000000471,
#     'temperature_equality_value': 69.72999999999996,
#     'pressure_equality_check': False,
#     'temperature_equality_check': False,
#     'message': 'Component propane-g at T=300.1 K and P=1000000.0 Pa is in liquid phase.'
# }


class ComponentEosRootResult(BaseModel):
    '''
    Component EOS root analysis result model
    '''
    component_name: str = Field(
        ...,
        description="Name of the component"
    )
    pressure: PropertyValue = Field(
        ...,
        description="Pressure value with unit"
    )
    temperature: PropertyValue = Field(
        ...,
        description="Temperature value with unit"
    )
    root: int = Field(
        ...,
        description="Selected root number"
    )
    root_no: str = Field(
        ...,
        description="Description of the number of real roots"
    )
    phase: str = Field(
        ...,
        description="Predicted phase of the component"
    )
    vapor_pressure: PropertyValue = Field(
        ...,
        description="Vapor pressure value with unit"
    )
    critical_temperature: PropertyValue = Field(
        ...,
        description="Critical temperature value"
    )
    critical_pressure: PropertyValue = Field(
        ...,
        description="Critical pressure value"
    )
    tolerance: float = Field(
        ...,
        description="Tolerance used in the analysis"
    )
    vapor_pressure_check: float = Field(
        ...,
        description="Difference between system pressure and vapor pressure"
    )
    temperature_equality_value: float = Field(
        ...,
        description="Difference between system temperature and saturation temperature"
    )
    pressure_equality_check: bool = Field(
        ...,
        description="Check if system pressure equals vapor pressure within tolerance"
    )
    temperature_equality_check: bool = Field(
        ...,
        description="Check if system temperature equals saturation temperature within tolerance"
    )
    message: str = Field(
        ...,
        description="Summary message of the EOS root analysis"
    )


# NOTE: mixture eos root results model
# {
#     'component_name': 'carbon dioxide-g | n-butane-g',
#     'pressure': 1000000.0,
#     'pressure_unit': 'Pa',
#     'temperature': 444.0,
#     'temperature_unit': 'K',
#     'bubble_pressure': 18796812.05906,
#     'bubble_pressure_unit': 'Pa',
#     'dew_point_pressure': 5845352.760145266,
#     'dew_point_pressure_unit': 'Pa',
#     'bubble_point_temperature': 444.0,
#     'bubble_point_temperature_unit': 'K',
#     'dew_point_temperature': 444.0,
#     'dew_point_temperature_unit': 'K',
#     'phase': 'VAPOR',
#     'tolerance': 0.1,
#     'message': 'Mixture carbon dioxide-g | n-butane-g at T=444.0 K and P=1000000.0 Pa is in vapor phase as P < Dew Point Pressure (5845352.760145266 Pa).'
# }

class MixtureEosRootResult(BaseModel):
    '''
    Mixture EOS root analysis result model
    '''
    mixture_name: str = Field(
        ...,
        description="Names of the components in the mixture"
    )
    pressure: PropertyValue = Field(
        ...,
        description="Pressure value with unit"
    )
    temperature: PropertyValue = Field(
        ...,
        description="Temperature value with unit"
    )
    bubble_pressure: PropertyValue = Field(
        ...,
        description="Bubble point pressure value"
    )
    dew_point_pressure: PropertyValue = Field(
        ...,
        description="Dew point pressure value"
    )
    bubble_point_temperature: PropertyValue = Field(
        ...,
        description="Bubble point temperature value"
    )
    dew_point_temperature: PropertyValue = Field(
        ...,
        description="Dew point temperature value"
    )
    phase: str = Field(
        ...,
        description="Predicted phase of the mixture"
    )
    tolerance: float = Field(
        ...,
        description="Tolerance used in the analysis"
    )
    message: str = Field(
        ...,
        description="Summary message of the EOS root analysis for the mixture"
    )
