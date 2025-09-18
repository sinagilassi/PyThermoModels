# import libs
from pydantic import BaseModel, Field
from typing import Any, Dict, List, Optional, Tuple, Type, Union, Literal
from pythermodb_settings.models import Temperature, Pressure


class OperatingConditions(BaseModel):
    '''
    Operating conditions model for thermodynamic calculations

    Attributes
    ----------
    temperature: list
        Temperature and a unit
    pressure: list
        Pressure and a unit
    phase: str, optional
        Phase type, options are:
            - `VAPOR`: vapor phase
            - `LIQUID`: liquid phase
            - `VAPOR-LIQUID`: vapor-liquid phase
            - `SUPERCRITICAL`: supercritical phase
    '''
    temperature: Temperature = Field(
        ...,
        description="Temperature and a unit`"
    )
    pressure: Pressure = Field(
        ...,
        description="Pressure and a unit`"
    )
    phase: Optional[Literal[
        'VAPOR', 'LIQUID', 'VAPOR-LIQUID', 'SUPERCRITICAL'
    ]] = Field(
        None,
        description="Phase type, options are: `VAPOR`: vapor phase, `LIQUID`: liquid phase, `VAPOR-LIQUID`: vapor-liquid phase, `SUPERCRITICAL`: supercritical phase"
    )
