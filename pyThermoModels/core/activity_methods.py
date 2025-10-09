# import libs
import logging
from typing import Any, Dict, Literal, List
from pythermodb_settings.models import Component, Temperature, Pressure
from pythermodb_settings.utils import set_component_id
from pyThermoLinkDB.models import ModelSource
# local
from ..docs import ThermoModelCore
from ..utils import set_feed_specification

# NOTE: logger
logger = logging.getLogger(__name__)
