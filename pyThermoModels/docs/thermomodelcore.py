# import package/modules
from typing import Dict, List, Union, Literal, Optional, Tuple
import pycuc
import json
# local
from .eoscore import EOSCore
from .fugacity import FugacityClass
from .thermodb import ThermoDB
from .thermolinkdb import ThermoLinkDB
from ..plugin import ReferenceManager
from .fugacitycore import FugacityCore
from ..utils import eos_model_name
from .eosutils import EOSUtils


class ThermoModelCore(ThermoDB, ThermoLinkDB, ReferenceManager):
    """
    ThermoModelCore class for thermodynamic model calculations.
    """

    _input = {}
    _references = {}
    # results
    _results = {}

    def __init__(self):
        # init ThermoDB
        ThermoDB.__init__(self)
        ThermoLinkDB.__init__(self)
        ReferenceManager.__init__(self)

        # reference plugin
        self._references = self.load_reference()

    @property
    def input(self):
        return self._input

    def check_thermodb(self):
        '''
        Check thermo databook (thermodb file created by pyThermoDB)

        Parameters
        ----------
        None

        Returns
        -------
        thermodb : dict
            dict of components in the thermodb
        '''
        try:
            # check name
            return self.thermodb
        except Exception as e:
            raise Exception('Checking thermodb failed! ', e)

    def check_fugacity_reference(self, eos_model: str,
                                 res_format: Literal['dict',
                                                     'json', 'string'] = 'dict'
                                 ) -> Union[Dict, str]:
        '''
        Check fugacity reference

        Parameters
        ----------
        eos_model : str
            equation of state name (SRK, PR, etc.)
        res_format : str
            result format, `dict`: dictionary, `json`: json format, `string`: string format

        Returns
        -------
        res : dict or str
            fugacity reference
        '''
        try:
            # check empty eos_model
            if not eos_model:
                raise Exception('Empty equation of state name!')

            # eos model
            eos_model = eos_model.upper()
            # reference
            reference = self._references.get(eos_model, None)

            # check
            if not reference:
                raise Exception('Invalid equation of state name!')

            # NOTE:
            if res_format == 'dict':
                return reference
            elif res_format == 'json':
                return json.dumps(reference, indent=4)
            elif res_format == 'string':
                return str(reference)
            else:
                raise Exception('Invalid result format!')
        except Exception as e:
            raise Exception("Calculating the Fugacity failed!, ", e)

    def cal_fugacity(self,
                     model_name: Literal['SRK', 'PR', 'VdW'],
                     model_input: Dict,
                     model_source: Dict,
                     solver_method: Literal['ls',
                                            'newton', 'fsolve'] = 'ls',
                     liquid_fugacity_mode: Literal['EOS'] = 'EOS'):
        '''
        Starts calculating fugacity for the single component

        Parameters
        ----------
        model_name: str
            eos model name, `SRK`: Soave-Redlich-Kwong, `PR`: Peng-Robinson, `VdW`: Van der Waals
        model_input: dict
            model input
                - phase (Optional): str, `VAPOR`: vapor phase, `LIQUID`: liquid phase, `VAPOR-LIQUID`: vapor-liquid phase, `SUPERCRITICAL`: supercritical phase
                - pressure: list, pressure in SI unit, such as `[1.2*1e5, 'Pa']`
                - temperature: list, temperature in SI unit, such as `[300, 'K']`
        solver_method: str
            solver method, `ls`: least square method, `newton`: newton method, `fsolve`: fsolve method
        model_source: dict
            datasource and equationsource needed for fugacity calculation
        liquid_fugacity_mode: str
            liquid fugacity method, `Poynting`: Poynting method (soon), `EOS`: Equation of state (lowest Z)

        Returns
        -------
        res: tuple
            compressibility factor (Z), fugacity coefficient (phi), eos parameters, and fugacity parameters

        Notes
        -----
        ### solver_method:
            - `ls`: least square method
            - `newton`: newton method
            - `fsolve`: fsolve method

        ### root_analysis_set:

        root analysis set, `None`: default (calculation performed according to phase provided), `1`: 3 roots (VAPOR-LIQUID), `2`: 1 root (LIQUID), `3`: 1 root (VAPOR), `4`: 1 root (SUPERCRITICAL). Otherwise, eos root analysis is performed to determine the phase of the component at the provided operating conditions including temperature and pressure.

            - set [1]: 3 roots (VAPOR-LIQUID)
                - At `T < Tc` and `P = Psat`, 3 real roots → smallest is liquid, largest is vapor.
            - set [2]: 1 root (LIQUID)
                - At `T < Tc` and `P > Psat`, EOS may give 1 or 3 roots → use smallest (liquid).
            - set [3]: 1 root (VAPOR)
                - At `T < Tc` and `P < Psat`, EOS may give 1 or 3 roots → use largest (vapor).
            - set [4]: 1 root (SUPERCRITICAL)
                - At `T > Tc`, only 1 real root → fluid is supercritical (vapor-like or liquid-like).

        Examples
        --------
        `model_input` defines as:

        ```python
        # single component
        component = 'CO2

        # model input
        # eos model
        eos_model = 'SRK'

        # component phase (Optional)
        phase = "VAPOR"

        # temperature [K] (other units are automatically converted to K)
        T = 300

        # pressure [Pa] (other units are automatically converted to Pa)
        P = 1.2*1e5

        # model input
        model_input = {
            "component": component,
            "phase": phase, # optional
            "pressure": [P, 'Pa'],
            "temperature": [T, 'K'],
        }

        # model source
        model_source = {
            "datasource": datasource,
            "equationsource": equationsource
        }

        # fugacity calculation
        res = tm.cal_fugacity(
            model_name=eos_model, model_input=model_input, model_source=model_source)
        ```
        '''
        try:
            # SECTION: set input parameters
            # eos
            eos_model = model_name.upper()
            eos_model = eos_model_name(eos_model)

            # SECTION: phase
            phase = model_input.get('phase', None)

            # check
            if phase is not None:
                phase = phase.upper()
                # check phase
                if phase not in ['VAPOR', 'LIQUID', 'VAPOR-LIQUID', 'SUPERCRITICAL']:
                    raise Exception('Invalid phase provided!')

            # SECTION: calculation mode
            calculation_mode = 'single'

            # SECTION: component list
            component = model_input.get('component', None)
            # check
            if component is None or component == 'None':
                raise Exception('Component name is not provided!')

            # set
            components = [component.strip()]

            # SECTION: mole fraction
            mole_fraction = [1]

            # SECTION: operating conditions
            # NOTE: check temperature and pressure
            if 'pressure' not in model_input.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in model_input.keys():
                raise Exception('No temperature in operating conditions!')

            # operating conditions
            operating_conditions = {
                "pressure": model_input["pressure"],
                "temperature": model_input["temperature"]
            }

            # SECTION: check eos roots
            if phase is None:
                # check eos roots for single component
                eos_roots_analysis = self.check_eos_roots_single_component(
                    model_name=eos_model, model_input=model_input,
                    model_source=model_source)

                # NOTE: eos parms
                # set phase
                phase = str(eos_roots_analysis['phase'])

            # check phase
            if phase is None:
                raise Exception('Phase is not provided!')

            # user set phase for the calculation
            eos_parms = {
                'phase': phase,
                'eos-model': eos_model,
                'mode': calculation_mode,
                'liquid-fugacity-mode': liquid_fugacity_mode
            }

            # SECTION: set datasource and equationsource
            # NOTE: check if datasource and equationsource are provided in model_input
            # datasource
            datasource = model_source.get('datasource', {})
            # equationsource
            equationsource = model_source.get('equationsource', {})
            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # SECTION: reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.set_datasource(components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                components, reference)

            # SECTION: init fugacity core
            FugacityCoreC = FugacityCore(
                component_datasource,
                equation_equationsource,
                components,
                operating_conditions,
                eos_parms
            )

            # SECTION: calculation mode
            res = FugacityCoreC.fugacity_cal(
                mole_fraction,
                solver_method=solver_method)

            return res
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)

    def check_eos_roots_single_component(self,
                                         model_name: Literal['SRK', 'PR', 'VdW'],
                                         model_input: Dict,
                                         model_source: Dict,
                                         **kwargs) -> Dict:
        '''
        Check eos roots for the single component at different temperature and pressure.

        Parameters
        ----------
        model_name: str
            eos model name, `SRK`: Soave-Redlich-Kwong, `PR`: Peng-Robinson,
        model_input: Dict
            model input as:
                - component: str, component name
                - phase: Optional[str], `VAPOR`: vapor phase, `LIQUID`: liquid phase, `VAPOR-LIQUID`: vapor-liquid phase, `SUPERCRITICAL`: supercritical phase
                - pressure: list, pressure and unit, such as `[1.2*1e5, 'Pa']`
                - temperature: list, temperature and unit, such as `[300, 'K']`
        model_source: Dict
            datasource and equationsource needed for fugacity calculation as:
                - datasource: dict, datasource for the component (`generated by PyThermoDB`)
                - equationsource: dict, equationsource for the component (`generated by PyThermoDB`)
        **kwargs: Optional[Dict]
            additional arguments
                - tolerance: float, tolerance for the calculation (default: 1e-4)

        Returns
        -------
        res: dict
            eos root analysis

        Notes
        -----
        The phase of the component is determined based on the temperature and pressure provided in the model_input.

        1. At T < Tc and P = Psat, 3 real roots → smallest is liquid, largest is vapor (vapor-liquid).
        2. At T < Tc and P > Psat, EOS may give 1 or 3 roots → use smallest (liquid).
        3. At T < Tc and P < Psat, EOS may give 1 or 3 roots → use largest (vapor).
        4. At T = Tc, one real root (critical point) → fluid is at critical state.
        5. At T > Tc, only 1 real root → fluid is supercritical (vapor-like or liquid-like).

        Pressure and temperature are in SI units (Pa, K), if not provided in the model_input, they are automatically converted.

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics
        '''
        try:
            # SECTION: keywords
            # tolerance
            tolerance = kwargs.get('tolerance', 1e-4)

            # SECTION: set input parameters
            # NOTE: eos model
            eos_model = model_name.upper()
            eos_model = eos_model_name(eos_model)

            # NOTE: component list
            component = model_input.get('component', None)
            if component is None or component == 'None':
                raise Exception('Component name is not provided!')

            # set
            components = [component]

            # NOTE: check temperature and pressure
            if 'pressure' not in model_input.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in model_input.keys():
                raise Exception('No temperature in operating conditions!')

            # NOTE: operating conditions
            operating_conditions = {
                "pressure": model_input["pressure"],
                "temperature": model_input["temperature"]
            }

            # SECTION: set datasource and equationsource
            # NOTE: check if datasource and equationsource are provided in model_input
            # datasource
            datasource = model_source.get('datasource', {})
            # equationsource
            equationsource = model_source.get('equationsource', {})
            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # SECTION: reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.set_datasource(components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                components, reference)

            # SECTION: operating conditions
            # pressure [Pa]
            P = pycuc.to(
                operating_conditions["pressure"][0],
                f"{operating_conditions['pressure'][1]} => Pa")
            # temperature [K]
            T = pycuc.to(
                operating_conditions["temperature"][0],
                f"{operating_conditions['temperature'][1]} => K")

            # SECTION: init
            EOSUtilsC = EOSUtils(component_datasource, equation_equationsource)

            # SECTION: eos root analysis
            res = EOSUtilsC.eos_root_analysis(
                P, T, components, tolerance=tolerance)

            return res[component]
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)

    def cal_fugacity_mixture(self,
                             model_name: Literal['SRK', 'PR'],
                             model_input: Dict,
                             model_source: Dict,
                             solver_method: Literal['ls',
                                                    'newton', 'fsolve'] = 'ls',
                             root_analysis_set: Optional[int] = None,
                             liquid_fugacity_mode: Literal['EOS'] = 'EOS'):
        '''
        Starts calculating fugacity for the single and multi-component systems

        Parameters
        ----------
        model_name: str
            eos model name, `SRK`: Soave-Redlich-Kwong, `PR`: Peng-Robinson,
        model_input: dict
            model input
                - phase: str, `VAPOR`: vapor phase, `LIQUID`: liquid phase, `VAPOR-LIQUID`: vapor-liquid phase, `SUPERCRITICAL`: supercritical phase
                - feed-specification: dict, such as `{'CO2': 1.0}` or `{'CO2': 0.5, 'N2': 0.5}`
                - pressure: list, pressure in SI unit, such as `[1.2*1e5, 'Pa']`
                - temperature: list, temperature in SI unit, such as `[300, 'K']`
        solver_method: str
            solver method, `ls`: least square method, `newton`: newton method, `fsolve`: fsolve method
        model_source: dict
            datasource and equationsource needed for fugacity calculation
        root_analysis_set: Optional[int]
            root analysis set, `None`: default (calculation performed according to phase provided), `1`: 3 roots (VAPOR-LIQUID), `2`: 1 root (LIQUID), `3`: 1 root (VAPOR), `4`: 1 root (SUPERCRITICAL)
        liquid_fugacity_mode: str
            liquid fugacity method, `Poynting`: Poynting method (soon), `EOS`: Equation of state (lowest Z)

        Returns
        -------
        res: tuple
            compressibility factor (Z), fugacity coefficient (phi), eos parameters, and fugacity parameters

        Notes
        -----
        ### solver_method:
            - `ls`: least square method
            - `newton`: newton method
            - `fsolve`: fsolve method

        ### root_analysis_set:
            - set [1]: 3 roots (VAPOR-LIQUID)
                - At `T < Tc` and `P = Psat`, 3 real roots → smallest is liquid, largest is vapor.
            - set [2]: 1 root (LIQUID)
                - At `T < Tc` and `P > Psat`, EOS may give 1 or 3 roots → use smallest (liquid).
            - set [3]: 1 root (VAPOR)
                - At `T < Tc` and `P < Psat`, EOS may give 1 or 3 roots → use largest (vapor).
            - set [4]: 1 root (SUPERCRITICAL)
                - At `T > Tc`, only 1 real root → fluid is supercritical (vapor-like or liquid-like).

        Examples
        --------
        `model_input` defines as:

        ```python
        # single component
        N0s = {
            "CO2": 1.0
        }

        # multi-component
        N0s = {
            "CO2": 0.5,
            "N2": 0.5
        }

        # Others
        # model input
        # eos model
        eos_model = 'SRK'
        # component phase
        phase = "VAPOR"
        # temperature [K]
        T = 300
        # pressure [Pa]
        P = 1.2*1e5

        # model input
        model_input = {
            "phase": phase,
            "feed-specification": N0s,
            "pressure": [P, 'Pa'],
            "temperature": [T, 'K'],
        }

        # model source
        model_source = {
            "datasource": datasource,
            "equationsource": equationsource
        }
        ```
        '''
        try:
            # SECTION: set input parameters
            # eos
            eos_model = model_name.upper()
            eos_model = eos_model_name(eos_model)

            # SECTION: phase
            phase = model_input.get('phase', 'VAPOR')
            phase = phase.upper()
            # check phase
            if phase not in ['VAPOR', 'LIQUID', 'VAPOR-LIQUID', 'SUPERCRITICAL']:
                raise Exception('Invalid phase provided!')

            # NOTE: calculation mode
            calculation_mode = 'single'
            # component number
            component_num = 0

            # input
            feed_spec = model_input.get('feed-specification', None)
            # check
            if feed_spec is None or feed_spec == 'None':
                raise Exception('Feed specification is not provided!')

            # component list
            components = []
            # mole fraction
            mole_fraction = []
            # looping through feed-specification
            for key, value in feed_spec.items():
                components.append(key)
                mole_fraction.append(value)

            # check
            if isinstance(components, list):
                # set
                component_num = len(components)
                # single
                if len(components) == 1:
                    calculation_mode = 'single'
                # mixture
                else:
                    calculation_mode = 'mixture'
            else:
                raise Exception('Components list not provided!')

            # check if multi
            if calculation_mode == 'mixture':
                # check
                if len(mole_fraction) != component_num:
                    raise Exception('Mole fraction list not provided!')

            # NOTE: check temperature and pressure
            if 'pressure' not in model_input.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in model_input.keys():
                raise Exception('No temperature in operating conditions!')

            # operating conditions
            operating_conditions = {
                "pressure": model_input["pressure"],
                "temperature": model_input["temperature"]
            }

            # SECTION: check eos roots
            eos_roots_analysis = self.check_eos_roots_single_component(
                model_name=eos_model, model_input=model_input, model_source=model_source)

            # NOTE: eos parms
            if len(eos_roots_analysis) == 1:
                eos_parms = {
                    'phase': eos_roots_analysis['phase'],
                    'eos-model': eos_model,
                    'mode': calculation_mode,
                    'liquid-fugacity-mode': liquid_fugacity_mode
                }
            else:
                eos_parms = {
                    'phase': phase,
                    'eos-model': eos_model,
                    'mode': calculation_mode,
                    'liquid-fugacity-mode': liquid_fugacity_mode
                }

            # SECTION: set datasource and equationsource
            # NOTE: check if datasource and equationsource are provided in model_input
            # datasource
            datasource = model_source.get('datasource', {})
            # equationsource
            equationsource = model_source.get('equationsource', {})
            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # SECTION: reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.set_datasource(components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                components, reference)

            # SECTION: init fugacity core
            FugacityCoreC = FugacityCore(
                component_datasource, equation_equationsource, components, operating_conditions, eos_parms)

            # NOTE: root analysis mode
            root_analysis_set = FugacityCoreC.root_analysis_mode(phase)

            # SECTION: calculation mode
            res = FugacityCoreC.fugacity_cal(
                mole_fraction, solver_method=solver_method, root_analysis_set=root_analysis_set)

            return res
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)

    def check_eos_roots_multi_component(self, model_name: Literal['SRK', 'PR'], model_input: Dict, model_source: Dict) -> List[Dict]:
        '''
        Check eos roots for the multi-components at different temperature and pressure.

        Parameters
        ----------
        model_name: str
            eos model name, `SRK`: Soave-Redlich-Kwong, `PR`: Peng-Robinson,
        model_input: Dict
            model input
        model_source: Dict
            datasource and equationsource needed for fugacity calculation

        Returns
        -------
        res: list
            eos root analysis

        Notes
        -----
        1. At T < Tc and P = Psat, 3 real roots → smallest is liquid, largest is vapor.
        2. At T < Tc and P > Psat, EOS may give 1 or 3 roots → use smallest (liquid).
        3. At T < Tc and P < Psat, EOS may give 1 or 3 roots → use largest (vapor).
        4. At T = Tc, one real root (critical point) → fluid is at critical state.
        5. At T > Tc, only 1 real root → fluid is supercritical (vapor-like or liquid-like).

        Pressure and temperature are in SI units (Pa, K), if not provided in the model_input, they are automatically converted.

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics
        '''
        try:
            # SECTION: set input parameters
            # NOTE: eos model
            eos_model = model_name.upper()
            eos_model = eos_model_name(eos_model)

            # NOTE: calculation mode
            calculation_mode = 'single'
            # component number
            component_num = 0

            # NOTE: feed specification
            feed_spec = model_input.get('feed-specification')
            # check
            if feed_spec is None or feed_spec == 'None':
                raise Exception('Feed specification is not provided!')

            # NOTE: component list
            components = []
            # mole fraction
            mole_fraction = []
            # looping through feed-specification
            for key, value in feed_spec.items():
                components.append(key)
                mole_fraction.append(value)

            # check
            if isinstance(components, list):
                # set
                component_num = len(components)
                # single
                if len(components) == 1:
                    calculation_mode = 'single'
                # mixture
                else:
                    calculation_mode = 'mixture'
            else:
                raise Exception('Components list not provided!')

            # NOTE: check if multi
            if calculation_mode == 'mixture':
                # check
                if len(mole_fraction) != component_num:
                    raise Exception('Mole fraction list not provided!')

            # NOTE: check temperature and pressure
            if 'pressure' not in model_input.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in model_input.keys():
                raise Exception('No temperature in operating conditions!')

            # NOTE: operating conditions
            operating_conditions = {
                "pressure": model_input["pressure"],
                "temperature": model_input["temperature"]
            }

            # SECTION: set datasource and equationsource
            # NOTE: check if datasource and equationsource are provided in model_input
            # datasource
            datasource = model_source.get('datasource', {})
            # equationsource
            equationsource = model_source.get('equationsource', {})
            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # SECTION: reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.set_datasource(components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                components, reference)

            # SECTION: operating conditions
            # pressure [Pa]
            P = pycuc.to(
                operating_conditions["pressure"][0], f"{operating_conditions['pressure'][1]} => Pa")
            # temperature [K]
            T = pycuc.to(
                operating_conditions["temperature"][0], f"{operating_conditions['temperature'][1]} => K")

            # SECTION: init
            EOSUtilsC = EOSUtils(component_datasource, equation_equationsource)

            # SECTION: eos root analysis
            res = EOSUtilsC.eos_root_analysis(P, T, components)

            return res
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)

    def cal_activity(self, model_name: Literal['NRTL', 'UNIQUAC'], model_input: Dict):
        '''
        Initializes activity coefficient calculation

        Parameters
        ----------
        model_name: str
            activity coefficient model name, `NRTL`: Non-Random Two-Liquid, `UNIQUAC`: Universal Quasi-Chemical
        model_input: dict
            model input

        Returns
        -------
        res : dict
            activity coefficient
        '''
        try:
            pass
        except Exception as e:
            raise Exception("Activity coefficient calculation failed!", e)
