# import package/modules
import pandas as pd
import pycuc
# local
from .eoscore import EOSCore
from .fugacity import FugacityClass
from .thermodb import ThermoDB
from .thermolinkdb import ThermoLinkDB
from ..plugin import ReferenceManager
from .fugacitycore import FugacityCore
from ..utils import eos_model_name
from .eosutils import EOSUtils


class Manager(ThermoDB, ThermoLinkDB, ReferenceManager):

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

    def fugacity_check_reference(self, eos_model, dataframe=False):
        '''
        Check fugacity reference

        Parameters
        ----------
        eos_model : str
            equation of state name
        dataframe : bool, optional
            return dataframe, default False

        Returns
        -------
        ref : dict
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

            # check
            if dataframe:
                return pd.DataFrame(reference)
            else:
                # return
                return reference
        except Exception as e:
            raise Exception("Calculating the Fugacity failed!, ", e)

    def fugacity_cal_init(self, model_input):
        '''
        Initialize fugacity calculation

        Parameters
        ----------
        model_input: dict
            model input

        Returns
        -------
        res : list
            fugacity and fugacity coefficient
        '''
        try:
            # eos
            eos_model = model_input.get('eos-model', 'PR')
            eos_model = eos_model.upper()
            # phase
            phase = model_input.get('phase', 'gas')
            phase = phase.upper()
            # component
            components = model_input["components"]
            # mole fraction
            mole_fraction = model_input["mole-fraction"]
            # operating conditions
            operating_conditions = model_input["operating_conditions"]

            # reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.build_datasource(components, reference)

            # # * init eos class
            EOSCoreC = EOSCore(
                component_datasource, components, eos_model, mole_fraction, operating_conditions)

            # ! calculate compressibility factor Z
            # select method
            select_eos_model = {
                "PR": lambda: EOSCoreC._eosPR()
            }

            # check
            if eos_model not in select_eos_model.keys():
                raise Exception("EOS model not found!")

            # res
            # Z
            res1 = select_eos_model.get(eos_model)()

            # ! init fugacity class
            fugacityC = FugacityClass(
                component_datasource, components, res1, operating_conditions)

            # select method
            selectFugacity = {
                "PR": lambda x: fugacityC.FugacityPR(x)
            }

            # res
            res2 = selectFugacity.get(eos_model)(phase)

            return res1, res2

        except Exception as e:
            raise Exception("Initializing fugacity calculation failed!, ", e)

    def cal_fugacity_coefficient(self, model_input: dict, solver_method='ls', root_analysis_set=None, liquid_fugacity_calculation_method='Poynting'):
        '''
        Calculate fugacity coefficient for the single and multi-component systems

        Parameters
        ----------
        model_input: dict
            model input
        solver_method: str
            solver method
        root_analysis_set: int
            root analysis set
        liquid_fugacity_calculation_method: str
            liquid fugacity method, `Poynting`: Poynting method, `EOS`: Equation of state (lowest Z)

        Returns
        -------
        res: tuple
            compressibility factor (Z), fugacity coefficient (phi), and eos parms

        Notes
        -----
        ### solver_method:
            - ls: least square method
            - newton: newton method
            - fsolve: fsolve method

        ### root_analysis_set:
            - [1]: 3 roots (VAPOR-LIQUID)
            - [2]: 1 root (LIQUID)
            - [3]: 1 root (VAPOR)
            - [4]: 1 root (SUPERCRITICAL)

        ### `model_input` must be defined as:
            ```python
            ## single component
            # component list
            comp_list = ["CO2"]
            # mole fraction
            MoFri = []

            ## multi-component
            # component list
            comp_list = ["EtOH", "MeOH"]
            # mole fraction
            MoFri = [0.75, 0.25]

            ## Others
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
                "eos-model": eos_model,
                "phase": phase,
                "components": comp_list,
                "mole-fraction": MoFri,
                "operating_conditions": {
                    "pressure": [P, 'Pa'],
                    "temperature": [T, 'K'],
                },
                "datasource": datasource,
                "equationsource": equationsource
            }
            ```
        '''
        try:
            # eos
            eos_model = model_input.get('eos-model', 'SRK')
            eos_model = eos_model.upper()
            eos_model = eos_model_name(eos_model)

            # phase
            phase = model_input.get('phase', 'VAPOR')
            phase = phase.upper()

            # calculation mode
            calculation_mode = 'single'
            # component number
            component_num = 0

            # input
            feed_spec = model_input.get('feed-spec')
            # check
            if feed_spec is None or feed_spec == 'None':
                raise Exception('Feed specification is not provided!')

            # component list
            components = []
            # mole fraction
            mole_fraction = []
            # looping through feed-spec
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

            # operating conditions
            operating_conditions = model_input["operating-conditions"]
            # check temperature and pressure
            if 'pressure' not in operating_conditions.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in operating_conditions.keys():
                raise Exception('No temperature in operating conditions!')

            # eos parms
            eos_parms = {
                'phase': phase,
                'eos-model': eos_model,
                'mode': calculation_mode,
                'liquid-fugacity-calculation-method': liquid_fugacity_calculation_method
            }

            # datasource
            datasource = model_input.get('datasource', {})
            # equationsource
            equationsource = model_input.get('equationsource', {})
            # set thermodb link
            link_status = self.set_thermodb_link(datasource, equationsource)
            # check
            if not link_status:
                raise Exception('Thermodb link failed!')

            # reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.set_datasource(components, reference)
            # build equation source
            equation_equationsource = self.set_equationsource(
                components, reference)

            # init
            FugacityCoreC = FugacityCore(
                component_datasource, equation_equationsource, components, operating_conditions, eos_parms)

            # root analysis mode
            if root_analysis_set is None:
                root_analysis_set = FugacityCoreC.root_analysis_mode()

            # calculation
            res = FugacityCoreC.fugacity_cal(
                mole_fraction, solver_method=solver_method, root_analysis_set=root_analysis_set)

            return res
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)

    def check_eos_roots(self, model_input):
        '''
        Checking eos roots

        Parameters
        ----------
        model_input: dict
            model input
        solver_method: str
            solver method
        root_analysis_set: int
            root analysis set

        Returns
        -------
        res: list
            eos root analysis

        Notes
        -----

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics
        '''
        try:
            # eos
            eos_model = model_input.get('eos-model', 'SRK')
            eos_model = eos_model.upper()
            eos_model = eos_model_name(eos_model)

            # phase
            phase = model_input.get('phase', 'VAPOR')
            phase = phase.upper()

            # calculation mode
            calculation_mode = 'single'
            # component number
            component_num = 0

            # component
            components = model_input["components"]

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

            # mole fraction
            mole_fraction = model_input["mole-fraction"]

            # check if multi
            if calculation_mode == 'mixture':
                # check
                if len(mole_fraction) != component_num:
                    raise Exception('Mole fraction list not provided!')

            # operating conditions
            operating_conditions = model_input["operating_conditions"]
            # check temperature and pressure
            if 'pressure' not in operating_conditions.keys():
                raise Exception('No pressure in operating conditions!')

            if 'temperature' not in operating_conditions.keys():
                raise Exception('No temperature in operating conditions!')

            # reference for eos
            reference = self._references.get(eos_model, None)

            # build datasource
            component_datasource = self.build_datasource(components, reference)
            # build equation source
            equation_equationsource = self.build_equationsource(
                components, reference)

            # pressure [Pa]
            P = pycuc.to(
                operating_conditions["pressure"][0], f"{operating_conditions['pressure'][1]} => Pa")
            # temperature [K]
            T = pycuc.to(
                operating_conditions["temperature"][0], f"{operating_conditions['temperature'][1]} => K")

            # init
            EOSUtilsC = EOSUtils(component_datasource, equation_equationsource)

            # eos root analysis
            res = EOSUtilsC.eos_root_analysis(P, T, components)

            return res
        except Exception as e:
            raise Exception("Fugacity calculation failed!, ", e)
