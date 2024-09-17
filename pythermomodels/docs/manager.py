# import package/modules
import pandas as pd
# local
from .eoscore import EOSCore
from .fugacity import FugacityClass
from .thermodb import ThermoDB
from ..plugin import ReferenceManager


class Manager(ThermoDB, ReferenceManager):

    _input = {}
    _references = {}

    def __init__(self):
        # init ThermoDB
        ThermoDB.__init__(self)
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
            eos_model = model_input.get('eos-model', 'Peng_Robinson')
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
                "PENG_ROBINSON": lambda: EOSCoreC._eosPR()
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
                "PENG_ROBINSON": lambda x: fugacityC.FugacityPR(x)
            }

            # res
            res = selectFugacity.get(eos_model)(phase)

            return res

        except Exception as e:
            raise Exception("Initializing fugacity calculation failed!, ", e)
