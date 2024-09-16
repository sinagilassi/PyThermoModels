# import package/modules
import yaml
import os
import csv
# local
from .eoscore import eosCoreClass
from .fugacity import FugacityClass
from .thermodb import ThermoDB


class Manager(ThermoDB):

    _input = {}

    def __init__(self):
        # init ThermoDB
        ThermoDB.__init__(self)

    @property
    def input(self):
        return self._input

    @staticmethod
    def load_yml(yml_file):
        '''
        Load a yml file

        Parameters
        ----------
        yml_file: str
            yml file path
        '''
        try:
            # check file exists
            if os.path.exists(yml_file):
                # load yml
                with open(yml_file, 'r') as f:
                    return yaml.load(f, Loader=yaml.FullLoader)
        except Exception as e:
            raise Exception("Loading yml failed!, ", e)

    def csv_loader(self):
        '''
        list of available components in the database
        '''
        try:
            # input file path
            input_file_path = self.model_input['data_source']['path']

            # check file not exists
            if not os.path.exists(input_file_path):
                raise Exception("Loading data_source failed!")

            # component data
            comp_data = []

            csv.register_dialect('myDialect', delimiter=',',
                                 skipinitialspace=True, quoting=csv.QUOTE_MINIMAL)

            with open(input_file_path, 'r') as file:
                reader = csv.reader(file)

                # ignore header
                next(reader, None)
                next(reader, None)

                for row in reader:
                    comp_data.append([row[1], row[2]])

            if len(comp_data) > 0:
                return comp_data
            else:
                raise
        except Exception as e:
            raise Exception("Loading data filed!, ", e)

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

    def fugacity_cal_init(self, model_input):
        '''
        Initialize fugacity calculation
        '''
        try:
            # eos
            eos_model = model_input.get('eos-model', 'Peng–Robinson')
            # phase
            phase = model_input.get('phase', 'gas')
            # component
            components = model_input["components"]
            # mole fraction
            mole_fraction = model_input["mole-fraction"]
            # operating conditions
            operating_conditions = model_input["operating_conditions"]

            # thermodb
            thermodb = self.thermodb
            # thermodb rule
            thermodb_rule = self.thermodb_rule

            # check value
            a = thermodb['CO2'].check_property(
                'GENERAL').get_property('dHf_IG')['value']
            print(type(a))

            # # * init eos class
            _eosCoreClass = eosCoreClass(
                compData, components, eos_model, mole_fraction, operating_conditions)

            # # select method
            # selectEOS = {
            #     "PR": lambda: _eosCoreClass._eosPR()
            # }
            # # res
            # _eosRes = selectEOS.get(eosModel)()

            # # * init fugacity class
            # _fugacityClass = FugacityClass(compData, compList, _eosRes, phase)

            # # select method
            # selectFugacity = {
            #     "PR": lambda: _fugacityClass.FugacityPR()
            # }

            # # res
            # fugacity = selectFugacity.get(eosModel)()

            fugacity = 1

            return fugacity

        except Exception as e:
            raise Exception("Initializing fugacity calculation failed!, ", e)
