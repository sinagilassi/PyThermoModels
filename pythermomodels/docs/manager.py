# import package/modules
import yaml
import os
import csv
# local
from .eoscore import eosCoreClass
from .fugacity import FugacityClass


class ManagerClass:
    def __init__(self, model_input):
        self.model_input = model_input

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

    def fugacity_cal_init(self):
        '''
        Initialize fugacity calculation
        '''
        try:
            # model input
            model_setting = self.model_input["model_setting"]

            eos = model_setting.get('eos', 'Peng–Robinson')
            phase = model_setting.get('phase', 'gas')

            # component
            components = self.model_input["components"]

            # operating conditions
            operating_conditions = self.model_input["operating_conditions"]

            # data source
            data_source = self.model_input["data_source"]
            # load

            # # * init eos class
            # _eosCoreClass = eosCoreClass(
            #     compData, compList, eosModel, moleFraction, params)

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
