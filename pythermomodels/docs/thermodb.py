# THERMODB

# import packages/modules
import os
import yaml
# local


class ThermoDB:
    # vars
    _thermodb = {}
    _thermodb_rule = {}

    def __init__(self):
        pass

    @property
    def thermodb(self):
        return self._thermodb

    @property
    def thermodb_rule(self):
        return self._thermodb_rule

    def config_thermodb(self, name, config_file):
        '''
        Load yml config file

        Parameters
        ----------
        name : str
            name of the record
        config_file: str
            config file path
        '''
        try:
            # check
            if config_file is None:
                # current folder
                current_folder = os.path.dirname(__file__)
                # parent folder
                parent_folder = os.path.dirname(current_folder)
                # data folder
                data_folder = os.path.join(parent_folder, 'data')
                # config file for all
                config_file = os.path.join(data_folder, 'thermodb_config.yml')

            # check file
            if not os.path.exists(config_file):
                raise Exception('Configuration file not found!')
            # load
            with open(config_file, 'r') as f:
                _ref = yaml.load(f, Loader=yaml.FullLoader)
                # set
                self._thermodb_rule[str(name).strip()] = _ref['thermodb']
        except Exception as e:
            raise Exception('Configuration failed!, ', e)

    def add_thermodb(self, name, data):
        '''
        Add new thermodb such as: CO2_thermodb

        Parameters
        ----------
        name: str
            name of the record
        data: pyThermoDB.docs.compbuilder.CompBuilder
            data of the record

        Returns
        -------
        None
        '''
        try:
            self._thermodb[name] = data
        except Exception as e:
            raise Exception('Adding new record failed!, ', e)

    def update_thermodb(self, name, data):
        '''
        Update existing record

        Parameters
        ----------
        name: str
            name of the record
        data: pyThermoDB.docs.compbuilder.CompBuilder
            data of the record

        Returns
        -------
        None
        '''
        try:
            self._thermodb[name] = data
        except Exception as e:
            raise Exception('Updating record failed!, ', e)

    def delete_thermodb(self, name):
        '''
        Delete existing record

        Parameters
        ----------
        name: str
            name of the record

        Returns
        -------
        None
        '''
        try:
            del self._thermodb[name]
        except Exception as e:
            raise Exception('Deleting record failed!, ', e)

    def build_datasource(self, components):
        '''
        Build datasource

        Parameters
        ----------
        components : list
            list of components

        Returns
        -------
        None
        '''
        try:
            for component in components:
                if component in self._thermodb:
                    # res
                    return 1
        except Exception as e:
            raise Exception('Building datasource failed!, ', e)
