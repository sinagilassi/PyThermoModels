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
        # load reference
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
            # load external local
            self.config_thermodb_local(name)

            # check
            if config_file:
                # check file
                if not os.path.exists(config_file):
                    raise Exception('Configuration file not found!')

                # set name
                name = str(name).strip()
                # load
                with open(config_file, 'r') as f:
                    _ref = yaml.load(f, Loader=yaml.FullLoader)
                    # check name exists
                    if name in _ref.keys():
                        # get record
                        record_thermodb_rule = _ref[name]
                        # looping through
                        self._thermodb_rule[name].update(record_thermodb_rule)
                    else:
                        raise Exception('Record not found!')

        except Exception as e:
            raise Exception('Configuration failed!, ', e)

    def config_thermodb_local(self, name):
        '''
        Load local yml config file

        Parameters
        ----------
        name : str
            name of the record

        '''
        try:
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

            # set name
            name = str(name).strip()
            # load
            with open(config_file, 'r') as f:
                _ref = yaml.load(f, Loader=yaml.FullLoader)

                # general data
                record_thermodb = _ref['THERMODB']['ALL']
                # set
                self._thermodb_rule[name] = {}
                self._thermodb_rule[name].update(record_thermodb)
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

    def build_datasource(self, components, reference):
        '''
        Build datasource

        Parameters
        ----------
        components : list
            list of components
        reference : dict
            dict data of the reference

        Returns
        -------
        None
        '''
        try:
            # check reference
            if reference is None or reference == 'None':
                raise Exception('Empty reference!')

            # reference
            # get all dependent data
            dependent_data_src = reference['DEPENDANT-DATA']
            dependent_data = []
            # check
            if dependent_data_src is not None:
                for item, value in dependent_data_src.items():
                    _item_symbol = value['symbol']
                    dependent_data.append(_item_symbol)

            # datasource
            datasource = {}
            for component in components:
                if component in self._thermodb:
                    # set
                    datasource[component] = {}
                    # parms
                    for item in dependent_data:
                        # check property src
                        check_property_src = 'GENERAL'
                        # check
                        if component in self.thermodb_rule.keys():
                            # property source
                            check_property_src = self.thermodb_rule[component][item][0]

                        # val
                        _val = self.thermodb[component].check_property(
                            check_property_src).get_property(str(item).strip())
                        # unit conversion

                        datasource[component][item] = _val
            # res
            return datasource
        except Exception as e:
            raise Exception('Building datasource failed!, ', e)

    def build_equationsource(self, components, reference):
        '''
        Build datasource

        Parameters
        ----------
        components : list
            list of components
        reference : dict
            dict data of the reference

        Returns
        -------
        None
        '''
        try:
            # check reference
            if reference is None or reference == 'None':
                raise Exception('Empty reference!')

            # reference
            # get all dependent data
            dependent_data_src = reference['DEPENDANT-EQUATIONS']
            dependent_data = []
            # check
            if dependent_data_src is not None and dependent_data_src != 'None':
                for item, value in dependent_data_src.items():
                    _item_symbol = value['symbol']
                    dependent_data.append(_item_symbol)

            # datasource
            datasource = {}
            for component in components:
                if component in self._thermodb:
                    # set
                    datasource[component] = {}
                    # parms
                    for item in dependent_data:
                        # check property src
                        check_property_src = None
                        # check
                        if component in self.thermodb_rule.keys():
                            # property source
                            check_property_src = self.thermodb_rule[component][item][0]

                            # val
                            _val = self.thermodb[component].check_function(
                                check_property_src)
                            # unit conversion

                            datasource[component][item] = _val
            # res
            return datasource
        except Exception as e:
            raise Exception('Building datasource failed!, ', e)
