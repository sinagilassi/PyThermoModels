# THERMODB

# import packages/modules

# local


class ThermoLinkDB:
    # vars
    _datasource = {}
    _equationsource = {}
    #
    _thermodb_component = {}

    def __init__(self):
        # load reference
        pass

    @property
    def datasource(self):
        return self._datasource

    @property
    def equationsource(self):
        return self._equationsource

    @property
    def thermodb_component(self):
        return self._thermodb_component

    def set_thermodb_link(self, datasource, equationsource):
        '''
        Sets thermodb link

        Parameters
        ----------
        datasource : dict
            dict data of the datasource
        equationsource : dict
            dict data of the equationsource

        Returns
        -------
        bool
            True if successful, False otherwise
        '''
        try:
            # set datasource | equationsource
            self._datasource = {}
            self._datasource = datasource
            self._equationsource = {}
            self._equationsource = equationsource

            # datasource component
            # check
            if self._datasource is not None:
                ds_components = datasource.keys()
            else:
                ds_components = []

            # equationsource component
            # check
            if self._equationsource is not None:
                eq_components = equationsource.keys()
            else:
                eq_components = []

            # thermodb component
            self._thermodb_component = list(
                set(ds_components) & set(eq_components))
            # res
            return True
        except Exception as e:
            raise Exception('Loading record failed!, ', e)

    def set_datasource(self, components, reference):
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
        datasource : dict
            dict data of the datasource
        '''
        try:
            # check reference
            if reference is None or reference == 'None':
                raise Exception('Empty reference!')

            # reference
            # get all dependent data
            dependent_data_src = reference.get('DEPENDANT-DATA')
            dependent_data = []
            # check
            if dependent_data_src is not None and dependent_data_src != 'None':
                for item, value in dependent_data_src.items():
                    _item_symbol = str(value['symbol'])
                    # add symbol
                    dependent_data.append(_item_symbol)

            # datasource
            datasource = {}

            # NOTE: checking components
            for component in components:
                if component in self._thermodb_component:
                    # set
                    datasource[component] = {}
                    # parms
                    for item in dependent_data:
                        # get value
                        _val = self.datasource[component].get(item, None)

                        # check
                        if _val is not None and _val != 'None':
                            # save
                            datasource[component][item] = _val

            # NOTE: checking general data defined as 'NRTL' and 'UNIQUAC'
            for records in self._thermodb_component:
                # set upper
                records_ = records.upper()

                # check item
                if records_ in ['UNIQUAC', 'NRTL']:
                    # set
                    datasource[records_] = {}
                    # parms
                    for item in dependent_data:
                        # get value
                        _val = self.datasource[records].get(item, None)

                        # check
                        if _val is not None and _val != 'None':
                            # save
                            datasource[records_][item] = _val

            # res
            return datasource
        except Exception as e:
            raise Exception('Building datasource failed!, ', e)

    def set_equationsource(self, components, reference):
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
        equationsource : dict
            dict data of the equationsource
        '''
        try:
            # check reference
            if reference is None or reference == 'None':
                raise Exception('Empty reference!')

            # reference
            # get all dependent data
            dependent_data_src = reference.get('DEPENDANT-EQUATIONS')
            dependent_data = []
            # check
            if dependent_data_src is not None and dependent_data_src != 'None':
                for item, value in dependent_data_src.items():
                    _item_symbol = value['symbol']
                    dependent_data.append(_item_symbol)

            # equationsource
            equationsource = {}
            for component in components:
                if component in self._thermodb_component:
                    # set
                    equationsource[component] = {}
                    # parms
                    for item in dependent_data:
                        # get value
                        _val = self.equationsource[component][item]
                        # save
                        equationsource[component][item] = _val
            # res
            return equationsource
        except Exception as e:
            raise Exception('Building equationsource failed!, ', e)
