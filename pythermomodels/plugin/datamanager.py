# DATA MANAGER CLASS
# -------------------

# import packages/modules
import os
import yaml


class DataManager:

    def __init__(self):
        pass

    def load_reference(self):
        '''
        Load reference yml files having details about equations

        Parameters
        ----------
        None

        Returns
        -------
        reference : dict
            reference yml files
        '''
        try:
            # current folder relative
            current_folder = os.path.dirname(__file__)
            # parent folder
            parent_folder = os.path.dirname(current_folder)
            # plugin folder
            plugin_folder = os.path.join(parent_folder, 'plugin')
            # reference yml file
            reference_file = os.path.join(plugin_folder, 'reference.yml')

            # check file exists
            if os.path.exists(reference_file):
                # load yml
                with open(reference_file, 'r') as f:
                    return yaml.load(f, Loader=yaml.FullLoader)
            else:
                raise Exception('Reference file not found!')
        except Exception as e:
            raise Exception('Loading reference failed! ', e)

    def select_reference(self, name):
        '''
        Select reference yml file

        Parameters
        ----------
        name: str
            reference name

        Returns
        -------
        reference : dict
            selected reference yml file
        '''
        try:
            # reference
            reference = self.load_reference()
            # selected reference
            return reference[name]
        except Exception as e:
            raise Exception('Selecting reference failed! ', e)
