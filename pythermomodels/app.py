# import packages/modules

# local
from .configs import __description__, __version__, packageName, \
    packageShortName
from .docs import eosCoreClass, FugacityClass, ManagerClass


def intro():
    '''
    Package description
    '''
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def calculate_fugacity(input_file):
    '''
    Calculate fugacity for gas/liquid/solid phases

    Parameters
    ----------
    input_file: str
        input_file file path

    Returns
    -------
    fugacity: list
        fugacity for gas/liquid/solid phase
    '''
    try:
        # load input file
        model_input = ManagerClass.load_yml(input_file)

        fugacity = 1

        # return
        return fugacity
    except Exception as e:
        raise Exception("Calculating the Fugacity failed!, ", e)
