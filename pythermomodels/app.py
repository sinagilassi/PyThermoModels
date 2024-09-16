# import packages/modules

# local
from .configs import __description__, __version__, packageName, \
    packageShortName
from .docs import eosCoreClass, FugacityClass, Manager


def intro():
    '''
    Package description
    '''
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def fugacity():
    '''
    Initialize fugacity calculation

    Parameters
    ----------
    None

    Returns
    -------
    fugacity: list
        fugacity for gas/liquid/solid phase
    '''
    try:
        # manager
        ManagerC = Manager()
        # return
        return ManagerC
    except Exception as e:
        raise Exception("Calculating the Fugacity failed!, ", e)
