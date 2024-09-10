# import packages/modules

# local
from .config import __description__, __version__, packageName, packageShortName


def intro():
    '''
    Package description
    '''
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)
