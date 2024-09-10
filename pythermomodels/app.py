# import packages/modules

# local
from .configs import __description__, __version__, packageName, packageShortName
from .docs import eosCoreClass, FugacityClass


def intro():
    '''
    Package description
    '''
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def calculate_fugacity(modelInput):
    '''
    Calculate fugacity for gas/liquid/solid phase

    args:
        modelInput:
            eos-name: eos equation name
            phase: component phase (gas/liquid/solid)
            compList: component list
            MoFr: mole fraction
            params: 
                pressure [Pa]
                temperature [K]
            unit: set unit (SI: default, cgs)  
    '''
    try:
        # get primary info
        compList = modelInput.get("components")
        # eos method
        eosModel = modelInput.get('eos-model')
        # phase
        phase = modelInput.get("phase")
        # ole fraction
        moleFraction = modelInput.get('MoFr', 1)
        # params
        params = modelInput.get('params')

        # check component list
        compListUnique = dUtilityClass.buildComponentList(compList)

        # load all data
        compData = loadDataEOS(compListUnique)

        # * init eos class
        _eosCoreClass = eosCoreClass(
            compData, compList, eosModel, moleFraction, params)

        # select method
        selectEOS = {
            "PR": lambda: _eosCoreClass._eosPR()
        }
        # res
        _eosRes = selectEOS.get(eosModel)()

        # * init fugacity class
        _fugacityClass = FugacityClass(compData, compList, _eosRes, phase)

        # select method
        selectFugacity = {
            "PR": lambda: _fugacityClass.FugacityPR()
        }

        # res
        _fugacityRes = selectFugacity.get(eosModel)()

        # return
        return _fugacityRes
    except Exception as e:
        raise Exception("Fugacity calculation failed!, ", e)
