# THERMODYNAMIC RELATIONS
# ------------------------

# import packages/modules
import numpy as np
import re
from math import sqrt, exp, pow, log
# internals
from ..configs import Tref, R_CONST


def extractDataCol(colName, setData, equation):
    """
    build a tuple 
    args:
        colName: column name
    """
    try:
        # choose equation
        if equation == 1:
            # heat capacity at constant pressure [kJ/kmol.K]
            CpExprData = tuple(
                [{"symbol": item['component-symbol'], "Cp": item[colName], "unit": item['unit']} for item in setData])
            # res
            return CpExprData
        elif equation == -1:
            raise Exception(
                "equation name is not found, check available equation list.")
        else:
            return -1

    except Exception as e:
        raise Exception("Extracting data failed!, ", e)


def calHeatCapacityAtConstantPressure(compList, T, params, equation):
    """ 
    calculate gas diffusivity [m2/s]
    args:
        params: changes with equation
        eq1: Chapman-Enskog (default)
    """
    # config data
    setData = params['setData']

    # choose equation
    if equation == 1:
        # heat capacity at constant pressure [kJ/kmol.K]
        CpExprData = extractDataCol("expr", setData, 1)
        # return
        return calCpEq1(compList, T, CpExprData)
    elif equation == -1:
        raise Exception(
            "equation name is not found, check available equation list.")
    else:
        return -1


def calCpEq1(comList, T, loadData):
    """
        cal: heat capacity at constant pressure
        unit: [kJ/kmol.K] 

        args:
            comList: component name list
            T: temperature [K]
            loadData: expr to be evaluated 
    """
    # try/except
    try:
        # heat capacity
        _Cpi = []

        # heat capacity data
        # loadData

        for i in comList:
            # fun expr
            CpiData = [item['Cp'] for item in loadData if i == item['symbol']]
            # fun
            def cpFun(T): return eval(CpiData[0])
            # fun exe
            CpiVal = cpFun(T)
            # store
            _Cpi.append(CpiVal)

        # convert to numpy array
        Cpi = np.array(_Cpi)
        # print("Cpi: ", Cpi)
        # res
        return Cpi
    except Exception as e:
        print(e)


def calMeanHeatCapacityAtConstantPressure(comList, params, equation, T2, T1=Tref):
    """
        cal: mean heat capacity at constant pressure 
        unit: [kJ/kmol.K] 
        args:
            comList: name of components
            T2: final temperature [K]
            T1: reference temperature [K]
    """
    # try/except
    try:
        # cp at T1 [kJ/kmol.K]
        CpT1 = calHeatCapacityAtConstantPressure(comList, T1, params, equation)
        # cp at T2 [kJ/kmol.K]
        CpT2 = calHeatCapacityAtConstantPressure(comList, T2, params, equation)
        # cp average
        CpAvg = (CpT1 + CpT2)*0.50
        # print("CpAvg: ", CpAvg)
        # res
        return CpAvg
    except Exception as e:
        print(e)
        raise


def calMixtureHeatCapacityAtConstantPressure(MoFri, HeCaCoPri):
    """
    cal: heat capacity at constant pressure of mixture 
    unit: [kJ/kmol.K] 

    args:
        MoFri: mole fraction of components 
        HeCaCoPri: heat capacity at constant pressure of components [kJ/kmol.K]
    """
    # try/except
    try:
        # check
        MoFriSize = np.size(MoFri)
        HeCaCoPriSize = np.size(HeCaCoPri)

        if MoFriSize != HeCaCoPriSize:
            raise

        # dot multiplication
        CpMix = np.dot(MoFri, HeCaCoPri)
        # res
        return CpMix
    except Exception as e:
        print(e)


def calEnthalpyChange(comList, T2, T1=Tref):
    """
        cal: enthalpy change
        unit: [kJ/kmol] 

        args:
            comList: component name list
            T2: final temperature [K]
            T1: reference temperature [K]
        return:
            dH: enthalpy change between T1 and T2
    """
    # try/except
    try:
        # cp average [kJ/kmol.K]
        CpAvg = calMeanHeatCapacityAtConstantPressure(comList, T2, T1)

        # enthalpy change [kJ/kmol]
        dH = CpAvg*(T2 - T1)
        # res
        return dH
    except Exception as e:
        print(e)


def calStandardEnthalpyOfReaction(reaExpr, standardHeatOfFormationList):
    """
        cal: standard enthalpy of reaction at 25C
        unit: [kJ/kmol]

        args:
            reaExpr: reaction expression:  
                    A + B <=> C + D

        return:
            standard heat of reaction [kJ/kmol]
    """
    # try/except
    try:
        # analyze reactions
        reaType = reaExpr.replace("<", "").replace(">", "")
        # reactant/products list
        compR = reaType.replace(r" ", "").split("=")

        # componets
        reactantList = re.findall(r"([0-9.]*)([a-zA-Z0-9.]+)", compR[0])
        productList = re.findall(r"([0-9.]*)([a-zA-Z0-9.]+)", compR[1])

        # standard heat of formation at 25
        _dHf25iReactantList = []
        _dHf25iProductList = []

        # load data [kJ/mol]
        loadData = standardHeatOfFormationList

        # reactant
        for i in reactantList:
            # fun expr
            dHf25iData = [item['dHf25']*float(i[0]) if len(i[0]) != 0 else item['dHf25']*1
                          for item in loadData if i[1] == item['symbol']]
            # store
            _dHf25iReactantList.append(dHf25iData)

        # products
        for i in productList:

            # fun expr
            dHf25iData = [item['dHf25']*float(i[0]) if len(i[0]) != 0 else item['dHf25']*1
                          for item in loadData if i[1] == item['symbol']]
            # store
            _dHf25iProductList.append(dHf25iData)

        # conversion
        dHf25iReactantList = np.array(_dHf25iReactantList).flatten()
        dHf25iProductList = np.array(_dHf25iProductList).flatten()

        # standard heat of formation at 25 [kJ/mol]
        dHf25iProductListSum = np.sum(dHf25iProductList)
        dHf25iReactantListSum = np.sum(dHf25iReactantList)

        # convert kj/mol => kJ/kmol
        dHf25 = (dHf25iProductListSum-dHf25iReactantListSum)*1000.00
        # res
        return dHf25
    except Exception as e:
        print(e)


def calHeatOfReaction(dHf25, dH):
    """
        cal: enthalpy of reaction at T
        unit:[kJ/kmol]

        args:
            dHf25: standard heat of reaction at 25C [kJ/kmol]
            dH: enthalpy of reaction [kJ/kmol]
    """
    # try/except
    try:
        # heat of reaction
        dHr = dHf25 + dH
        # res
        return dHr
    except Exception as e:
        print(e)


def calSpaceVelocity(VoFlRa, ReVo):
    """
        cal: space velocity [1/s]

        args:
            VoFlRa: volumetric flowrate [m^3/s]
            ReVo: reactor volume [m^3]

    """
    # try/except
    try:
        SpVe = VoFlRa/ReVo
        # res
        return SpVe
    except Exception as e:
        print(e)


def calGasHourlySpaceVelocity(VoFlRa, ReVo):
    """
        cal: gas hourly space velocity [1/h]

        args:
            VoFlRa: volumetric flowrate [m^3/h]
            ReVo: reactor volume [m^3]

    """
    # try/except
    try:
        GaHoSpVe = VoFlRa/ReVo
        # res
        return GaHoSpVe
    except Exception as e:
        print(e)

# NOTE


def calEnthalpyChangeOfReaction(reactionListSorted, T):
    """
    cal: standard enthalpy of reaction at 25C [kJ/kmol]
    args:
        reactionListSorted: reaction expression dict
        T: temperature [K]
    """
    # try/except
    try:
        # reaction list
        # print(f"reactionListSorted {reactionListSorted}")

        # enthalpy change list
        EnChList = []

        # reactant coefficient
        for item in reactionListSorted:
            # reactants
            _reactants = [i['symbol'] for i in item['reactants']]
            _reactantCpMeanList = calMeanHeatCapacityAtConstantPressure(
                _reactants, T)
            # reactant coeff
            _reactantCoeff = [i['coeff'] for i in item['reactants']]
            # convertion
            _loop1 = np.array(_reactantCpMeanList)
            _loop2 = np.array(_reactantCoeff)
            _loop3 = np.dot(_loop1, _loop2)

            # products
            _products = [i['symbol'] for i in item['products']]
            _productCpMeanList = calMeanHeatCapacityAtConstantPressure(
                _products, T)
            # product coeff
            _productCoeff = [i['coeff'] for i in item['products']]
            # convertion
            _loop5 = np.array(_productCpMeanList)
            _loop6 = np.array(_productCoeff)
            _loop7 = np.dot(_loop5, _loop6)

            # Cp mean of reaction
            CpMean = _loop7 + _loop3
            # print(f"CpMean: {CpMean}")

            # enthalpy change between Tref and T [kJ/kmol]
            EnChT = CpMean*(T - Tref)
            # print(f"EnChT: {EnChT}")

            # store
            EnChList.append(EnChT)

        # res
        return EnChList
    except Exception as e:
        print(e)
        raise


def calVolumetricFlowrateIG(P, T, MoFlRai):
    """
    calculate: volumetric flowrate of ideal gas (IG) [m^3/s]
    args:
        P: pressure [Pa]
        T: temperature [K]
        MoFlRai: component molar flowrate [mol/m^3]
    """
    VoFlRa = (R_CONST*T/P)*np.sum(MoFlRai)
    return VoFlRa


def calConcentrationIG(MoFlRai, VoFlRa):
    """
    calculate: concentration species species of ideal gas (IG) [mol/m^3]
    args: 
        MoFlRai: component molar flowrate [mol/m^3]
        VoFlRa: total volumetric flowrate [m^3/s]
    """
    CoSpi = MoFlRai/VoFlRa
    return CoSpi


def calDensityIG(MW, CoSp):
    """ 
    calculate: density of ideal gas (IG) [kg/m^3]
    args:
        MW: molecular weight [kg/mol]
        CoSp: concentration species [mol/m^3]
    """
    try:
        # density
        den = MW*CoSp
        return den
    except Exception as e:
        pass


def calDensityIGFromEOS(P, T, MixMW):
    """ 
    calculate: density of ideal gas (IG) [kg/m^3]
    args:
        P: pressure [Pa]
        T: temperature [K]
        MixMW: mixture molecular weight [kg/mol] 
    """
    # try/exception
    try:
        # Rg [J/kg.K]
        Rg = R_CONST/MixMW
        # density
        den = P/(Rg*T)
        return den
    except Exception as e:
        raise


def calMolarFlowRate(SpCo, SuGaVe, CrSeAr):
    """
    calculate molar flowrate
    args:
        SpCo: species concentration [kmol/m^3] | [mol/m^3]
        SuGaVe: superficial gas velocity [m/s]
        CrSeAr: cross sectional area [m^2]
    output: 
        MoFlRa: molar flowrate [kmol/s] | [mol/s]
    """
    # try/exception
    try:
        # [kmol/m^3]*[m/s]*[m^2]
        MoFlRa = SpCo*SuGaVe*CrSeAr
        return MoFlRa
    except Exception as e:
        raise

# NOTE
### viscosity ###


def calGasVisEq1(params, T):
    """ 
    gas viscosity equation 1 - Pa.s
    args:
        params: 
            equation parameters list [A,B,C,D]
        T: temperature [K]
    """
    # try/except
    try:
        A = float(params[0])
        B = float(params[1])
        C = float(params[2])
        D = float(params[3])
        _res = A*1e-6*(T**B)/(1+C*(1/T)+D*(T**-2))
        return _res
    except Exception as e:
        raise


def calGasVisEq2(eqExpr, T):
    """ 
    gas viscosity equation - Pa.s
    args:
        eqExpr: equation expression
        T: temperature [K]
    """
    # try/except
    try:
        return eval(eqExpr)
    except Exception as e:
        raise


def calGasViscosity(comList, T, loadData):
    """
        cal: gas viscosity at low pressure 
        unit: [Pa.s]
        args:
            comList: component name list
            T: temperature [K]
    """
    # try/except
    try:
        # heat capacity
        _Vii = []

        # load data

        for i in comList:
            # get id
            eqIdData = [item['id']
                        for item in loadData if i == item['component-symbol']]
            # get eq parameters
            eqData = [{"eqParams": [item['A'], item['B'], item['C'], item['D']], "eqExpr": item['expr']}
                      for item in loadData if i == item['component-symbol']]
            # check
            _eqLen = len(eqIdData) + len(eqData)
            if _eqLen > 0:
                _eqIdSet = eqIdData[0]
                _eqData = eqData[0]
                if _eqIdSet == "eq1":
                    # REVIEW
                    # eq1
                    _eqParams = _eqData.get('eqParams')
                    _res = calGasVisEq1(_eqParams, T)
                    _Vii.append(_res)
                elif _eqIdSet == "eq2":
                    # REVIEW
                    # eq2
                    _eqExpr = _eqData.get('eqExpr')
                    # build fun
                    _res = calGasVisEq2(_eqExpr, T)
                    _Vii.append(_res)
                else:
                    print('viscosity data not found, update app database!')
                    raise
            else:
                print("component not found, update the app database!")
                raise

        # convert to numpy array
        Vii = np.array(_Vii)

        # res
        return Vii
    except Exception as e:
        print(e)

# NOTE
### thermal conductivity ###


def calGasThermalConductivity(comList, T, loadData):
    """
        cal: gas thermal conductivity at low pressure 
        unit: [W/m.K]

        args:
            comList: component name list
            T: temperature [K]
    """
    # try/except
    try:
        # thermal conductivity list
        _ThCoi = []

        # load data

        for i in comList:
            # get id
            eqIdData = [item['id']
                        for item in loadData if i == item['component-symbol']]
            # get eq parameters
            eqData = [{"eqParams": [item['A'], item['B'], item['C'], item['D']], "eqExpr": item['expr']}
                      for item in loadData if i == item['component-symbol']]
            # check
            _eqLen = len(eqIdData) + len(eqData)
            if _eqLen > 0:
                _eqIdSet = eqIdData[0]
                _eqData = eqData[0]
                if _eqIdSet == "eq1":
                    # REVIEW
                    # eq1
                    _eqParams = _eqData.get('eqParams')
                    _eqExpr = _eqData.get('eqExpr')
                    _res = calGasTherCondEq1(_eqExpr, _eqParams, T)
                    _ThCoi.append(_res)
                elif _eqIdSet == "eq2":
                    # REVIEW
                    # eq2
                    _eqExpr = _eqData.get('eqExpr')
                    # build fun
                    _res = calGasVisEq2(_eqExpr, T)
                    _ThCoi.append(_res)
                else:
                    print('viscosity data not found, update app database!')
                    raise
            else:
                print("component not found, update the app database!")
                raise

        # convert to numpy array
        ThCoi = np.array(_ThCoi)

        # res
        return ThCoi
    except Exception as e:
        print(e)


def calGasTherCondEq1(expr, params, T):
    """ 
    gas thermal conductivity equation 1 - W/m.K
    args:
        params: 
            equation parameters list [C1, C2, C3, C4]
        T: temperature [K]
    """
    # try/except
    try:
        # params
        C1 = float(params[0])
        C2 = float(params[1])
        C3 = float(params[2])
        C4 = float(params[3])
        # expr
        def ThCoFun(T, C1, C2, C3, C4): return eval(expr)
        _var1 = C1*(T**C2)
        _var2 = 1 + (C3/T) + C4/(T**2)
        _res0 = _var1/_var2
        _res = ThCoFun(T, C1, C2, C3, C4)
        return _res
    except Exception as e:
        raise


def calGasTherCondEq1V1(params, T):
    """ 
    gas thermal conductivity equation 1 - W/m.K
    args:
        params: 
            equation parameters list [C1, C2, C3, C4]
        T: temperature [K]
    """
    # try/except
    try:
        C1 = params[0]
        C2 = params[1]
        C3 = params[2]
        C4 = params[3]
        _var1 = C1*(T**C2)
        _var2 = 1 + (C3/T) + C4/(T**2)
        _res = _var1/_var2
        return _res
    except Exception as e:
        raise


def calGasTherCondEq2(eqExpr, T):
    pass

# NOTE
### mixture property ###


def calMixPropertySelection(params, equation):
    '''
    select method for mixture property
    args:
        params: changes with respect of eq.
        equation:
            1: Method of Wilke (default)
    '''
    try:
        # choose equation
        if equation == 1:
            res = calMixturePropertyM1(params)
            # return
            return res
        else:
            raise Exception("the equation not found, ")

    except Exception as e:
        print(e)


def calMixturePropertyM1(params):
    '''
    calculate mixture property M1
        Method of Wilke
    args:
        compNo: component number
        Xi: property name []
        MoFri: mole fraction [-]
        MWi: molecular weight [g/mol]
    '''
    try:
        # params
        Xi = params['Xi']
        MoFri = params['MoFri']
        MWi = params['MWi']

        # component no
        compNo = len(MoFri)

        # wilke res
        wilkeCo = np.zeros((compNo, compNo))
        for i in range(compNo):
            for j in range(compNo):
                if i == j:
                    wilkeCo[i, j] = 1
                else:
                    if i < j:
                        # wilke coefficient mix
                        A = 1 + sqrt(Xi[i]/Xi[j])*((MWi[j]/MWi[i])**(1/4))
                        AA = A**2
                        B = 8*(1+(MWi[i]/MWi[j]))
                        BB = sqrt(B)
                        wilkeCo[i, j] = AA/BB
                    else:
                        C = (Xi[i]/Xi[j])*(MWi[j]/MWi[i]) * wilkeCo[j, i]
                        wilkeCo[i, j] = C
        # vars
        A = np.zeros(compNo)
        B = np.zeros((compNo, compNo))
        # mixture property
        mixProp = np.zeros(compNo)
        for i in range(compNo):
            A[i] = Xi[i]*MoFri[i]
            for j in range(compNo):
                B[i, j] = MoFri[j]*wilkeCo[i, j]
            # set
            mixProp[i] = A[i]/np.sum(B[i, :])

        mixPropVal = np.sum(mixProp)
        # res
        return mixPropVal
    except Exception as e:
        print(e)

# NOTE
# vapour pressure


def calVapourPressureEq1(params, T):
    '''
    Antoine Equation for Vapor Pressures of Pure Species

    args:
        params: Antoine equation params:
            A
            B
            C
        T: temperature [K]

    output:
        res: vapour pressure [Pa] - needs conversion from bar to Pa
    '''
    # try/except
    try:
        A = float(params[0])
        B = float(params[1])
        C = float(params[2])

        # vapour pressure [kPa]
        res = exp(A - (B/(T+C)))

        # res
        return res*1e5
    except Exception as e:
        raise


def calVapourPressureEq2(params, T):
    '''
    Antoine Equation for Vapor Pressures of Pure Species

    TABLE 2-8 Vapor Pressure of Inorganic and Organic Liquids, ln P = C1 + C2/T + C3 ln T + C4T C5, P in Pa

    args:
        params: Antoine equation params:
            C1
            C2
            C3
            C4
            C5
        T: temperature [K]

    output:
        res: vapour pressure [Pa] 
    '''
    # try/except
    try:
        C1 = float(params[0])
        C2 = float(params[1])
        C3 = float(params[2])
        C4 = float(params[3])
        C5 = float(params[4])

        # vapour pressure [kPa]
        res = exp(C1 + (C2/T) + C3*log(T) + C4*pow(T, C5))

        # res
        return res
    except Exception as e:
        raise


def calVapourPressure(comList, T, loadData):
    """
        calculate vapor pressure of pure component [Pa]

        args:
            comList: component name list (symbol) *** array ***
            T: temperature [K]
            loadData: database *** array ***
    """
    # try/except
    try:
        # res
        _Vp = []

        # load data

        for i in comList:
            # get id
            eqIdData = [item['id']
                        for item in loadData if str(i) == item['component-symbol']]
            # check
            _eqLen = len(eqIdData)
            if _eqLen > 0:
                _eqIdSet = eqIdData[0]
                # ! check
                if _eqIdSet == "eq1":
                    # * get eq parameters
                    eqData = [{"eqParams": [item['A'], item['B'], item['C']], "eqExpr": item['expr']}
                              for item in loadData if i == item['component-symbol']]
                    if len(eqData) > 0:
                        # set
                        _eqData = eqData[0]
                        # REVIEW
                        # eq1
                        _eqParams = _eqData.get('eqParams')
                        _eqExpr = _eqData.get('eqExpr')
                        _res = calVapourPressureEq1(_eqParams, T)
                        _Vp.append(_res)
                elif _eqIdSet == "eq2":
                    # * get eq parameters
                    eqData = [{"eqParams": [item['C1'], item['C2'], item['C3'], item['C4'], item['C5']]}
                              for item in loadData if i == item['component-symbol']]
                    if len(eqData) > 0:
                        # set
                        _eqData = eqData[0]
                        # REVIEW
                        # eq1
                        _eqParams = _eqData.get('eqParams')
                        _res = calVapourPressureEq2(_eqParams, T)
                        _Vp.append(_res)
                    pass
                else:
                    print('data not found, update app database!')
                    raise Exception()
            else:
                print("component not found, update the app database!")
                raise Exception()

        # convert to numpy array
        Vp = np.array(_Vp)

        # ! check
        if len(comList) == 1:
            return Vp[0]
        else:
            return Vp

    except Exception as e:
        print(e)


def calMolarVolume(P, T, Z):
    '''
    calculate molar-volume [m^3/mol]

    args:
        Z: compressibility factor [-]

    P: pressure [Pa]
    T: temperature [K]
    R: universal gas constant [J/mol.K]

    output:
        Vm: molar volume [m^3/mol]
    '''
    return Z * ((R_CONST * T) / P)


def RackettEquation(Vc, Zc, Tr):
    '''
    estimation of saturated liquid volume

    args:
        Vc: critical molar-volume
        Zc: critical compressibility factor
        Tr: reduced temperature
    '''
    _ZcPower = pow(1-Tr, 0.2857)
    return Vc*pow(Zc, _ZcPower)


def ModifiedRackettEquation(T, Pc, Tc, w):
    '''
    estimation of saturated liquid volume

    args:
        T: temperature [K]
        Pc: critical pressure [Pa]
        Tc: critical temperature [K]
        w: acentric factor
    '''
    _c0 = R_CONST*Tc/Pc
    ZRA = 0.2956 - 0.08775*w
    Tr = T/Tc
    _c1 = 1 + pow(1-Tr, 2/7)
    _c2 = _c0*pow(ZRA, _c1)

    return _c2


def calAcentricFactor():
    pass


def calVaporPressureV1():
    '''
    calculate vapor pressure using the clasusius-clapeyron equation
    '''
    pass


def calVaporPressureV2(T, Tc, w):
    '''
    calculate vapor pressure using the shortcut equation

    Note: 
        The shortcut vapor pressure equation must be regarded as an approximation for rapid estimates. 
        The  approximation  is  generally  good  above  P  =  0.5  bar; the  percent  error  can  become significant at lower pressures (and temperatures). 
        Keep in mind that its estimates are based on the critical pressure which is generally 40-50 bar and acentric factor (at Tr = 0.7).

    args:
        T: desired temperature [K]
        Pc: critical pressure [Pa]
        Tc: critical temperature [K]
        w: acentric factor [-]

    return:
        VaPe: vapor pressure [Pa]
    '''
    # reduced temperature
    Tr = T/Tc
    # vapor pressure
    _VaPe = (7/3)*(1+w)*(1-(1/Tr))
    VaPe = pow(10, _VaPe)
    # -> conversion [bar -> Pa]
    VaPe = VaPe*1e5
    # res
    return VaPe


def SetPhase(state):
    '''
    set phase 
    '''
    # define phase
    _phaseSelection = {
        "l": "liquid",
        "g": "gas",
        "s": "solid"
    }

    return _phaseSelection.get(state)


if __name__ == "__main__":
    main()
