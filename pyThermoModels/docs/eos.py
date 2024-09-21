# eos class
# parameters for three equation of states
# van der Waals
# Redlich-Kwong and Soave
# Peng-Robinson

# import packages/modules
import numpy as np
from scipy.optimize import fsolve
import math
# internals
from ..configs import constants as CONST
from ..configs import EOS_ROOT_ACCURACY
from ..utils import roundNum, removeDuplicatesList


class EOS:
    # init
    def __init__(self, P, T, eosName, moleFraction=[]):
        self.P = P
        self.T = T
        self.eosName = eosName
        self.moleFraction = moleFraction

    # component no
    def componentNoSet(self):
        return len(self.moleFraction)

    def eos_A(self, a):
        # A value
        # var
        res = 0
        # eos name
        eosNameSet = self.eosName
        # select method
        selectEOS = {
            "PR": lambda: (a * self.P) / np.power(CONST.R_CONST * self.T, 2),
            "RK": lambda: (a * self.P) / (np.power(CONST.R_CONST, 2) * np.power(self.T, 2.5)),
            "VDW": lambda: (a * self.P) / np.power(CONST.R_CONST * self.T, 2),
        }
        # exe
        res = selectEOS.get(eosNameSet)()
        # return
        return res

    def eos_B(self, b):
        """ calculate B """
        return (b * self.P) / (CONST.R_CONST * self.T)

    def eos_Z(self, V):
        """ calculate specific volume """
        return (self.P * V) / (CONST.R_CONST * self.T)

    def eos_alpha(self, B):
        """ calculate alpha in f(Z) """
        # var
        res = 0
        # eos name
        eosNameSet = self.eosName

        # select eos
        selectEOS = {
            "VDW": lambda B: -1 - B,
            "SKR": lambda B: -1,
            "PR": lambda B: -1 + B,
        }
        # res
        res = selectEOS.get(eosNameSet)(B)
        # return
        return res

    def eos_beta(self, A, B):
        """ calculate parameter beta """
        # var
        res = 0
        # eos name
        eosNameSet = self.eosName

        # select eos
        selectEOS = {
            "VDW": lambda A, B: A,
            "SKR": lambda A, B: A - B - np.power(B, 2),
            "PR": lambda A, B: A - 3 * np.power(B, 2) - 2 * B,
        }
        # res
        res = selectEOS.get(eosNameSet)(A, B)
        # return
        return res

    def eos_gamma(self, A, B):
        """ calculate parameter gamma """
        # var
        res = 0
        # eos name
        eosNameSet = self.eosName

        # select eos
        selectEOS = {
            "VDW": lambda A, B: -A * B,
            "SKR": lambda A, B: -A * B,
            "PR": lambda A, B: -A * B + np.power(B, 2) + np.power(B, 3),
        }
        # res
        res = selectEOS.get(eosNameSet)(A, B)
        # return
        return res

    # f(Z)
    def fZ(self, x, *data):
        """ f(Z) gets parameters: alpha, beta, gamma """
        # print(data)
        alpha, beta, gamma = data
        fZSet = np.power(x, 3) + alpha*np.power(x, 2) + beta*x + gamma
        return fZSet

    # find fZ root
    def findRootfZ(self, alpha, beta, gamma):
        '''
        roots 

            1. P=P*, 3 real roots
            2. T<Tc, P>P*, 1 real root (liquid)
            3. T<Tc, P<P*, 1 real root (superheated vapor)
            4. T>Tc, 1 real root (supercritical fluid varies between `vapor-like` and `liquid-like`)
        '''
        # vars
        data = (alpha, beta, gamma)
        #
        zList = []
        zGuess = np.linspace(0, 2, 21)
        for item in zGuess:
            zLoop = fsolve(self.fZ, [item], args=data)
            zList.append(zLoop[0])

        # list -> array
        zListArray = np.array(zList)
        # find real values
        zListRealValues = np.isreal(zListArray)
        # index to value
        zListLimit = zListArray[np.where(zListRealValues == True)]
        # set accuracy
        zListLimit = self.sortRootfZ(zListLimit)
        #  remove duplicate items
        zListNet = np.array(removeDuplicatesList(zListLimit))

        # return
        return zListNet

    # net fZ root
    def sortRootfZ(self, data):
        zList = roundNum(data, EOS_ROOT_ACCURACY)
        return zList

    # k[i,j]
    def kijFill(self):
        """ 
            interaction parameter should be taken from experimental data, otherwise have to set to zero! 
        """
        # components number
        componentsNo = self.componentNoSet()
        print(f"componentsNo {componentsNo}")
        # square matrix
        matrixShape = (componentsNo, componentsNo)
        kijMatrix = np.zeros(matrixShape)
        # set
        for i in range(componentsNo):
            for j in range(componentsNo):
                if i == j:
                    kijMatrix[i, j] = 0.0
                else:
                    kijMatrix[i, j] = 0.0
        # res
        return kijMatrix

    # a[i,j]
    def aijFill(self, ai, kij):
        # components number
        componentsNo = self.componentNoSet()
        # square matrix
        matrixShape = (componentsNo, componentsNo)
        aijMatrix = np.zeros(matrixShape)
        # set
        for i in range(componentsNo):
            for j in range(componentsNo):
                aijMatrix[i, j] = (1-kij[i, j])*np.sqrt(ai[i]*ai[j])
        # res
        return aijMatrix

    # mixing rule b
    def bMixing(self, bi, xi):
        # components number
        componentsNo = self.componentNoSet()
        # square matrix
        matrixShape = (componentsNo)
        bMatrix = np.zeros(matrixShape)
        # set
        for i in range(componentsNo):
            bMatrix[i] = xi[i]*bi[i]
        # sum
        res = np.sum(bMatrix)
        return res

    # mixing rule a
    def aMixing(self, aij, xi):
        # components number
        componentsNo = self.componentNoSet()
        # square matrix
        matrixShape = (componentsNo, componentsNo)
        aMatrix = np.zeros(matrixShape)
        # set
        for i in range(componentsNo):
            for j in range(componentsNo):
                aMatrix[i, j] = xi[i]*xi[j]*aij[i, j]
        # sum
        res = np.sum(aMatrix)
        return res

    def molarVolume(self, Z):
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
        return Z * ((CONST.R_CONST * self.T) / self.P)

    def aPR(self, Pc, Tc, w):
        '''
        Calculate peng-robinson a constant

        Parameters
        ----------
        Pc : float
            critical pressure [Pa] (its unit is compatible with R [J/mol.K]/[Pa.m^3/mol.K])
        Tc : float
            critical temperature [K]
        w : float
            acentric factor [-]

        Returns
        -------
        a : float
            PR constant [Pa.(m3^2)/(mol^2)]

        Notes
        -----
        check for w
        while universal gas constant [J/mol.K]
        '''
        # check for w
        if w < 0.49:
            k = 0.37464 + 1.54226 * w - 0.26993 * (w**2)
        else:
            k = 0.379642 + 1.48503*w - 0.164423*(w**2) + 0.016666*(w**3)
        alpha = math.pow(1 + k * (1 - math.sqrt(self.T / Tc)), 2)
        a0 = (0.45723553 * (math.pow(CONST.R_CONST, 2) * math.pow(Tc, 2))) / Pc
        res = a0*alpha
        return res

    def bPR(self, Pc, Tc):
        '''
        calculate peng-robinson b constant

        args:
            Pc: critical pressure [Pa]
                its unit is compatible with R [J/mol.K]/[Pa.m^3/mol.K]
            Tc: critical temperature [K]

        while universal gas constant [J/mol.K]

        output:
            b: PR constant [m^3/mol]  
        '''
        res = (0.07779607 * CONST.R_CONST * Tc) / Pc
        return res

    @classmethod
    def aVDW(cls, Pc, Tc):
        '''
        calculate van der Waals equation a constant

        args:
            Pc: critical pressure [Pa]
                its unit is compatible with R [J/mol.K]/[Pa.m^3/mol.K]
            Tc: critical temperature [K]

        while universal gas constant [J/mol.K]

        output:
            a: constant [Pa.(m3^2)/(mol^2)]
        '''
        return (27/64)*((CONST.R_CONST**2)*(Tc**2)/(Pc))

    @classmethod
    def bVDW(cls, Pc, Tc):
        '''
        calculate van der Waals equation b constant

        args:
            Pc: critical pressure [Pa]
                its unit is compatible with R [J/mol.K]/[Pa.m^3/mol.K]
            Tc: critical temperature [K]

        while universal gas constant [J/mol.K]

        output:
            b: constant [m^3/mol]  
        '''
        return (CONST.R_CONST/8)*(Tc/Pc)

    @classmethod
    def abVDM(cls, components):
        '''
        calculate van der Waals equation a and b constants

        return:
            a:
            b:
        '''
        # component no
        compNo = len(components)

        ai = []
        bi = []

        for i in range(compNo):
            _Pc = components[i].Pc
            _Tc = components[i].Tc
            _a = cls.aVDW(_Pc, _Tc)
            _b = cls.bVDW(_Pc, _Tc)
            # save
            ai.append(_a)
            bi.append(_b)

        # res
        return ai, bi
