# import packages/modules
import numpy as np
from math import pow, sqrt, exp, log
# local
from ..configs import R_CONST


class Activity():

    def __init__(self):
        pass

    def excess_molar_gibbs_free_energy(self, xi, AcCoi):
        '''
        Calculates excess molar Gibbs free energy from component activity coefficient

        Parameters
        ----------
        xi : list
            liquid mole fraction
        AcCoi : list
            activity coefficient

        Returns
        -------
        ExGiFrEn : float
            excess molar Gibbs free energy
        '''
        try:
            # G(E)/RT
            c0 = np.zeros(self.xi)

            for i in range(self.xi):
                c0[i] = xi[i]*log(AcCoi[i])

            # excess molar gibbs energy
            ExGiFrEn = np.sum(c0)

            # res
            return ExGiFrEn
        except Exception as e:
            raise Exception("Calculating the Activity failed!, ", e)

    def Margules_activity_coefficient(self, xi, Aij):
        '''
        Determines activity coefficient for a `binary system` using the Gibbs
        excess function defined as a one and two-parameter Margules equation.

        Parameters
        ----------
        xi : float
            mole fraction
        Aij : numpy array
            Margules constant

        Returns
        -------
        AcCo : float
            activity coefficient

        Notes
        -----
        1. Aij: `Margules constant` taken from experimental results.

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics
        2. Fundamental of Chemical Engineering Thermodynamics
        '''
        try:
            # check
            if self.compNo > 2:
                return [0, 0]

            # activity coefficient
            AcCo = []

            # ! check method
            if Aij.size == 1:
                # 1-parameter
                for i in range(self.compNo):
                    AcCo.append(exp(Aij[0]*pow(1-xi[i], 2)))
            elif Aij.size == 2:
                # 2-parameter
                # set
                A12 = Aij[0]
                A21 = Aij[1]
                for i in range(self.compNo):
                    if i == 0:
                        _AcCo = exp(pow(1-xi[i], 2) *
                                    (A12 + 2*(A21 - A12)*xi[i]))
                    elif i == 1:
                        _AcCo = exp(pow(1-xi[i], 2) *
                                    (A21 + 2*(A12 - A21)*xi[i]))
                    # save
                    AcCo.append(_AcCo)
            # res
            return AcCo
        except Exception as e:
            raise Exception(
                'Calculating Margules activity coefficient failed!, ', e)

    def NRTL_activity_coefficient(self, xi, T, aij, gij):
        '''
        Determines activity coefficient for a multi-component system using Non-random two-liquid (NRTL) model

        Parameters
        ----------
        xi : list
            mole fraction
        T : float
            temperature [K]
        aij : list
            non-randomness parameter (a[i,j]=a[j,i])
        gij : list
            interaction energy parameter (calculated with temperature)

        Returns
        -------
        AcCo : list
            activity coefficient

        Notes
        -----
        1. taij: temperature dependent parameters (ta[i,i]=ta[j,j]=0)
        '''
        # component no
        compNo = xi.shape[0]

        # temperature dependent parameter
        taij = np.zeros((compNo, compNo))
        for i in range(compNo):
            for j in range(compNo):
                if j != i:
                    taij[i, j] = gij[i, j]/(R_CONST*T)

        # dependent parameters
        Gij = np.ones((compNo, compNo))
        k = 0
        for i in range(compNo):
            for j in range(compNo):
                if i != j:
                    Gij[i, j] = exp(-1*aij[i, j]*taij[i, j])
                    k += 1
                else:
                    Gij[i, j] = 1

        # activity coefficient
        AcCoi = np.zeros(compNo)

        # activity coefficient
        C0 = np.zeros((compNo, compNo))

        for i in range(compNo):
            _c0 = 0
            for j in range(compNo):
                _c0 = taij[j, i]*Gij[j, i]*xi[j] + _c0

            _c1 = 0
            for k in range(compNo):
                _c1 = Gij[k, i]*xi[k] + _c1

            for j in range(compNo):
                _c2 = xi[j]*Gij[i, j]

                _c3 = 0
                for k in range(compNo):
                    _c3 = Gij[k, j]*xi[k] + _c3

                _c4 = 0
                for n in range(compNo):
                    _c4 = xi[n]*taij[n, j]*Gij[n, j] + _c4

                _c5 = taij[i, j] - (_c4/_c3)

                # set
                C0[i, j] = (_c2/_c3)*_c5

            _c6 = (_c0/_c1) + np.sum(C0[i, :])
            AcCoi[i] = exp(_c6)

        # res
        return AcCoi
