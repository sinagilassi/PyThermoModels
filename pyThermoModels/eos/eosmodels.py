# EOS MODELS
# ----------
# import libs
import logging
import numpy as np
from math import pow, exp, log, sqrt
from typing import Optional, Any
import pycuc
# local
from ..configs import R_CONST, PREDEFINED_PARAMETERS

# NOTE: logger
logger = logging.getLogger(__name__)


class EOSModels():
    # init
    def __init__(
        self,
        datasource,
        equationsource,
        **kwargs
    ):
        '''
        Initialize the EOSModels class

        Parameters
        ----------
        datasource : dict
            datasource for the component parameters
        equationsource : dict
            equationsource for the equation of state parameters
        **kwargs : dict
            additional parameters for the class
        '''
        # NOTE: datasource and equationsource are dictionaries
        self.datasource = datasource
        self.equationsource = equationsource

        # NOTE: custom parameters
        self.k_ij = kwargs.get('k_ij', None)

    def eos_parameter_selection(self, method: str):
        '''
        Determine the parameters of equation of states

        Parameters
        ----------
        method : str
            equation of state model, default SRK

        Returns
        -------
        res : dict
            equation of state parameters
                1. method: str
                2. sigma: float
                3. epsilon: float
                4. omega: float
                5. psi: float
                6. alpha: float
        '''
        try:
            # NOTE: eos parameters
            # Reference: Introduction to Chemical Engineering Thermodynamics (2018)
            # Table 3.1: Parameter Assignments for Equations of State
            eos_params = {
                "vdW": {
                    "sigma": 0,
                    "epsilon": 0,
                    "omega": 0.12500,
                    'psi': 0.42188,
                    'alpha': 1,
                },
                "RK": {
                    "sigma": 1,
                    "epsilon": 0,
                    "omega": 0.08664,
                    'psi': 0.42748,
                    'alpha': lambda Tr: pow(Tr, -0.50)
                },
                "SRK": {
                    "sigma": 1,
                    "epsilon": 0,
                    "omega": 0.08664,
                    'psi': 0.42748,
                    'alpha': lambda Tr, omega: pow(
                        1 + (0.480 + 1.574 * omega - 0.176 * pow(omega, 2)) *
                        (1 - pow(Tr, 0.5)),
                        2
                    )
                },
                "PR": {
                    "sigma": 1 + sqrt(2),
                    "epsilon": 1 - sqrt(2),
                    "omega": 0.07780,
                    'psi': 0.45724,
                    'alpha': lambda Tr, omega: pow(
                        1+(0.37464 + 1.54226*omega - 0.26992*pow(omega, 2)) *
                        (1 - pow(Tr, 0.5)),
                        2
                    )
                },
            }

            # model parameters
            sigma = eos_params[method]['sigma']
            epsilon = eos_params[method]['epsilon']
            omega = eos_params[method]['omega']
            psi = eos_params[method]['psi']
            alpha = eos_params[method]['alpha']

            # save
            res = {
                "method": method,
                "sigma": sigma,
                "epsilon": epsilon,
                "omega": omega,
                "psi": psi,
                "alpha": alpha,
            }

            # res
            return res
        except Exception as e:
            raise Exception(f"Error in eos_parameter_estimation: {e}")

    def eos_parameters(
        self,
        P: float,
        T: float,
        component_name: str,
        method="SRK"
    ):
        '''
        Determine the parameters of equation of states

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        component_name : str
            name of the component
        method : str
            equation of state model, default SRK

        Returns
        -------
        res : dict
            equation of state parameters

        References
        ----------
        The method is taken from Introduction to Chemical Engineering Thermodynamics (2018), Determination of Equation-of-State Parameters (page 98)
        and Roots of the Generic Cubic Equation of State (page 99)

        Notes
        -----
        Required data for each component:
        - Pc: critical pressure [Pa]
        - Tc: critical temperature [K]
        - AcFa: acentric factor [-]
        - Zc: critical compressibility factor [-]
        - VaPr: vapor parameters [Pa]
        - MW: molecular weights [g/mol]
        '''
        # SECTION: set component datasource
        component_datasource = self.datasource.get(component_name, {})

        # check
        if component_datasource is None:
            raise Exception("component datasource not found!")

        # init
        params = {
            'Pc': {
                'value': None,
                'unit': 'Pa'
            },
            'Tc': {
                'value': None,
                'unit': 'K'
            },
        }

        # SECTION: get component parameters
        for item, value in component_datasource.items():
            # val unit checking
            _val = float(value['value'])
            _unit = value['unit']

            # ! conversion
            if item == 'Pc':
                # convert to Pa
                _unit_block = f"{_unit} => Pa"
                _val_conv = pycuc.to(_val, _unit_block)
                params['Pc']['value'] = _val_conv
            elif item == 'Tc':
                # convert to K
                _unit_block = f"{_unit} => K"
                _val_conv = pycuc.to(_val, _unit_block)
                params['Tc']['value'] = _val_conv

        # check data updated
        if params['Pc']['value'] is None or params['Tc']['value'] is None:
            raise Exception("component parameters not found!")

        # set values
        Pc = params['Pc']['value']
        Pc = float(Pc)
        Tc = params['Tc']['value']
        Tc = float(Tc)

        # NOTE: universal gas constant [J/mol.K]
        R = R_CONST

        # SECTION: EOS Parameters selection
        eos_parameter_selection_ = self.eos_parameter_selection(method)

        # NOTE: model parameters
        sigma: float = eos_parameter_selection_['sigma']
        epsilon: float = eos_parameter_selection_['epsilon']
        omega: float = eos_parameter_selection_['omega']
        psi: float = eos_parameter_selection_['psi']
        alpha_: Any = eos_parameter_selection_['alpha']

        # Tr
        Tr = T/Tc
        # Pr
        Pr = P/Pc

        # NOTE: alpha function
        # >> alpha
        alpha = 1.0

        # check type
        if method == "SRK" or method == 'PR':
            alpha: float = alpha_(Tr, omega)
        elif method == 'RK':
            alpha: float = alpha_(Tr)
        elif method == 'vdW':
            alpha: float = alpha_
        else:
            raise Exception("Unknown equation of state method!")

        # SECTION: Determination of Equation-of-State Parameters (page 98)
        # a(T)
        a = psi*alpha*pow(R, 2)*pow(Tc, 2)/Pc
        # b
        b = omega*R*Tc/Pc

        # SECTION: Roots of the Generic Cubic Equation of State (page 99)
        # beta
        beta0 = b*(P)/(R*T)
        beta1 = omega*(Pr/Tr)
        # q
        q0 = a/(b*R*T)
        q1 = psi*alpha/(omega*Tr)

        # NOTE: Chemical and engineering thermodynamics, Sandler
        # B
        B = b*(P)/(R*T)
        # A
        A = 0
        if method == 'vdW' or method == 'RK' or method == 'PR':
            A = a*(P)/pow(R*T, 2)
        elif method == 'SRK':
            A = a*(P)/(pow(R, 2)*pow(T, 2.5))
        else:
            raise Exception("Unknown equation of state method!")

        # res
        res = {
            "eos-model": method,
            "component": component_name,
            "sigma": sigma,
            "epsilon": epsilon,
            "omega": omega,
            "psi": psi,
            "Tr": Tr,
            "Pr": Pr,
            "alpha": alpha,
            "a": a,
            "b": b,
            "beta0": beta0,
            "q0": q0,
            "beta": beta1,
            "q": q1,
            "A": A,
            "B": B,
            "P": P,
            "T": T
        }

        # res
        return res

    def eos_parameters_mixture(
        self,
        P,
        T,
        amix,
        bmix,
        aij,
        A_mix,
        B_mix,
        mixture_name: str,
        eos_model: str
    ):
        '''
        Updates the single params with mixing value of a and b

        Parameters
        ----------
        P : float
            system pressure [Pa]
        T : float
            system temperature [K]
        amix : float
            mixing a factor
        bmix : float
            mixing b factor
        aij : float
            mixing a[i,j]
        A_mix : float
            mixing A factor
        B_mix : float
            mixing B factor
        mixture_name : str
            name of the mixture
        eos_model : str
            equation of state model

        Returns
        -------
        params : dict
            equation of state parameters
        '''
        try:
            # res
            params = {}

            # universal gas constant [J/mol.K]
            R = R_CONST

            # update
            params['amix'] = amix
            params['bmix'] = bmix

            # beta
            beta1 = bmix*(P)/(R*T)
            params['beta1'] = beta1

            # q
            q1 = amix/(bmix*R*T)
            params['q'] = q1

            # aij *** new key ***
            params['aij'] = aij

            # NOTE: check method
            # A
            if (
                eos_model == 'vdW' or
                eos_model == 'RK' or
                eos_model == 'PR'
            ):
                # A
                A = amix*(P)/pow(R*T, 2)
            elif eos_model == 'SRK':
                # A
                A = amix*(P)/(pow(R, 2)*pow(T, 2.5))
            else:
                raise Exception("Unknown equation of state method!")

            params['A'] = A
            params['A_mix'] = A_mix

            # B
            B = bmix*(P)/(R*T)
            params['B'] = B
            params['B_mix'] = B_mix

            # method
            params['eos-model'] = eos_model
            # mixture name
            params['component'] = mixture_name

            # NOTE: eos parameters
            params['alpha'] = self.eos_alpha(B, eos_model)
            params['beta'] = self.eos_beta(A, B, eos_model)
            params['gamma'] = self.eos_gamma(A, B, eos_model)

            # res
            return params
        except Exception as e:
            raise Exception(f"Error in eos_parameters_mixture: {e}")

    def eos_mixing_rule(
        self,
        xi,
        params_list,
        k_ij: Optional[np.ndarray | list] = None
    ):
        '''
        Mixing rule to determine mixture a and b parameters

        Parameters
        ----------
        xi : float
            mole fraction
        params_list : list
            list of dict of params
        k_ij : numpy array, optional
            2D array of binary interaction parameter (BIP), default is empty

        Returns
        -------
        a_mix : float
            mixing a
        b_mix : float
            mixing b
        aij : numpy array
            mixing a[i,j]

        Notes
        -----
        Based on Van der Waals Mixing Rules (Classical Quadratic Mixing Rules)

        - aij, bij are the pure component attraction parameters
        - kj is the binary interaction parameter (BIP), which accounts for deviations from ideal mixing
        '''
        # record no
        rNo = len(params_list)

        # ki
        if k_ij is None:
            k_ij = np.zeros((rNo, rNo))
        else:
            # check
            if isinstance(k_ij, list):
                k_ij = np.array(k_ij)

        # ai,bi, Ai,Bi
        ai = np.zeros(rNo)
        bi = np.zeros(rNo)
        Ai = np.zeros(rNo)
        Bi = np.zeros(rNo)

        # extract data
        for i in range(rNo):
            ai[i] = params_list[i]['a']
            bi[i] = params_list[i]['b']
            Ai[i] = params_list[i]['A']
            Bi[i] = params_list[i]['B']

        # NOTE: Attraction parameter amix
        a_ij = self.__aij(ai, k_ij)
        A_ij = self.__aij(Ai, k_ij)

        # NOTE: Calculate a_mix
        a_mix = 0.0
        A_mix = 0.0

        # looping through the matrix
        for i in range(rNo):
            for j in range(rNo):
                a_mix += xi[i] * xi[j] * a_ij[i, j]
                A_mix += xi[i] * xi[j] * A_ij[i, j]

        # NOTE: Covolume parameter
        # bmix
        b_mix = np.dot(xi, bi)
        # Bmix
        B_mix = np.dot(xi, Bi)

        # res
        return a_mix, b_mix, a_ij, A_mix, B_mix

    def __aij(
        self,
        ai: np.ndarray,
        k_ij: np.ndarray
    ):
        '''
        calculate aij for mixture using Van der Waals mixing rules

        Parameters
        ----------
        ai : numpy array
            1D array of pure component attraction parameter (a), default is empty
        k_ij : numpy array
            2D array of binary interaction parameter (BIP), default is empty

        Returns
        -------
        aij : numpy array
            mixing aij[i,j]
        '''
        # record no
        rNo = len(ai)

        # NOTE: Attraction parameter
        # aij
        aij = np.zeros((rNo, rNo))

        # looping through the matrix
        for i in range(rNo):
            for j in range(rNo):
                aij[i, j] = (1 - k_ij[i, j])*sqrt(ai[i]*ai[j])

        # res
        return aij

    def eos_alpha(self, B, eosNameSet):
        """ calculate alpha in f(Z) """
        try:
            # select eos
            selectEOS = {
                "VDW": lambda B: -1 - B,
                "SRK": lambda B: -1,
                "RK": lambda B: -1,
                "PR": lambda B: -1 + B,
            }

            # >> check eosNameSet
            if eosNameSet not in selectEOS.keys():
                logger.error(f"eos_alpha: Unknown eosNameSet {eosNameSet}")
                raise

            # res
            res = selectEOS[eosNameSet](B)
            # return
            return res
        except Exception as e:
            raise Exception(f"Error in eos_alpha: {e}")

    def eos_beta(self, A, B, eosNameSet):
        """ calculate parameter beta """
        try:
            # select eos
            selectEOS = {
                "VDW": lambda A, B: A,
                "SRK": lambda A, B: A - B - np.power(B, 2),
                "RK": lambda A, B: A - B - np.power(B, 2),
                "PR": lambda A, B: A - 3 * np.power(B, 2) - 2 * B,
            }

            # >> check eosNameSet
            if eosNameSet not in selectEOS.keys():
                logger.error(f"eos_beta: Unknown eosNameSet {eosNameSet}")
                raise

            # res
            res = selectEOS[eosNameSet](A, B)
            # return
            return res
        except Exception as e:
            raise Exception(f"Error in eos_beta: {e}")

    def eos_gamma(self, A, B, eosNameSet):
        """ calculate parameter gamma """
        try:
            # select eos
            selectEOS = {
                "VDW": lambda A, B: -A * B,
                "SRK": lambda A, B: -A * B,
                "RK": lambda A, B: -A * B,
                "PR": lambda A, B: -A * B + np.power(B, 2) + np.power(B, 3),
            }
            # >> check eosNameSet
            if eosNameSet not in selectEOS.keys():
                logger.error(f"eos_gamma: Unknown eosNameSet {eosNameSet}")
                raise

            # res
            res = selectEOS[eosNameSet](A, B)
            # return
            return res
        except Exception as e:
            raise Exception(f"Error in eos_gamma: {e}")

    def eos_equation(self, x, params):
        '''
        Build a polynomial 3rd degree

        Parameters
        ----------
        x : float
            variable
        params : dict
            parameters

        Returns
        -------
        fZ : float
            function

        Notes
        -----
        params:
        - sigma:
        - epsilon:
        - omega:
        - psi:
        - Tr:
        - Pr:
        - alpha:
        - a:
        - b:
        - beta:
        - q:
        - P:
        - T:

        References
        ----------
        1. Introduction to Chemical Engineering Thermodynamics
        '''
        # model parameters
        sigma = params['sigma']
        epsilon = params['epsilon']
        omega = params['omega']
        psi = params['psi']
        beta = params['beta']
        q = params['q']

        # function coefficient
        a0 = 1
        a1 = (sigma+epsilon)*beta - (1+beta)
        a2 = beta*(q + epsilon*sigma*beta - (1+beta)*(sigma+epsilon))
        a3 = (beta**2)*(q + (1+beta)*epsilon*sigma)

        fZ = a0*(x**3) + a1*(x**2) + a2*(x) - a3

        return fZ

    def eos_equation_coefficient(self, params):
        '''
        Build a list of coefficients for the cubic equation

        Parameters
        ----------
        params : dict
            parameters

        Returns
        -------
        list
            coefficients of the cubic equation

        Notes
        -----
        params:
        - sigma:
        - epsilon:
        - omega:
        - psi:
        - Tr:
        - Pr:
        - alpha:
        - a:
        - b:
        - beta:
        - q:
        - P:
        - T:
        '''
        # model parameters
        sigma = params['sigma']
        epsilon = params['epsilon']
        beta = params['beta']
        q = params['q']

        # function coefficient
        a0 = 1
        a1 = (sigma+epsilon)*beta - (1+beta)
        a2 = beta*(q + epsilon*sigma*beta - (1+beta)*(sigma+epsilon))
        a3 = (beta**2)*(q + (1+beta)*epsilon*sigma)

        return [a0, a1, a2, a3]

    def eos_equation_mixture(self, x, params):
        """
        Build a polynomial 3rd degree

        Parameters
        ----------
        x : float
            variable
        params : dict
            parameters

        Returns
        -------
        fZ : float
            function
        """
        # print(data)
        alpha, beta, gamma = params['alpha'], params['beta'], params['gamma']

        # set
        fZSet = x**3 + alpha*(x**2) + beta*x + gamma
        return fZSet

    def eos_equation_coefficient_mixture(self, params):
        """
        Build a list of coefficients for the cubic equation

        Parameters
        ----------
        params : dict
            parameters

        Returns
        -------
        list
            coefficients of the cubic equation
        """
        # print(data)
        alpha, beta, gamma = params['alpha'], params['beta'], params['gamma']

        return [1, alpha, beta, gamma]
