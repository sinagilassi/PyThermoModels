# EOS MODELS
# ----------

# import packages/modules
import numpy as np
from math import pow, exp, log, sqrt
import pycuc
# local
from ..configs import R_CONST


class EOSModels:
    # init
    def __init__(self, datasource):
        self.datasource = datasource

    def eos_parameters(self, P, T, component_name, method="SRK"):
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
        1. Introduction to Chemical Engineering Thermodynamics (2018)

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
        # set component datasource
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

        # get component parameters
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
        Tc = params['Tc']['value']

        # universal gas constant [J/mol.K]
        R = R_CONST

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
                'alpha': lambda x: x**-0.50,
            },
            "SRK": {
                "sigma": 1,
                "epsilon": 0,
                "omega": 0.08664,
                'psi': 0.42748,
                'alpha': lambda Tr, omega: pow(1+(0.480+1.574*omega-0.176*(omega**2))*(1-pow(Tr, 0.5)), 2)
            },
            "PR": {
                "sigma": 1+sqrt(2),
                "epsilon": 1-sqrt(2),
                "omega": 0.07780,
                'psi': 0.45724,
                'alpha': lambda Tr, omega: pow(1+(0.37464+1.54226*omega-0.26992*(omega**2))*(1-pow(Tr, 0.5)), 2)
            },
        }

        # model parameters
        sigma = eos_params[method]['sigma']
        epsilon = eos_params[method]['epsilon']
        omega = eos_params[method]['omega']
        psi = eos_params[method]['psi']

        # Tr
        Tr = T/Tc
        # Pr
        Pr = P/Pc

        # alpha
        alpha = -1
        if method == "SRK" or method == 'PR':
            alpha = eos_params[method]['alpha'](Tr, omega)
        elif method == 'RK':
            alpha = eos_params[method]['alpha'](Tr)
        elif method == 'vdW':
            alpha = eos_params[method]['alpha']
        else:
            alpha = -1

        # a
        a = psi*alpha*pow(R, 2)*pow(Tc, 2)/(Pc)
        # b
        b = omega*R*Tc/(Pc)

        # beta
        beta0 = b*(P)/(R*T)
        beta1 = omega*(Pr/Tr)
        # q
        q0 = a/(b*R*T)
        q1 = psi*alpha/(omega*Tr)

        # B
        B = b*(P)/(R*T)

        # res
        res = {
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
            "B": B,
            "P": P,
            "T": T
        }

        # res
        return res

    def eos_parameters_mixture(self, P, T, params, amix, bmix, aij):
        '''
        Updates the single params with mixing value of a and b

        Parameters
        ----------
        P : float
            system pressure [Pa]
        T : float
            system temperature [K]
        params : dict
            parameters
        amix : float
            mixing a factor
        bmix : float
            mixing b factor
        aij : float
            mixing a[i,j]
        '''
        # universal gas constant [J/mol.K]
        R = R_CONST

        # update
        params['a'] = amix
        params['b'] = bmix

        # beta
        beta1 = bmix*(P)/(R*T)
        params['beta'] = beta1

        # q
        q1 = amix/(bmix*R*T)
        params['q'] = q1

        # aij *** new key ***
        params['aij'] = aij

        # res
        return params

    def eos_mixing_rule(self, xi, params_list, k=[]):
        '''
        Mixing rule to determine mixture a and b

        Parameters
        ----------
        xi : float
            mole fraction
        params_list : list
            list of dict of params

        Returns
        -------
        amix : float
            mixing a
        bmix : float
            mixing b
        aij : numpy array
            mixing a[i,j]
        '''
        # record no
        rNo = len(params_list)

        # ki
        if len(k) == 0:
            kij = np.zeros((rNo, rNo))
        else:
            kij = k

        # ai/bi
        ai = np.zeros(rNo)
        bi = np.zeros(rNo)

        # extract data
        for i in range(rNo):
            ai[i] = params_list[i]['a']
            bi[i] = params_list[i]['b']

        # amix
        aij = np.zeros((rNo, rNo))
        for i in range(rNo):
            for j in range(rNo):
                if i != j:
                    aij[i, j] = (1-kij[i, j])*sqrt(ai[i]*ai[j])

        _xiaij_0 = xi*aij
        _xiaij_1 = np.sum(_xiaij_0, axis=1)
        amix = np.dot(xi, _xiaij_1)

        # bmix
        bmix = np.dot(xi, bi)

        # res
        return amix, bmix, aij

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
