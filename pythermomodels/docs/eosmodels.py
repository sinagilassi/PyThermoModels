# EOS MODELS
# ----------

# import packages/modules
import numpy as np
from math import sqrt
# local
from ..configs import R_CONST


class EOSModels:
    # init
    def __init__(self):
        pass

    def eos_parameters(self, P, T, params, method="SRK"):
        '''
        Determine the parameters of equation of states

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        params:
            parameters:
                Pc: critical pressure [Pa]
                Tc: critical temperature [K]
                w: acentric factor [-]
                Zc: critical compressiblity factor [-]
                antoine parameters [Pa]:
                    parameters:
                    unit:
                MW: molecular weights [g/mol]
        method : str
            equation of state model, default SRK

        Returns
        -------
        res : dict
            equation of state parameters

        References
        ----------
        1. Introduction to Chemical Engineering Thermodynamics (2018)
        '''
        # get
        Pc = params['Pc']['value']
        Tc = params['Tc']['value']
        w = params['w']['value']
        Zc = params['Zc']['value']
        antoine_parameters = params['antoine-equation']['parameters']
        A = antoine_parameters[0]
        B = antoine_parameters[1]
        C = antoine_parameters[2]
        MW = params['MW']['value']

        # universal gas constant [J/mol.K]
        R = 8.314

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
        a = psi*alpha*pow(R, 2)*pow(Tc, 2)/(Pc*1e5)
        # b
        b = omega*R*Tc/(Pc*1e5)

        # beta
        beta0 = b*(P*1e5)/(R*T)
        beta1 = omega*(Pr/Tr)
        # q
        q0 = a/(b*R*T)
        q1 = psi*alpha/(omega*Tr)

        # B
        B = b*(P*1e5)/(R*T)

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
        Update the single params with mixing value of a and b

        Parameters
        ----------
        P : float
            system pressure [Pa]
        T : float
            system temperature [K]
        params : dict
            parameters
        amix : float
            molar fraction
        bmix : float
            molar fraction
        aij : float
            molar fraction
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
