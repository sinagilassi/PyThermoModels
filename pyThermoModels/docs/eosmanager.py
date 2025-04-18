# EOS SOLVER CLASS
# ----------------

# import packages/modules
import numpy as np
from scipy import optimize
from math import pow, exp, log, sqrt
from typing import List, Dict
# local
from .eosmodels import EOSModels
from ..configs import R_CONST


class EOSManager(EOSModels):

    def __init__(self, datasource, equationsource):
        self.datasource = datasource
        self.equationsource = equationsource
        # init
        EOSModels.__init__(self, datasource, equationsource)

    def __call__(self):
        pass

    def eos_roots(self, P: float, T: float,
                  components: List[str],
                  root_analysis: int,
                  xi=[], eos_model: str = "SRK",
                  solver_method: str = "ls",
                  mode: str = "single"):
        '''
        Estimates fugacity coefficient at fixed temperature and pressure through finding Z (lowest: for liquid, largest: for vapor)

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        components : list
            list of components
        root_analysis : dict
            root analysis
        xi : list
            mole fraction of components
        eos_model : str
            equation of state model, default SRK
        solver_method : str
            solver method, default ls
        mode : str
            mode, default single

        Returns
        -------
        Z : float
            fugacity coefficient
        _eos_params : dict
            equation of state parameters
        _eos_params_comp : dict
            equation of state parameters for each component

        Notes
        -----
        ### Roots:
            1. P=P*, 3 real roots
            2. T<Tc, P>P*, 1 real root (liquid)
            3. T<Tc, P<P*, 1 real root (superheated vapor)
            4. T>Tc, 1 real root (supercritical fluid varies between `vapor-like` and `liquid-like`)

        ### Hints
            1. _eos_params updated for single and mixture
        '''
        # SRK params
        _eos_params = []
        _eos_params_comp = {}

        # SECTION: eos parameters for each component
        for component in components:
            # build eos params
            _params = self.eos_parameters(P, T, component, method=eos_model)

            # save eos params
            _eos_params.append(_params)
            # update component params
            _eos_params_comp[component] = _params

        # NOTE: root analysis (case no:1,2,3,4)
        _root = root_analysis['root']

        # SECTION: check root analysis
        if mode == 'mixture':
            # mixture name
            mixture_name = " | ".join(components)

            # mole fraction
            xi = np.array(xi)

            # mixture a and b
            amix, bmix, aij = self.eos_mixing_rule(xi, _eos_params)

            # new params *** mixture ***
            _params_mixture = self.eos_parameters_mixture(
                P, T, amix, bmix, aij, mixture_name, eos_model)
            # set
            _eos_params.append(_params_mixture)
            _eos_params_comp['mixture'] = _params_mixture

        # check mode
        if mode == 'single':
            # SRK params
            _eos_params_0 = _eos_params[0]
        elif mode == 'mixture':
            _eos_params_0 = _eos_params[-1]
        else:
            _eos_params_0 = []

        # root status
        _root_0 = _root[0]

        # find fZ roots
        guess_no = 50
        zList = np.zeros(guess_no)
        Z = []
        _bound = guess_no*[(0, 0)]
        fZ_cost = np.zeros(guess_no)
        k = 0

        # SECTION: functions
        if mode == 'single':
            fZ = self.eos_equation  # SRK_equation
        elif mode == 'mixture':
            fZ = self.eos_equation_mixture
        else:
            raise Exception("mode must be 'single' or 'mixture'")

        # NOTE: fZ prime and second
        fpZ = None  # eos_equation_prime
        fp2Z = None  # eos_equation_prime2

        # SECTION: *** least-square ***
        if solver_method == 'ls':
            # NOTE: check root analysis
            if _root_0 == 1:  # ! 3 roots (vapor-liquid)
                # initial guess
                zGuess, steps = np.linspace(-2, 2, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    # check both halves
                    if item < 0.5:
                        _bound[k] = (-2, 0.5)
                        k += 1
                    else:
                        _bound[k] = (0.5, 2)
                        k += 1

            elif _root_0 == 2:  # ! 1 root (liquid)
                # initial guess
                zGuess, steps = np.linspace(-2, 0.5, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (-2, 0.5)
                    k += 1

            elif _root_0 == 3:  # ! 1 root (vapor)
                # initial guess
                zGuess, steps = np.linspace(0.5, 2, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (0.5, 2)
                    k += 1

            elif _root_0 == 4:  # ! 1 root (superheat)
                # initial guess
                zGuess, steps = np.linspace(-2, 2, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (-2, 2)
                    k += 1

            else:
                raise Exception("root analysis failed!")

            # set
            k = 0

            # SECTION: find Z
            for item in zGuess:
                _x0 = item
                # least-square
                _res = optimize.least_squares(fZ, _x0, args=(
                    _eos_params_0,), bounds=_bound[k], ftol=1e-8, xtol=1e-8)
                # roots
                # _res = optimize.root(fZ, _x0, args=(_eos_params_0,))

                # result checks
                if (_res.success is True and _res.x > 0):
                    # x
                    # check
                    if len(_res.x) == 1:
                        zList[k] = _res.x[0]
                    else:
                        raise Exception("least_squares failed!")
                    # cost/fun
                    fZ_cost[k] = _res.cost

                # set
                k += 1

            # res analysis
            _zList = zList[zList > 0]
            # min/max
            _Z_min = np.min(_zList)
            _Z_max = np.max(_zList)

            # SECTION: check roots
            if _root_0 == 1:  # ! 3 roots
                Z.append(_Z_min)
                Z.append(_Z_max)
            elif _root_0 == 2:  # ! 1 root (liquid)
                Z.append(_Z_min)
            elif _root_0 == 3:  # ! 1 root (vapor)
                Z.append(_Z_max)
            elif _root_0 == 4:  # ! 1 root (superheat)
                Z.append(_Z_max)
            else:
                raise Exception("root analysis failed!")

        # SECTION: *** newton method ***
        elif solver_method == 'newton':
            # initial guess
            zGuess, steps = np.linspace(1e-5, 2, guess_no, retstep=True)

            k = 0
            for item in zGuess:
                _x0 = item
                _res = optimize.newton(fZ, _x0, fprime=fpZ,
                                       fprime2=fp2Z, args=(_eos_params_0,))
                zList[k] = _res
                k += 1

        # SECTION: *** fsolve method ***
        elif solver_method == 'fsolve':
            # initial guess
            zGuess, steps = np.linspace(1e-5, 2, guess_no, retstep=True)

            k = 0
            for item in zGuess:
                _x0 = item
                _res = optimize.fsolve(fZ, _x0, args=(_eos_params_0,))
                # result checks
                if _res:
                    zList[k] = _res
                # set
                k += 1

        # res
        return np.array(Z), _eos_params, _eos_params_comp

    def eos_fugacity(self, P: float, T: float, Z,
                     params,
                     components: List,
                     yi=[],
                     eos_model: str = "SRK",
                     mode: str = "single",
                     **kwargs):
        '''
        Determines fugacity coefficient

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        Z : list
            compressibility factor
        params : list
            equation of state parameters
        components : list
            list of components
        yi : list
            mole fraction of components
        eos_model : str
            equation of state model, default SRK
        mode : str
            mode, default single
        **kwargs : dict
            additional parameters
            - params_comp: dict

        Returns
        -------
        fugacity : list
            fugacity coefficient
        '''
        # NOTE: load params component
        params_comp = kwargs.get('params_comp', None)

        # NOTE: universal gas constant [J/mol.K]
        R = R_CONST

        # NOTE: molar volume [m^3/mol]
        V = Z*R*T/P

        # SECTION: single
        if mode == 'single':

            # model parameters
            sigma = params['sigma']
            epsilon = params['epsilon']
            omega = params['omega']
            psi = params['psi']
            alpha = params['alpha']
            beta = params['beta']
            q = params['q']
            a = params['a']
            b = params['b']
            B = params['B']

            # model selection
            if eos_model == "vdW":
                _phi = (Z-1) - log(Z*(1-(b/V))) - (a/(R*T*V))
            elif eos_model == "RK":
                _phi = (Z-1) - log(Z*(1-(b/V))) - (a/(b*R*T))*log(1+(b/V))
            elif eos_model == "SRK" or eos_model == "PR":
                _phi = (Z-1) - log(Z-B) - (alpha/(b*R*T*(sigma-epsilon))) * \
                    log((Z+sigma*B)/(Z+epsilon*B))
            else:
                raise Exception("eos_model must be 'vdW', 'RK', 'SRK' or 'PR'")

            phi = [exp(_phi)]

        # # SECTION: mixture
        elif mode == 'mixture':
            # NOTE: check model
            if eos_model == "SRK":
                # fugacity coefficient
                phi = self.SRK(P, T, Z, params, components, yi)
            elif eos_model == "RK":
                # fugacity coefficient
                phi = self.RK(P, T, Z, params, components, yi)
            elif eos_model == "PR":
                # fugacity coefficient
                phi = self.PR(P, T, Z, params, components, yi)
            else:
                raise Exception(f"{eos_model} not available!")
        else:
            raise Exception("mode must be 'single' or 'mixture'")

        # res
        return phi

    def SRK(self,
            P: float,
            T: float,
            Z: float,
            params,
            components: List,
            yi):
        """
        Calculate fugacity coefficients for each component in a vapor mixture using the SRK EOS.

        Parameters
        ----------
        P : float
            Pressure of the system [Pa].
        T : float
            Temperature of the system [K].
        Z : float
            Compressibility factor of the system.
        params : list
            List of dictionaries containing the EOS parameters for each component.
        components : list
            List of component names.
        yi : list
            Mole fractions of each component in the vapor phase.

        Returns
        -------
        phi : list
            Fugacity coefficients for each component in the vapor phase.
        """
        try:
            # component number
            N = len(components)

            # fugacity coefficient
            phi = []

            # SECTION: mixing parameters
            _params_mix = params[-1]
            # mix parameters
            A = _params_mix['A']
            B = _params_mix['B']
            # required parameters
            a = _params_mix['amix']
            b = _params_mix['bmix']
            aij = _params_mix['aij']  # array

            # SECTION: component parameters
            ai = []
            bi = []
            # looping through components
            for i in range(N):
                # select component parameters
                # NOTE: load params component
                _params_i = params[i]
                # required parameters
                ai_ = _params_i['a']
                bi_ = _params_i['b']
                # set
                ai.append(ai_)
                bi.append(bi_)

            for i in range(N):
                # fugacity coefficient
                sum_aij = yi@aij[i, :]

                # term 2
                ln_phi = (bi[i] / b) * (Z - 1) - np.log(Z - B) \
                    - (A / B) * (2 * sum_aij / a -
                                 bi[i] / b) * np.log(1 + B / Z)

                res_ = np.exp(ln_phi)
                phi.append(res_)

            return phi
        except Exception as e:
            raise Exception(f"Error in srk_fugacity_coefficients: {e}") from e

    def PR(self,
            P: float,
            T: float,
            Z: float,
            params,
            components: List,
            yi):
        """
        Calculate fugacity coefficients for each component in a vapor mixture using the PR EOS.

        Parameters
        ----------
        P : float
            Pressure of the system [Pa].
        T : float
            Temperature of the system [K].
        Z : float
            Compressibility factor of the system.
        params : list
            List of dictionaries containing the EOS parameters for each component.
        components : list
            List of component names.
        yi : list
            Mole fractions of each component in the vapor phase.

        Returns
        -------
        phi : list
            Fugacity coefficients for each component in the vapor phase.
        """
        try:
            # component number
            N = len(components)

            # fugacity coefficient
            phi = []

            # SECTION: mixing parameters
            _params_mix = params[-1]
            # mix parameters
            A = _params_mix['A']
            B = _params_mix['B']
            # required parameters
            a = _params_mix['amix']
            b = _params_mix['bmix']
            aij = _params_mix['aij']  # array

            # SECTION: component parameters
            ai = []
            bi = []
            # looping through components
            for i in range(N):
                # select component parameters
                # NOTE: load params component
                _params_i = params[i]
                # required parameters
                ai_ = _params_i['a']
                bi_ = _params_i['b']
                # set
                ai.append(ai_)
                bi.append(bi_)

            for i in range(N):
                # fugacity coefficient
                sum_aij = yi@aij[i, :]

                # term 2
                ln_phi = (bi[i]/b)*(Z - 1) - np.log(Z - B) \
                    - (A / (2 * np.sqrt(2) * B)) * (2 * sum_aij / a - bi[i]/b) \
                    * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))

                res_ = np.exp(ln_phi)
                phi.append(res_)

            return phi
        except Exception as e:
            raise Exception(f"Error in pr_fugacity_coefficients: {e}") from e

    def RK(self,
            P: float,
            T: float,
            Z: float,
            params,
            components: List,
            yi):
        """
        Calculate fugacity coefficients for each component in a vapor mixture using the RK EOS.

        Parameters
        ----------
        P : float
            Pressure of the system [Pa].
        T : float
            Temperature of the system [K].
        Z : float
            Compressibility factor of the system.
        params : list
            List of dictionaries containing the EOS parameters for each component.
        components : list
            List of component names.
        yi : list
            Mole fractions of each component in the vapor phase.

        Returns
        -------
        phi : list
            Fugacity coefficients for each component in the vapor phase.
        """
        try:
            # NOTE: universal gas constant [J/mol.K]
            R = R_CONST
            # NOTE: molar volume [m^3/mol]
            Vm = Z*R*T/P

            # component number
            N = len(components)

            # fugacity coefficient
            phi = []

            # SECTION: mixing parameters
            _params_mix = params[-1]
            # mix parameters
            A = _params_mix['A']
            B = _params_mix['B']
            # required parameters
            a = _params_mix['amix']
            b = _params_mix['bmix']
            aij = _params_mix['aij']  # array

            # SECTION: component parameters
            ai = []
            bi = []
            # looping through components
            for i in range(N):
                # select component parameters
                # NOTE: load params component
                _params_i = params[i]
                # required parameters
                ai_ = _params_i['a']
                bi_ = _params_i['b']
                # set
                ai.append(ai_)
                bi.append(bi_)

            for i in range(N):
                # fugacity coefficient
                sum_aij = yi@aij[i, :]

                # term 2
                ln_phi = (bi[i]/b)*(Z - 1) - np.log(Z - B) \
                    - (a / (b * R * T**1.5)) * (2 * sum_aij / a - bi[i]/b) \
                    * np.log((Vm + b) / Vm)

                res_ = np.exp(ln_phi)
                phi.append(res_)

            return phi
        except Exception as e:
            raise Exception(f"Error in pr_fugacity_coefficients: {e}") from e
