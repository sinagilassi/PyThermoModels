# EOS SOLVER CLASS
# ----------------

# import packages/modules
import numpy as np
from scipy import optimize
from math import exp, log
from typing import List, Dict, Literal, Callable
# local
from .eosmodels import EOSModels
from ..configs import R_CONST


class EOSManager(EOSModels):

    def __init__(self, datasource, equationsource, **kwargs):
        '''
        Initialize the EOSManager class.

        Parameters
        ----------
        datasource : object
            Data source object containing thermodynamic data.
        equationsource : object
            Equation source object containing equation of state data.
        **kwargs : dict, optional
            Additional keyword arguments.
        '''
        # SECTION: init variables
        self.datasource = datasource
        self.equationsource = equationsource

        # SECTION: custom parameters
        # k_ij, 2D array of binary interaction parameter (BIP)
        self.k_ij = kwargs.get('k_ij', None)

        # NOTE: init
        EOSModels.__init__(self, datasource, equationsource, **kwargs)

    def __call__(self):
        pass

    def eos_roots(
        self,
        P: float,
        T: float,
        components: List[str],
        root_analysis: dict,
        xi=[],
        eos_model: str = "SRK",
        solver_method: str = "ls",
        mode: str = "single",
        **kwargs
    ):
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
        # SECTION: SRK params
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

        # SECTION: root analysis (case no:1,2,3,4)
        # NOTE: root selection as:
        # 1. 3 roots (vapor-liquid)
        # 2. 1 root (liquid)
        # 3. 1 root (vapor) and (supercritical fluid)
        # 4. 1 root (solid)
        _root = root_analysis['root']

        # SECTION: check root analysis
        if mode == 'mixture':
            # mixture name
            mixture_name = " | ".join(components)

            # mole fraction
            xi = np.array(xi)

            # mixture a and b
            amix, bmix, aij, A_mix, B_mix = self.eos_mixing_rule(
                xi,
                _eos_params,
                k_ij=self.k_ij
            )

            # new params *** mixture ***
            _params_mixture = self.eos_parameters_mixture(
                P,
                T,
                amix,
                bmix,
                aij,
                A_mix,
                B_mix,
                mixture_name,
                eos_model
            )
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
        k = 0

        # SECTION: set functions & coefficients
        if mode == 'single':
            # main function
            fZ = self.eos_equation
            # coefficients
            fZ_coeff = self.eos_equation_coefficient(_eos_params_0)
        elif mode == 'mixture':
            fZ = self.eos_equation_mixture
            fZ_coeff = self.eos_equation_coefficient_mixture(_eos_params_0)
        else:
            raise Exception("mode must be 'single' or 'mixture'")

        # NOTE: fZ prime and second
        fpZ = None  # eos_equation_prime
        fp2Z = None  # eos_equation_prime2

        # SECTION: *** least-square ***
        if solver_method == 'ls':
            # check eos params
            if _eos_params_0 is None:
                raise Exception("eos_params is None!")

            if not isinstance(_eos_params_0, dict):
                raise Exception("eos_params must be a dictionary!")

            # ! root finding
            Z = self.root_ls(
                root_id=_root_0,
                fZ=fZ,
                eos_params=_eos_params_0
            )

            # check res
            if len(Z) == 0:
                raise Exception("root analysis failed!")

        # SECTION: *** newton method ***
        elif solver_method == 'newton':
            # initial guess
            zGuess, steps = np.linspace(1e-5, 2, guess_no, retstep=True)

            # loop through initial guess
            for k, item in enumerate(zGuess):
                _x0 = item
                # ! root finding
                _res = optimize.newton(
                    fZ,
                    _x0,
                    fprime=fpZ,
                    fprime2=fp2Z,
                    args=(_eos_params_0,)
                )

                zList[k] = _res

            # NOTE: root analysis
            # Filter real and positive roots, and remove duplicates
            Z_ = np.array([z for z in zList if np.isreal(z) and z > 0])
            # Remove duplicates (with tolerance for floating point)
            if len(Z_) > 0:
                # Use .real to handle complex numbers safely
                Z = np.unique(np.round(Z_.real, 5))
            else:
                raise Exception("No valid roots found in the Newton method.")

        # SECTION: *** fsolve method ***
        elif solver_method == 'fsolve':
            # initial guess
            zGuess, steps = np.linspace(1e-5, 2, guess_no, retstep=True)

            # loop through initial guess
            for k, item in enumerate(zGuess):
                _x0 = item
                # ! root finding
                _res = optimize.fsolve(
                    fZ,
                    _x0,
                    args=(_eos_params_0,)
                )
                # result checks
                if _res:
                    if len(_res) == 1:
                        zList[k] = _res[0]

            # NOTE: root analysis
            # Filter real and positive roots, and remove duplicates
            Z_ = np.array([z for z in zList if np.isreal(z) and z > 0])
            # Remove duplicates (with tolerance for floating point)
            if len(Z_) > 0:
                # Use .real to handle complex numbers safely
                Z = np.unique(np.round(Z_.real, 5))
            else:
                raise Exception("No valid roots found in the fsolve method.")

        # SECTION: *** root method ***
        elif solver_method == 'root':
            # ! root finding
            _res = np.roots(fZ_coeff)
            # save
            _res = list(_res[_res > 0])
            _res = [float(i) for i in _res]

            # root analysis
            Z = self.eos_root_analysis(int(_root_0), _res)
        else:
            raise Exception(
                "solver_method must be 'ls', 'newton', 'fsolve' or 'root'")

        # res
        return np.array(Z), _eos_params, _eos_params_comp

    def root_ls(
        self,
        root_id: Literal[1, 2, 3, 4],
        fZ: Callable,
        eos_params: Dict,
        guess_no: int = 50,
        ftol=1e-8,
        xtol=1e-8,
        **kwargs
    ):
        '''
        Finds the roots of a function using the least-squares method.

        Parameters
        ----------
        root_id : int
            The ID of the root analysis performed defined as follows:
                - 1: 3 roots (vapor-liquid)
                - 2: 1 root (liquid)
                - 3: 1 root (vapor)
                - 4: 1 root (superheat)
        fZ : function
            The function for which to find the roots.
        eos_params : list
            The parameters for the equation of state.
        guess_no : int, optional
            The number of initial guesses for the roots. Default is 50.
        ftol : float, optional
            Tolerance for termination by the change of the cost function. Default is 1e-8.
        xtol : float, optional
            Tolerance for termination by the change of the independent variables. Default is 1e-8.
        **kwargs : dict, optional
            Additional parameters for the optimization process.
            - bounds: a list of minimum, medium, and maximum bounds for the guesses, default list is [-2.0, 0.5, 5.0]

        Returns
        -------
        res : OptimizeResult
            The result of the optimization.
        '''
        try:
            # SECTION: initialize variables
            # z list
            zList = np.zeros(guess_no)
            # compressibility factor
            Z = []
            # bound
            bound = guess_no*[(0, 0)]
            # cost
            fZ_cost = np.zeros(guess_no)
            # set index
            k = 0

            # SECTION: set bound
            bounds = kwargs.get('bounds', None)

            # check
            if bounds is not None:
                # set bounds
                b_min = bounds[0]
                b_max = bounds[1]
                b_mid = bounds[2]
            else:
                b_min = -2.0
                b_max = 5.0
                b_mid = 0.5

            # SECTION: define space
            # NOTE: check root analysis
            if root_id == 1:  # ! 3 roots (vapor-liquid)
                # initial guess
                zGuess, steps = np.linspace(
                    b_min, b_max, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    # check both halves
                    if item < 0.5:
                        bound[k] = (b_min, b_mid)
                        k += 1
                    else:
                        bound[k] = (b_mid, b_max)
                        k += 1
            elif root_id == 2:  # ! 1 root (liquid)
                # initial guess
                zGuess, steps = np.linspace(
                    b_min, b_mid, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    bound[k] = (b_min, b_mid)
                    k += 1
            elif root_id == 3:  # ! 1 root (vapor)
                # initial guess
                zGuess, steps = np.linspace(
                    b_mid, b_max, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    bound[k] = (b_mid, b_max)
                    k += 1
            elif root_id == 4:  # ! 1 root (superheat)
                # initial guess
                zGuess, steps = np.linspace(
                    b_min, b_max, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    bound[k] = (b_min, b_max)
                    k += 1
            else:
                raise Exception("root analysis failed!")

            # SECTION: find Z
            # set
            k = 0

            # loop through initial guess
            for item in zGuess:
                # NOTE: set initial guess
                _x0 = item

                # least-square
                _res = optimize.least_squares(
                    fZ,
                    _x0,
                    args=(eos_params,),
                    bounds=bound[k],
                    ftol=ftol,
                    xtol=xtol
                )

                # NOTE: result checks
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
            _res = zList[zList > 0]
            _res = [float(i) for i in _res]

            # root analysis
            Z = self.eos_root_analysis(int(root_id), _res)

            # res
            return Z
        except Exception as e:
            raise Exception(f"Error in find_rooting_ls: {e}") from e

    def eos_root_analysis(
        self,
        root_id: int,
        Zi: List[float]
    ) -> List[float]:
        '''
        Determines the number of roots based on the root ID and Z values.
        The root ID indicates the type of root analysis performed, and the Z values are the results of the analysis.

        Parameters
        ----------
        root_id : int
            The ID of the root analysis performed.
        Zi : list
            The Z values obtained from the root analysis.

        Returns
        -------
        Z : list
            The roots based on the root ID and Z values.
        '''
        try:

            # res
            Z = []

            # SECTION: check roots
            if root_id == 1:  # ! 3 roots
                # min/max
                _Z_min = min(Zi)
                Z.append(_Z_min)
                _Z_max = max(Zi)
                Z.append(_Z_max)
            elif root_id == 2:  # ! 1 root (liquid)
                # min
                _Z_min = min(Zi)
                Z.append(_Z_min)
            elif root_id == 3:  # ! 1 root (vapor)
                # max
                _Z_max = max(Zi)
                Z.append(_Z_max)
            elif root_id == 4:  # ! 1 root (superheat)
                # max
                _Z_max = max(Zi)
                Z.append(_Z_max)
            else:
                raise Exception("root analysis failed!")

            # res
            return Z
        except Exception as e:
            raise Exception(f"Error in eos_root_analysis: {e}") from e

    def eos_fugacity(
        self,
        P: float,
        T: float,
        Z,
        params,
        components: List,
        yi=[],
        eos_model: str = "SRK",
        mode: str = "single",
        **kwargs
    ):
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
        # params_comp = kwargs.get('params_comp', None)

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

    def SRK(
        self,
        P: float,
        T: float,
        Z: float,
        params,
        components: List,
        yi
    ):
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

            # mole fraction
            yi = np.array(yi)

            # fugacity coefficient
            phi = []

            # SECTION: mixing parameters
            _params_mix = params[-1]
            # mix parameters
            A = _params_mix['A']
            B = _params_mix['B']
            A_mix = _params_mix['A_mix']
            B_mix = _params_mix['B_mix']
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

    def PR(
        self,
        P: float,
        T: float,
        Z: float,
        params,
        components: List,
        yi
    ):
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

            # mole fraction
            yi = np.array(yi)

            # fugacity coefficient
            phi = []

            # SECTION: mixing parameters
            _params_mix = params[-1]
            # mix parameters
            A = _params_mix['A']
            B = _params_mix['B']
            A_mix = _params_mix['A_mix']
            B_mix = _params_mix['B_mix']
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

    def RK(
        self,
        P: float,
        T: float,
        Z: float,
        params,
        components: List,
        yi
    ):
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

    def vdW(
        self,
        P: float,
        T: float,
        Z: float,
        params,
        components: List,
        yi
    ):
        """
        Calculate fugacity coefficients for each component in a vapor mixture using the vdW EOS.

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
                ln_phi = (bi[i]/(Vm - b)) - np.log((Z - B)/Z) - \
                    (2 * sum_aij) / (R * T * Vm)

                res_ = np.exp(ln_phi)
                phi.append(res_)

            return phi
        except Exception as e:
            raise Exception(f"Error in pr_fugacity_coefficients: {e}") from e
