# EOS SOLVER CLASS
# ----------------

# import packages/modules
import numpy as np
from scipy import optimize
# local
from .eosmodels import EOSModels


class EOSManager(EOSModels):

    def __init__(self, datasource):
        self.datasource = datasource
        # init
        EOSModels.__init__(self, datasource)

    def __call__(self):
        pass

    def eos_roots(self, P, T, components, thermo_data, root_analysis, xi=[], eos_model="SRK", solver_method="ls", mode="single"):
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

        Notes
        -----
        1. P=P*, 3 real roots
        2. T<Tc, P>P*, 1 real root (liquid)
        3. T<Tc, P<P*, 1 real root (superheated vapor)
        4. T>Tc, 1 real root (supercritical fluid varies between `vapor-like` and `liquid-like`)
        '''
        # SRK params
        _eos_params = []

        for component in components:
            # build eos params
            _params = self.eos_parameters(P, T, component, method=eos_model)
            _params['eos-model'] = eos_model
            _eos_params.append(_params)

        # root analysis (case no:1,2,3,4)
        _root = root_analysis['root']

        if mode == 'mixture':
            # mole fraction
            xi = np.array(xi)

            # mixture a and b
            amix, bmix, aij = self.eos_mixing_rule(xi, _eos_params)
            # log
            # print(aij)

            # new params *** mixture ***
            _params_mixture = self.eos_parameters_mixture(
                P, T, _eos_params[0], amix, bmix, aij)
            _eos_params.append(_params_mixture)

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
        _x0_min = 1e-6
        k = 0

        # functions
        fZ = self.eos_equation  # SRK_equation
        # fpZ = eos_equation_prime
        # fp2Z = eos_equation_prime2

        # ! method selection
        if solver_method == 'ls':  # *** least-square ***

            if _root_0 == 1:  # 3 roots
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

            elif _root_0 == 2:  # 1 root (liquid)
                # initial guess
                zGuess, steps = np.linspace(-2, 0.5, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (-2, 0.5)
                    k += 1

            elif _root_0 == 3:  # 1 root (vapor)
                # initial guess
                zGuess, steps = np.linspace(0.5, 2, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (0.5, 2)
                    k += 1

            elif _root_0 == 4:  # 1 root (superheat)
                # initial guess
                zGuess, steps = np.linspace(-2, 2, guess_no, retstep=True)

                k = 0
                for item in zGuess:
                    _bound[k] = (-2, 2)
                    k += 1

            # set
            k = 0

            # find Z
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
                    zList[k] = _res.x
                    # cost/fun
                    fZ_cost[k] = _res.cost

                # set
                k += 1

            # res analysis
            _zList = zList[zList > 0]
            # min/max
            _Z_min = np.min(_zList)
            _Z_max = max(_zList)

            if _root_0 == 1:  # 3 roots
                Z.append(_Z_min)
                Z.append(_Z_max)
            elif _root_0 == 2:  # 1 root (liquid)
                Z.append(_Z_min)
            elif _root_0 == 3:  # 1 root (vapor)
                Z.append(_Z_max)
            elif _root_0 == 4:  # 1 root (superheat)
                Z.append(_Z_max)

        elif solver_method == 'newton':  # *** newton method ***
            # initial guess
            zGuess, steps = np.linspace(1e-5, 2, guess_no, retstep=True)

            k = 0
            for item in zGuess:
                _x0 = item
                _res = optimize.newton(fZ, _x0, fprime=fpZ,
                                       fprime2=fp2Z, args=(_eos_params_0,))
                zList[k] = _res
                k += 1

        elif solver_method == 'fsolve':  # *** fsolve method ***
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

        return np.array(Z), _eos_params
