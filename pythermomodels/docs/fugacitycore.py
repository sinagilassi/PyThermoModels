# FUGACITY CORE CLASS
# ------------------

# package/module list
import numpy as np
# local
from ..configs import Tref, R_CONST
from .eosmanager import EOSManager


class FugacityCore(EOSManager):
    def __init__(self, compData, components, params):
        self.compData = compData
        self.components = components
        # set
        self.P = params.get("pressure", 0)
        self.T = params.get("temperature", 0)
        # comp no
        self.componentsNo = len(self.components)

        # init
        EOSManager.__init__(self, self.compData)

    def gas_fugacity(self, thermo_data, thermo_fun, yi=[], eos_model='SRK', solver_method="ls", mode="single", root_analysis_set=3):
        '''
        Estimate gas fugacity using eos (largest Z)

        Parameters
        ----------
        thermo_data : dict
            thermodynamic data
        thermo_fun : dict
            thermodynamic functions
        yi : list
            mole fraction of components
        eos_model : str
            equation of state model, default SRK
        solver_method : str
            solver method, default ls
        mode : str
            mode, default single
        root_analysis_set : int
            root analysis set

        Returns
        -------
        _phi : float
            gas fugacity
        _eos_params : dict
            eos parameters

        Notes
        -----
        1. root_analysis_set is set to 3 to check three roots 
        '''

        # universal gas constant [J/mol.K]
        # R = R_CONST

        # compressibility factor (for pure vapor phase)
        # set eos roots
        root_analysis_res = {}
        root_analysis_res['root'] = [root_analysis_set]

        # vars
        Zi = 0

        # check
        if mode == 'mixture':
            # Z
            _Zi, _eos_params = self.eos_roots(self.P, self.T, self.components, thermo_data, root_analysis_res,
                                              xi=yi, eos_model=eos_model, solver_method=solver_method, mode=mode)
            # gas phase (only max root)
            Zi = _Zi[0]

            # fugacity coefficient (vapor phase)
            _phi = self.eos_fugacity(self.P, self.T, Zi, _eos_params, self.components,
                                     thermo_data, yi=yi, eos_model=eos_model, mode=mode)

        elif mode == 'single':
            # Z
            _Zi, _eos_params = self.eos_roots(
                self.P, self.T, self.components, thermo_data, root_analysis_res, eos_model=eos_model, solver_method=solver_method, mode=mode)
            Zi = _Zi[0]

            # fugacity coefficient (vapor phase)
            _phi = self.eos_fugacity(
                self.P, self.T, Zi, _eos_params[0], self.components, thermo_data, eos_model=eos_model, mode=mode)

        else:
            raise Exception("mode must be 'mixture' or 'single'")

        # res
        return Zi, _phi, _eos_params

    def liquid_fugacity(self, thermo_data, thermo_fun, eos_model='SRK', solver_method="ls", mode="single", root_analysis_set=3):
        '''
        Estimate liquid fugacity using the Poynting term

        Parameters
        ----------
        thermo_data : dict
            thermodynamic data
        thermo_fun : dict
            thermodynamic functions
        eos_model : str
            equation of state model, default SRK
        solver_method : str
            solver method, default ls
        mode : str
            mode, default single
        root_analysis_set : int
            root analysis set

        Returns
        -------
        _phi : float
            liquid fugacity
        _eos_params : dict
            eos parameters
        '''
        # component no.
        compNo = len(self.components)

        # antoine equation
        f_antoine_equation = thermo_fun['VaPr']

        # vapor pressure [Pa]
        # at T
        VaPe = np.zeros(compNo)
        k = 0
        for i in self.components:
            # [Pa]
            _VaPe = f_antoine_equation.cal(T=self.T)
            VaPe[k] = _VaPe
            k += 1

        # compressibility factor (for pure vapor phase)
        # set eos roots
        root_analysis_res = {}
        root_analysis_res['root'] = [root_analysis_set]

        # vars
        Zi = np.zeros(compNo)
        phi = np.zeros(compNo)
        fug_l_sat = np.zeros(compNo)
        eos_params = []
        for i in range(compNo):
            _Zi, _eos_params = self.eos_roots(
                self.P, self.T, self.components, thermo_data, root_analysis_res, eos_model=eos_model, solver_method=solver_method, mode=mode)
            Zi[i] = _Zi
            eos_params.append(_eos_params[0])

            # fugacity coefficient (vapor phase)
            _phi = self.eos_fugacity(
                self.P, self.T, _Zi, _eos_params[0], self.components, thermo_data, eos_model=eos_model, mode=mode)
            phi[i] = _phi

            # fugacity of saturated vapor at T and Psat [Pa]
            _fug_l_sat = VaPe[i]*phi[i]
            fug_l_sat[i] = _fug_l_sat

        # res
        return Zi, phi, fug_l_sat, eos_params
