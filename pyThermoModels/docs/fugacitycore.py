# FUGACITY CORE CLASS
# ------------------

# package/module list
import numpy as np
import pycuc
# local
from ..configs import Tref, R_CONST
from .eosmanager import EOSManager


class FugacityCore(EOSManager):
    '''
    Fugacity core class for calculating fugacity using equation of state (EOS) models.

    This class is a subclass of the EOSManager class and provides methods for calculating
    fugacity coefficients for different phases (vapor, liquid, solid) using various EOS models.
    '''

    def __init__(self, datasource, equationsource, components, operating_conditions, eos_parms):
        '''
        Initialize the FugacityCore class.

        Parameters
        ----------
        datasource : object
            Data source object containing thermodynamic data.
        equationsource : object
            Equation source object containing thermodynamic equations.
        components : list
            List of component names.
        operating_conditions : dict
            Dictionary containing operating conditions such as pressure and temperature.
        eos_parms : dict
            Dictionary containing EOS parameters such as phase, eos-model, mode, and liquid-fugacity-calculation-method.

        Notes
        -----
        1. pressure and temperature are converted to SI units (Pa and K).
        2. The phase can be 'VAPOR', 'LIQUID', 'VAPOR-LIQUID', 'SUPERCRITICAL', or 'SOLID'.
        3. The eos-model can be 'SRK', 'PR', or other EOS models.
        4. The mode can be 'single' or 'mixture'.
        '''
        # data/equations
        self.datasource = datasource
        self.equationsource = equationsource
        # components
        self.components = components
        # set
        self.P = pycuc.to(
            operating_conditions["pressure"][0], f"{operating_conditions['pressure'][1]} => Pa")
        # temperature
        self.T = pycuc.to(
            operating_conditions["temperature"][0], f"{operating_conditions['temperature'][1]} => K")

        # phase
        self.phase = eos_parms['phase']
        # eos model
        self.eos_model = eos_parms['eos-model']
        # mode
        self.mode = eos_parms['mode']
        # liquid fugacity calculation method
        self.liquid_fugacity_calculation_method = eos_parms['liquid-fugacity-mode']

        # comp no
        self.componentsNo = len(self.components)

        # init
        EOSManager.__init__(self, datasource, equationsource)

    def root_analysis_mode(self):
        '''
        Sets root analysis mode based on the phase set in the model input

        Returns
        -------
        calculation_mode : int
            fugacity calculation mode
        '''
        try:
            # calculation mode
            calculation_mode = -1
            # check phase
            if self.phase == 'VAPOR-LIQUID':
                calculation_mode = 1
            elif self.phase == 'LIQUID':
                calculation_mode = 2
            elif self.phase == 'VAPOR':
                calculation_mode = 3
            elif self.phase == 'SUPERCRITICAL':
                calculation_mode = 2
            elif self.phase == 'SOLID':
                calculation_mode = 2
            else:
                raise Exception('Invalid phase!')

            return calculation_mode
        except Exception as e:
            raise Exception('Setting fugacity calculation mode failed!, ', e)

    def fugacity_cal(self, yi: list, solver_method, root_analysis_set):
        '''
        Calculate fugacity

        Parameters
        ----------
        yi : list
            mole fraction of components
        solver_method : str
            solver method, default ls
        root_analysis_set : int
            root analysis set

        Returns
        -------
        res : dict
            fugacity result


        Notes
        -----
        Zi : float | list
            compressibility coefficient
        phi : float | list
            fugacity
        eos_params : dict
            eos parameters
        phi_pack : dict
            fugacity package
        '''
        try:
            # check
            if self.phase == 'VAPOR':
                res = self.gas_fugacity(
                    yi=yi, solver_method=solver_method, root_analysis_set=root_analysis_set)
            elif self.phase == 'VAPOR-LIQUID':
                res = self.gas_fugacity(
                    yi=yi, solver_method=solver_method, root_analysis_set=root_analysis_set)
            elif self.phase == 'LIQUID':
                # check
                if self.liquid_fugacity_calculation_method == 'Poynting':
                    res = self.liquid_fugacity(
                        yi=yi, solver_method=solver_method, root_analysis_set=root_analysis_set)
                elif self.liquid_fugacity_calculation_method == 'EOS':
                    res = self.gas_fugacity(
                        yi=yi, solver_method=solver_method, root_analysis_set=root_analysis_set)
            elif self.phase == 'SOLID':
                res = self.solid_fugacity()
            else:
                raise Exception('Invalid phase!')

            return res
        except Exception as e:
            raise Exception('Fugacity calculation failed!, ', e)

    def gas_fugacity(self, yi=[], **kwargs):
        '''
        Estimate gas fugacity using eos (largest Z)

        Parameters
        ----------
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
        Zis : list
            compressibility coefficient
        _phi : list
            fugacity
        _eos_params : dict
            eos parameters
        _phi_pack : dict
            fugacity package

        Notes
        -----
        1. root_analysis_set is set to 1 to check three roots
        '''
        # universal gas constant [J/mol.K]
        # R = R_CONST

        # default values
        eos_model = self.eos_model
        mode = self.mode

        # setting
        solver_method = kwargs.get('solver_method', 'ls')
        root_analysis_set = kwargs.get('root_analysis_set', 1)

        # eos model
        eos_model = self.eos_model
        # mode
        mode = self.mode

        # compressibility factor (for pure vapor phase)
        # set eos roots
        root_analysis_res = {}
        root_analysis_res['root'] = [root_analysis_set]

        # vars
        Zi = 0

        # set results
        _phi_pack = {
            'VAPOR': {},
            'LIQUID': {},
            'SOLID': {},
            'SUPERCRITICAL': {}
        }

        # check
        if mode == 'mixture':
            # Z
            _Zi, _eos_params = self.eos_roots(self.P, self.T, self.components, root_analysis_res,
                                              xi=yi, eos_model=eos_model, solver_method=solver_method, mode=mode)

            # Zi loop
            Zis = []
            # check root analysis res
            if self.phase == 'VAPOR-LIQUID':
                # ! liquid
                _Zi_set = np.min(_Zi)
                Zis.append(_Zi_set)
                # ! vapor
                _Zi_set = np.max(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'LIQUID':
                _Zi_set = np.min(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'VAPOR':
                _Zi_set = np.max(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'SUPERCRITICAL':
                # supercritical fluid
                raise Exception("supercritical fluid!")
            else:
                raise Exception('Invalid root analysis set!')

            # fugacity coefficient
            _phi = []
            # looping through Zis
            for i in range(len(Zis)):
                # fugacity coefficient (vapor phase)
                _phi_res = self.eos_fugacity(self.P, self.T, Zis[i], _eos_params, self.components,
                                             yi=yi, eos_model=eos_model, mode=mode)
                # save
                _phi.append(_phi_res)

                if len(Zis) == 1:
                    for j in range(len(_phi_res)):
                        # pack
                        _phi_pack[self.phase][self.components[j]] = {
                            'Zi': Zis[i],
                            'phi': _phi_res[j],
                            'mode': self.mode,
                            'phase': self.phase
                        }
                elif len(Zis) == 2:
                    for j in range(len(_phi_res)):
                        # pack
                        # set phase
                        if i == 0:
                            _phase = 'LIQUID'
                        else:
                            _phase = 'VAPOR'
                        # set
                        _phi_pack[_phase][self.components[j]] = {
                            'Zi': Zis[i],
                            'phi': _phi_res[j],
                            'mode': self.mode,
                            'phase': self.phase
                        }

        elif mode == 'single':
            # Z
            _Zi, _eos_params = self.eos_roots(
                self.P, self.T, self.components, root_analysis_res,
                eos_model=eos_model, solver_method=solver_method, mode=mode)

            # Zi loop
            Zis = []
            # check phase
            if self.phase == 'VAPOR-LIQUID':
                # ! liquid
                _Zi_set = np.min(_Zi)
                Zis.append(_Zi_set)
                # ! vapor
                _Zi_set = np.max(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'LIQUID':
                _Zi_set = np.min(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'VAPOR':
                _Zi_set = np.max(_Zi)
                Zis.append(_Zi_set)
            elif self.phase == 'SUPERCRITICAL':
                # supercritical fluid
                raise Exception("supercritical fluid!")
            else:
                raise Exception('Invalid root analysis set!')

            # check _Zi num
            _Zi_num = len(Zis)

            # vapor, liquid, vapor-liquid phase
            # set
            _phi = []
            # looping through _Zi
            for i in range(_Zi_num):
                # gas-liquid phase
                # fugacity coefficient (vapor phase)
                _phi_res = self.eos_fugacity(
                    self.P, self.T, Zis[i], _eos_params[0], self.components, eos_model=eos_model, mode=mode)
                # save
                _phi.append(_phi_res)
                # pack
                _phi_pack[self.phase][self.components[i]] = {
                    'Zi': Zis[i],
                    'phi': _phi_res,
                    'mode': self.mode,
                    'phase': self.phase
                }

        else:
            raise Exception("mode must be 'mixture' or 'single'")

        # res
        return Zis, _phi, _eos_params, _phi_pack

    def liquid_fugacity(self, **kwargs):
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
        # setting
        solver_method = kwargs.get('solver_method', 'ls')
        root_analysis_set = kwargs.get('root_analysis_set', 1)

        # component no.
        compNo = len(self.components)

        # thermodynamic data
        thermo_data = self.datasource
        # thermodynamic functions
        thermo_fun = self.equationsource

        # eos model
        eos_model = self.eos_model
        # mode
        mode = self.mode

        # antoine equation
        f_VaPr_equation = thermo_fun['VaPr']

        # vapor pressure [Pa]
        # at T
        VaPe = np.zeros(compNo)
        k = 0
        for i in self.components:
            # [Pa]
            _VaPe = f_VaPr_equation.cal(T=self.T)
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

    def solid_fugacity(self):
        '''
        Estimate solid fugacity
        '''
        raise Exception("Not implemented yet!")
