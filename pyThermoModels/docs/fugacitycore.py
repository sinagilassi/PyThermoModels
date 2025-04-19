# FUGACITY CORE CLASS
# ------------------

# package/module list
import numpy as np
import pycuc
# local
from ..configs import Tref, R_CONST
from .eosmanager import EOSManager
from .eosutils import EOSUtils


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

        # root_analysis_set
        self.root_analysis_set = self.root_analysis_mode(self.phase)

        # init
        EOSManager.__init__(self, datasource, equationsource)

    @property
    def phase(self):
        '''
        Phase of the system.

        Returns
        -------
        phase : str
            phase of the system (VAPOR, LIQUID, VAPOR-LIQUID, SUPERCRITICAL, SOLID)
        '''
        return self._phase

    @phase.setter
    def phase(self, value):
        '''
        Set the phase of the system.

        Parameters
        ----------
        value : str
            phase of the system (VAPOR, LIQUID, VAPOR-LIQUID, SUPERCRITICAL, SOLID)
        '''
        if value not in ['VAPOR', 'LIQUID', 'VAPOR-LIQUID', 'SUPERCRITICAL', 'SOLID']:
            raise ValueError("Invalid phase!")
        self._phase = value

    def root_analysis_mode(self, phase: str):
        '''
        Sets root analysis mode based on the phase set in the model input

        Returns
        -------
        calculation_mode : int
            fugacity calculation mode

        Notes
        -----
        1. calculation_mode = -1: invalid phase
        2. calculation_mode = 1: VAPOR-LIQUID
        3. calculation_mode = 2: LIQUID
        4. calculation_mode = 3: VAPOR or SUPERCRITICAL
        5. calculation_mode = 4: SOLID
        '''
        try:
            # calculation mode
            calculation_mode = -1
            # check phase
            if phase == 'VAPOR-LIQUID':
                calculation_mode = 1
            elif phase == 'LIQUID':
                calculation_mode = 2
            elif phase == 'VAPOR':
                calculation_mode = 3
            elif phase == 'SUPERCRITICAL':
                calculation_mode = 1
            elif phase == 'SOLID':
                calculation_mode = 4
            else:
                raise Exception('Invalid phase!')

            return calculation_mode
        except Exception as e:
            raise Exception('Setting fugacity calculation mode failed!, ', e)

    def fugacity_cal(self, yi: list, solver_method: str):
        '''
        Calculate fugacity

        Parameters
        ----------
        yi : list
            mole fraction of components
        solver_method : str
            solver method, default ls

        Returns
        -------
        res : dict
            fugacity result

        Notes
        -----
        phi_pack : dict
            fugacity package containing fugacity coefficients and compressibility factors for each component
        '''
        try:
            # check
            if self.phase == 'VAPOR' or self.phase == 'SUPERCRITICAL':
                # SECTION: vapor
                res = self.gas_fugacity(
                    yi=yi, solver_method=solver_method,
                    root_analysis_set=self.root_analysis_set)
            elif self.phase == 'VAPOR-LIQUID':
                # SECTION: vapor-liquid
                res = self.gas_fugacity(
                    yi=yi, solver_method=solver_method,
                    root_analysis_set=self.root_analysis_set)
            elif self.phase == 'LIQUID':
                # SECTION: liquid
                # NOTE: check liquid fugacity calculation method
                if self.liquid_fugacity_calculation_method == 'Poynting':
                    res = self.liquid_fugacity(
                        yi=yi, solver_method=solver_method,
                        root_analysis_set=self.root_analysis_set)
                elif self.liquid_fugacity_calculation_method == 'EOS':
                    res = self.gas_fugacity(
                        yi=yi, solver_method=solver_method,
                        root_analysis_set=self.root_analysis_set)
            elif self.phase == 'SOLID':
                # SECTION: solid
                res = self.solid_fugacity()
            else:
                raise Exception('Invalid phase!')

            # res
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
        _phi_pack : dict
            fugacity package containing fugacity coefficients and compressibility factors for each component

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

        # SECTION: set results
        _phi_pack = {}

        # check
        if self.phase == 'VAPOR':
            _phi_pack['vapor'] = {}
        elif self.phase == 'LIQUID':
            _phi_pack['liquid'] = {}
        elif self.phase == 'SOLID':
            _phi_pack['solid'] = {}
        elif self.phase == 'SUPERCRITICAL':
            _phi_pack['supercritical'] = {}
        elif self.phase == 'VAPOR-LIQUID':
            _phi_pack['vapor'] = {}
            _phi_pack['liquid'] = {}
        else:
            raise Exception('Invalid phase!')

        # Zi loop
        Zis = []
        _Zis_comp = {}
        _phi = []
        _phi_comp = {}

        # SECTION: check mode
        if mode == 'mixture':
            # Z
            _Zi, _eos_params, _eos_params_comp = self.eos_roots(self.P, self.T, self.components, root_analysis_res,
                                                                xi=yi, eos_model=eos_model, solver_method=solver_method, mode=mode)

            # Zi loop
            Zis = []

            # check root analysis res
            if self.phase == 'VAPOR-LIQUID':
                # ! liquid
                _Zi_set = float(np.min(_Zi))
                Zis.append(_Zi_set)
                # ! vapor
                _Zi_set = float(np.max(_Zi))
                Zis.append(_Zi_set)
                # set phase
                _phi_pack['phase'] = ['liquid', 'vapor']
            elif self.phase == 'LIQUID':
                # ! liquid
                _Zi_set = float(np.min(_Zi))
                Zis.append(_Zi_set)
                # set phase
                _phi_pack['phase'] = ['liquid']
            elif self.phase == 'VAPOR':
                # ! vapor
                _Zi_set = float(np.max(_Zi))
                Zis.append(_Zi_set)
                # set phase
                _phi_pack['phase'] = ['vapor']
            elif self.phase == 'SUPERCRITICAL':
                # REVIEW: supercritical fluid
                # ! supercritical fluid
                _Zi_set = float(np.max(_Zi))
                Zis.append(_Zi_set)
                # set phase
                _phi_pack['phase'] = ['supercritical']
            elif self.phase == 'SOLID':
                raise Exception("Not implemented yet!")
            else:
                raise Exception('Invalid root analysis set!')

            # fugacity coefficient calculation
            _phi = []

            # set components
            _phi_pack['component'] = self.components

            # NOTE: looping through Zis
            for i in range(len(Zis)):
                # ! fugacity coefficient (vapor phase)
                _phi_res = self.eos_fugacity(self.P, self.T, Zis[i], _eos_params, self.components,
                                             yi=yi, eos_model=eos_model, mode=mode, params_comp=_eos_params_comp)
                # save
                _phi.append(_phi_res)

                # NOTE: check Zis
                if len(Zis) == 1:
                    for j in range(len(_phi_res)):
                        # NOTE: pack
                        _phase_key = self.phase.lower()
                        _phi_pack[_phase_key][self.components[j]] = {
                            "temperature": self.T,
                            "temperature_unit": "K",
                            "pressure": self.P,
                            "pressure_unit": "Pa",
                            'compressibility_coefficient': {
                                'value': Zis[0],
                                'unit': 'dimensionless',
                                'symbol': 'Z'
                            },
                            'fugacity_coefficient': {
                                'value': float(_phi_res[j]),
                                'unit': 'dimensionless',
                                'symbol': 'phi'
                            },
                            'mode': (self.mode).upper(),
                            'phase': self.phase
                        }

                        # NOTE: save phi_comp
                        # _phi_comp[self.components[j]] = _phi_res[j]

                    # NOTE: save Zis_comp
                    mixture_comp = "-".join(self.components)
                    _Zis_comp[mixture_comp] = [float(Zis[0])]
                elif len(Zis) == 2:  # SECTION
                    # NOTE: 2 phases (liquid and vapor)
                    for j in range(len(_phi_res)):
                        # pack
                        # set phase
                        if i == 0:
                            _phase = 'liquid'
                        else:
                            _phase = 'vapor'

                        # NOTE: set
                        _phi_pack[_phase][self.components[j]] = {
                            "temperature": self.T,
                            "temperature_unit": "K",
                            "pressure": self.P,
                            "pressure_unit": "Pa",
                            'compressibility_coefficient': {
                                'value': Zis[i],
                                'unit': 'dimensionless',
                                'symbol': 'Z'
                            },
                            'fugacity_coefficient': {
                                'value': float(_phi_res[j]),
                                'unit': 'dimensionless',
                                'symbol': 'phi'
                            },
                            'mode': (self.mode).upper(),
                            'phase': _phase
                        }

                        # NOTE: save Zis_comp
                        _Zis_comp[self.components[j]] = [Zis[i]]
                        # NOTE: save phi_comp
                        # _phi_comp[self.components[j]] = _phi_res[j]
                else:
                    raise Exception('Invalid root analysis set!')
        elif mode == 'single':  # SECTION: check mode
            # NOTE: component
            component_ = self.components[0]

            # NOTE: find roots
            _Zi, _eos_params, _eos_params_comp = self.eos_roots(
                self.P, self.T, self.components, root_analysis_res,
                eos_model=eos_model, solver_method=solver_method, mode=mode)

            # Zi loop
            Zis = []
            # phase loop
            phases = []
            # init
            _phi_pack = {}
            # set phase
            _phi_pack['phase'] = []
            # set component
            _phi_pack['component'] = component_

            # NOTE: check phase
            if self.phase == 'VAPOR-LIQUID':
                # check root note: 2 phases (liquid and vapor)
                # min and max Zi
                _Zi_set = [float(np.min(_Zi)), float(np.max(_Zi))]
                # NOTE: *** liquid ***
                Zis.append(_Zi_set[0])
                # NOTE: *** vapor ***
                Zis.append(_Zi_set[1])
                # save
                _Zis_comp[self.components[0]] = [float(x) for x in _Zi_set]
                # set phase
                phases = ['LIQUID', 'VAPOR']
            elif self.phase == 'LIQUID':
                # NOTE: *** liquid ***
                _Zi_set = float(np.min(_Zi))
                Zis.append(_Zi_set)
                # save
                _Zis_comp[self.components[0]] = [float(_Zi_set)]
                # set
                phases = ['LIQUID']
            elif self.phase == 'VAPOR':
                # NOTE: *** vapor ***
                _Zi_set = float(np.max(_Zi))
                Zis.append(_Zi_set)
                # save
                _Zis_comp[self.components[0]] = [float(_Zi_set)]
                # set
                phases = ['VAPOR']
            elif self.phase == 'SUPERCRITICAL':
                # supercritical fluid
                # NOTE: *** supercritical fluid ***
                _Zi_set = float(np.max(_Zi))
                Zis.append(_Zi_set)
                # save
                _Zis_comp[self.components[0]] = [float(_Zi_set)]
                # set
                phases = ['SUPERCRITICAL']
            elif self.phase == 'SOLID':
                raise Exception("Not implemented yet!")
            else:
                raise Exception('Invalid root analysis set!')

            # SECTION: check _Zi num
            _Zi_num = len(Zis)

            # SECTION: define phases
            # NOTE vapor, liquid, vapor-liquid, liquid-liquid phase
            # set
            _phi = []

            # looping through _Zi
            for i in range(_Zi_num):
                # NOTE: gas-liquid phase
                # ! fugacity coefficient (vapor phase)
                _phi_res = self.eos_fugacity(
                    self.P, self.T, Zis[i], _eos_params[0], self.components, eos_model=eos_model, mode=mode)
                # save
                if _phi_res is None:
                    raise Exception("fugacity coefficient calculation failed!")

                if len(_phi_res) == 1:
                    # NOTE: single phase
                    _phi.append(_phi_res[0])
                else:
                    # NOTE: multiple phase
                    _phi.append(_phi_res)
                # _phi_comp[self.components[i]] = _phi_res

                # pack
                res__ = {
                    "temperature": self.T,
                    "temperature_unit": "K",
                    "pressure": self.P,
                    "pressure_unit": "Pa",
                    'compressibility_coefficient': {
                        'value': Zis[i],
                        'unit': 'dimensionless',
                        'symbol': 'Z'
                    },
                    'fugacity_coefficient': {
                        'value': _phi_res[0],
                        'unit': 'dimensionless',
                        'symbol': 'phi'
                    },
                    'mode': (self.mode).upper(),
                    'phase': self.phase,
                    'eos_model': eos_model,
                }

                # NOTE: save
                _phase_key = phases[i].lower()
                _phi_pack[_phase_key] = res__
                _phi_pack['phase'].append(_phase_key)

        else:
            raise Exception("mode must be 'mixture' or 'single'")

        # res
        return _phi_pack

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
        try:
            # setting
            solver_method = kwargs.get('solver_method', 'ls')
            root_analysis_set = kwargs.get('root_analysis_set', 1)

            # component no.
            compNo = len(self.components)

            # NOTE: init eos utils
            EOSUtilsC = EOSUtils(self.datasource, self.equationsource)

            # eos model
            eos_model = self.eos_model
            # mode
            mode = self.mode

            # SECTION
            # method 2
            VaPr_i = np.zeros(compNo)
            # set k
            k = 0
            # loop through components
            for component in self.components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = self.equationsource[str(component)]['VaPr']
                # args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = EOSUtilsC.check_args(
                    VaPr_args, self.datasource[str(component)])

                # build args
                _VaPr_args = EOSUtilsC.build_args(
                    VaPr_args_required, self.datasource[str(component)])

                # NOTE: update P and T
                _VaPr_args['P'] = self.P
                _VaPr_args['T'] = self.T

                # NOTE: execute
                _VaPr_res = VaPr_eq.cal(**_VaPr_args)
                # extract
                _VaPr_value = _VaPr_res['value']
                _VaPr_unit = _VaPr_res['unit']
                # unit conversion
                # NOTE: unit conversion
                _unit_block = f"{_VaPr_unit} => Pa"
                _VaPr = pycuc.to(_VaPr_value, _unit_block)
                # set
                VaPr_i[k] = _VaPr
                # increment k
                k += 1

            # SECTION: compressibility factor (for pure vapor phase)
            # set eos roots
            root_analysis_res = {}
            root_analysis_res['root'] = [root_analysis_set]

            # NOTE: vars
            Zi = np.zeros(compNo)
            Zi_comp = {}
            phi = []
            phi_comp = {}
            fug_l_sat = []
            eos_params = []

            # set results
            _phi_pack = {
                'VAPOR': {},
                'LIQUID': {},
                'SOLID': {},
                'SUPERCRITICAL': {}
            }

            # SECTION: liquid fugacity calculation
            # note: At T < Tc and P > Psat, EOS may give 1 or 3 roots → use smallest (liquid).
            # looping through components
            for i in range(compNo):
                # NOTE: find roots
                _Zi, _eos_params, _eos_params_comp = self.eos_roots(
                    self.P, self.T, self.components, root_analysis_res, self.datasource, eos_model=eos_model, solver_method=solver_method, mode=mode)

                # NOTE: check root analysis res
                # set
                Zi[i] = _Zi
                Zi_comp[self.components[i]] = _Zi.tolist()
                eos_params.append(_eos_params[0])

                # ! fugacity coefficient (vapor phase)
                _phi = self.eos_fugacity(
                    self.P, self.T, _Zi, _eos_params[0], self.components, self.datasource, eos_model=eos_model, mode=mode)
                # check
                if _phi is None:
                    raise Exception("fugacity coefficient calculation failed!")
                # set
                phi.append(_phi[0])
                phi_comp[self.components[i]] = [_phi[0]]

                # NOTE: fugacity of saturated vapor at T and Psat [Pa]
                _fug_l_sat = VaPr_i[i]*phi[i]
                fug_l_sat.append(_fug_l_sat)

                # save
                # set
                _phi_pack[self.phase][self.components[i]] = {
                    "temperature": self.T,
                    "temperature_unit": "K",
                    "pressure": self.P,
                    "pressure_unit": "Pa",
                    'compressibility_coefficient': {
                        'value': Zi[i],
                        'unit': 'dimensionless',
                        'symbol': 'Z'
                    },
                    'fugacity_coefficient': {
                        'value': fug_l_sat[i],
                        'unit': 'dimensionless',
                        'symbol': 'phi'
                    },
                    'mode': (self.mode).upper(),
                    'phase': self.phase
                }

            # res
            return Zi, phi, eos_params, _phi_pack
        except Exception as e:
            raise Exception('Liquid fugacity calculation failed!, ', e)

    def solid_fugacity(self):
        '''
        Estimate solid fugacity
        '''
        raise Exception("Not implemented yet!")
