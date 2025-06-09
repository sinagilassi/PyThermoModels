# EOS UTILS
# ----------

# import packages/modules
import numpy as np
from math import exp
from typing import List, Dict, Any, Union
import pycuc
# local
from ..configs import R_CONST


class EOSUtils:
    # init
    def __init__(self, datasource, equationsource):
        self.datasource = datasource
        self.equationsource = equationsource

    def rackett(self, Zc, Pc, Tc, T):
        '''
        Determine molar volume (assuming incompressible fluid)

        Parameters
        ----------
        Zc : float
            compressibility factor [-]
        Pc : float
            critical pressure [Pa]
        Tc : float
            critical temperature [K]
        T : float
            desired temperature [K]

        Returns
        -------
        Vsat : float
            saturated molar volume [m^3/mol]

        Note
        ----
        1. saturated molar volume is estimated using this equation within a slight percent error
        2. this equation is derived from the ideal gas law and the critical point conditions

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics, 2nd ed, 2012
        '''
        # universal gas constant [J/mol.K]
        R = R_CONST

        # critical molar volume [m^3/mol]
        Vc = Zc*R*T/Pc

        # reduced temperature
        Tr = T/Tc

        # saturated molar volume [m^3/mol]
        Vsat = Vc*pow(Zc, pow(1-Tr, 0.2857))

        # res
        return Vsat

    def poynting(self, V, Psat, P, T):
        '''
        Calculate Poynting term

        Parameters
        ----------
        V : float
            molar volume [m^3/mol]
        Psat : float
            saturated vapor pressure (at desired temperature) [Pa]
        P : float
            system pressure [Pa]
        T : float
            system temperature [K]

        Returns
        -------
        res : float
            Poynting term [m^3/mol]

        Note
        ----
        1. normally saturated molar volume is calculated for incompressible fluids

        References
        ----------
        1. Introductory Chemical Engineering Thermodynamics, 2nd ed, 2012
        '''
        # universal gas constant [J/mol.K]
        R = R_CONST
        # poynting
        res = exp(V*(P-Psat)/(R*T))
        # res
        return res

    def eos_root_analysis(self,
                          P: float,
                          T: float,
                          components: List[str],
                          tolerance: float = 1e-3,
                          **kwargs
                          ) -> Dict[str, Any]:
        '''
        Determine root numbers for each component

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        components : list[str]
            list of components
        tolerance : float
            tolerance to compare values

        Returns
        -------
        res : dict
            root analysis

        Notes
        -----
        The following rules are used to determine the number of roots for a given temperature and pressure:

        1. At T < Tc and P = Psat, 3 real roots → smallest is liquid, largest is vapor (`vapor-liquid`).
        2. At T < Tc and P > Psat, EOS may give 1 or 3 roots → use smallest (`liquid`).
        3. At T < Tc and P < Psat, EOS may give 1 or 3 roots → use largest (`vapor`).
        4. At T = Tc, one real root (critical point) → fluid is at critical state.
        5. At T > Tc, only 1 real root → fluid is supercritical (`vapor-like` or `liquid-like`).
        '''
        try:
            # NOTE: bubble and dew point
            # bubble point (liquid) and dew point (vapor) calculations are not implemented in this version
            bubble_point_pressure_mode = kwargs.get(
                'bubble_point_pressure_mode', None)
            dew_point_pressure_mode = kwargs.get(
                'dew_point_pressure_mode', None)

            # NOTE: mole-fraction
            mole_fraction = kwargs.get('mole_fraction', None)

            # NOTE: vars
            _vapor_pressure = []
            _critical_temperature = []
            _critical_pressure = []
            _root_analysis = []
            _root_no = []

            # res
            res = {}

            # phase
            phase = ""

            # set
            k = 0

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = self.equationsource[str(component)]['VaPr']
                # args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = self.check_args(
                    VaPr_args, self.datasource[str(component)])

                # build args
                _VaPr_args = self.build_args(
                    VaPr_args_required, self.datasource[str(component)])

                # NOTE: update P and T
                _VaPr_args['P'] = P
                _VaPr_args['T'] = T

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
                _vapor_pressure.append(_VaPr)

                # NOTE: data source
                # critical temperature
                _Tc_val = self.datasource[str(component)]['Tc']['value']
                _Tc_unit = self.datasource[str(component)]['Tc']['unit']
                # unit conversion
                _unit_block = f"{_Tc_unit} => K"
                _Tc = pycuc.to(_Tc_val, _unit_block)
                _critical_temperature.append(_Tc)

                # critical pressure
                _Pc_val = self.datasource[str(component)]['Pc']['value']
                _Pc_unit = self.datasource[str(component)]['Pc']['unit']
                # unit conversion
                _unit_block = f"{_Pc_unit} => Pa"
                _Pc = pycuc.to(_Pc_val, _unit_block)
                _critical_pressure.append(_Pc)

                # check
                # equality check
                pressure_equality_value = _VaPr - P
                temperature_equality_value = _Tc - T
                # check
                pressure_equality_check = abs(
                    pressure_equality_value) < tolerance
                temperature_equality_check = abs(
                    temperature_equality_value) < tolerance

                if pressure_equality_check and T < _Tc:  # ! P == _VaPr
                    _root_analysis.append(1)
                    _root_no.append("3 real roots")
                    # set phase
                    phase = "VAPOR-LIQUID"
                elif P >= _VaPr and T < _Tc:
                    _root_analysis.append(2)
                    _root_no.append("1 real root (liquid)")
                    # set phase
                    phase = "LIQUID"
                elif P <= _VaPr and T < _Tc:
                    _root_analysis.append(3)
                    _root_no.append("1 real root (vapor)")
                    # set phase
                    phase = "VAPOR"
                elif T > _Tc:
                    _root_analysis.append(4)
                    _root_no.append("1 real root (supercritical fluid)")
                    # set phase
                    phase = "SUPERCRITICAL"
                elif temperature_equality_check:  # ! T == _Tc
                    _root_analysis.append(5)
                    _root_no.append("1 real root (critical point)")
                    # set phase
                    phase = "CRITICAL"
                else:
                    raise Exception('Unknown root analysis!')

                # res
                _res = {
                    "component_name": component,
                    "pressure": P,
                    "pressure_unit": "Pa",
                    "temperature": T,
                    "temperature_unit": "K",
                    "root": _root_analysis[k],
                    "root-no": _root_no[k],
                    "phase": phase,
                    "vapor_pressure": _VaPr,
                    "vapor_pressure_unit": "Pa",
                    "critical_temperature": _Tc,
                    "critical_temperature_unit": "K",
                    "critical_pressure": _Pc,
                    "critical_pressure_unit": "Pa",
                    "tolerance": tolerance,
                    "vapor_pressure_check": pressure_equality_value,
                    "temperature_equality_value": temperature_equality_value,
                    "pressure_equality_check": pressure_equality_check,
                    "temperature_equality_check": temperature_equality_check,
                }

                # save
                res[str(component)] = _res

                # set
                k += 1

            # SECTION: calculate bubble and dew point pressure
            if bubble_point_pressure_mode is not None and dew_point_pressure_mode is not None:
                # mixture
                mixture = " | ".join(components)
                # input
                Xi = np.array(mole_fraction)
                Yi = np.array(mole_fraction)
                VaPri = np.array(_vapor_pressure)

                # NOTE: bubble point pressure calculation
                BuPoPr = float(Xi@VaPri)

                # NOTE: dew point pressure calculation
                DePoPr = float(1/(Yi@(1/VaPri)))

                # NOTE: compare with system pressure
                if P > BuPoPr:
                    # set phase
                    phase = "LIQUID"
                elif P < DePoPr:
                    # set phase
                    phase = "VAPOR"
                elif P > DePoPr and P < BuPoPr:
                    # set phase
                    phase = "VAPOR-LIQUID"
                elif abs(P - BuPoPr) < tolerance:
                    # set phase
                    phase = "Bubble Point (start of boiling)"
                elif abs(P - DePoPr) < tolerance:
                    # set phase
                    phase = "Dew Point (start of condensation)"
                else:
                    raise Exception('Unknown root analysis!')

                # res
                res_ = {
                    "component_name": mixture,
                    "pressure": P,
                    "pressure_unit": "Pa",
                    "temperature": T,
                    "temperature_unit": "K",
                    "bubble_pressure": BuPoPr,
                    "bubble_pressure_unit": "Pa",
                    "dew_point_pressure": DePoPr,
                    "dew_point_pressure_unit": "Pa",
                    "bubble_point_temperature": T,
                    "bubble_point_temperature_unit": "K",
                    "dew_point_temperature": T,
                    "dew_point_temperature_unit": "K",
                    "phase": phase,
                }

                # save
                res['mixture'] = res_

            return res
        except Exception as e:
            raise Exception('analyzing roots failed!, ', e)

    def check_args(self, args, component_datasource):
        '''
        Checks equation args

        Parameters
        ----------
        args : tuple
            equation args
        component_datasource : dict
            component datasource
        '''
        try:
            # required args
            required_args = []

            # datasource list
            datasource_component_list = list(component_datasource.keys())
            datasource_component_list.append("P")
            datasource_component_list.append("T")

            # check args within datasource
            for arg_key, arg_value in args.items():
                # symbol
                if arg_value['symbol'] in datasource_component_list:
                    # update
                    required_args.append(arg_value)
                else:
                    raise Exception('Args not in datasource!')

            # res
            return required_args

        except Exception as e:
            raise Exception('Finding args failed!, ', e)

    def build_args(self, args, component_datasource):
        '''
        Builds args

        Parameters
        ----------
        args : tuple
            equation args
        component_datasource : dict
            component datasource
        '''
        # res
        res = {}
        for arg in args:
            # symbol
            symbol = arg['symbol']
            # check not P and T
            if symbol != "P" and symbol != "T":
                # check in component database
                for key, value in component_datasource.items():
                    if symbol == key:
                        res[symbol] = value
        return res
