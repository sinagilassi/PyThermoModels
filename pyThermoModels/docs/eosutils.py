# EOS UTILS
# ----------

# import packages/modules
from math import exp
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

    def eos_root_analysis(self, P, T, components):
        '''
        Determine root numbers for each component

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        components : list
            list of components
        thermo_data : dict
            thermodynamic data
        thermo_fun : dict
            thermodynamic functions

        Returns
        -------
        res : dict
            root analysis

        Notes
        -----
        1. At T < Tc and P = Psat, 3 real roots → smallest is liquid, largest is vapor.
        2. At T < Tc and P > Psat, EOS may give 1 or 3 roots → use smallest (liquid).
        3. At T < Tc and P < Psat, EOS may give 1 or 3 roots → use largest (vapor).
        4. At T = Tc, one real root (critical point) → fluid is at critical state.
        5. At T > Tc, only 1 real root → fluid is supercritical (vapor-like or liquid-like).
        '''
        # vars
        _vapor_pressure = []
        _critical_temperature = []
        _critical_pressure = []
        _root_analysis = []
        _root_no = []

        # res
        res = []

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
            if P == _VaPr and T < _Tc:
                _root_analysis.append(1)
                _root_no.append("3 real roots")
            elif P >= _VaPr and T < _Tc:
                _root_analysis.append(2)
                _root_no.append("1 real root (liquid)")
            elif P <= _VaPr and T < _Tc:
                _root_analysis.append(3)
                _root_no.append("1 real root (vapor)")
            elif T > _Tc:
                _root_analysis.append(4)
                _root_no.append("1 real root (supercritical fluid)")
            elif T == _Tc:
                _root_analysis.append(5)
                _root_no.append("1 real root (critical point)")
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
                "vapor_pressure": _VaPr,
                "vapor_pressure_unit": "Pa",
                "critical_temperature": _Tc,
                "critical_temperature_unit": "K",
                "critical_pressure": _Pc,
                "critical_pressure_unit": "Pa",
            }

            # save
            res.append(_res)

            # set
            k += 1

        return res

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
