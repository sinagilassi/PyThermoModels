# EOS UTILS
# ----------

# import packages/modules
from math import exp
# local
from ..configs import R_CONST


class EOSUtils:
    # init
    def __init__(self):
        pass

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

    def eos_root_analysis(self, P, T, components, thermo_data, thermo_fun):
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
        1. P=P*, 3 real roots
        2. T<Tc, P>P*, 1 real root (liquid)
        3. T<Tc, P<P*, 1 real root (superheated vapor)
        4. T>Tc, 1 real root (supercritical fluid varies between `vapor-like` and `liquid-like`)
        '''
        # vars
        _vapor_pressure = []
        _critical_temperature = []
        _critical_pressure = []
        _root_analysis = []
        _root_no = []

        # antoine equation
        f_antoine_equation = thermo_fun['antoine-equation']['fun']

        for i in components:
            # antoine equations [bar]
            _antoine_parameters = thermo_data[str(
                i)]['antoine-equation']['parameters']
            _vp = f_antoine_equation(_antoine_parameters, T)
            _vapor_pressure.append(_vp)

            # critical temperature
            _Tc = thermo_data[str(i)]['Tc']['value']
            _critical_temperature.append(_Tc)

            # critical pressure
            _Pc = thermo_data[str(i)]['Pc']['value']
            _critical_pressure.append(_Pc)

            if P == _vp and T < _Tc:
                _root_analysis.append(1)
                _root_no.append("3 real roots")
            elif P >= _vp and T < _Tc:
                _root_analysis.append(2)
                _root_no.append("1 real root (liquid)")
            elif P <= _vp and T < _Tc:
                _root_analysis.append(3)
                _root_no.append("1 real root (vapor)")
            elif T > _Tc:
                _root_analysis.append(4)
                _root_no.append("1 real root (supercritical fluid)")

        # res
        res = {
            "components": components[0],
            "P": P,
            "T": T,
            "root": _root_analysis,
            "root-no": _root_no,
            "VaPr": _vapor_pressure,
            "Tc": _critical_temperature,
            "Pc": _critical_pressure
        }

        return res
