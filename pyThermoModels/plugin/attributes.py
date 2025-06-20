# import libs

# local

# equation of state models
EQUATION_OF_STATE_MODELS = {
    "PR": {
        'DEPENDENT_DATA': {
            'critical_temperature': {'unit': 'K', 'symbol': 'Tc'},
            'critical_pressure': {'unit': 'Pa', 'symbol': 'Pc'},
            'acentric_factor': {'unit': 'None', 'symbol': 'AcFa'},
        },
        'DEPENDENT_EQUATIONS': {
            'vapor_pressure': {'unit': 'Pa', 'symbol': 'VaPr'}
        }
    },
    "SRK": {
        'DEPENDENT_DATA': {
            'critical_temperature': {'unit': 'K', 'symbol': 'Tc'},
            'critical_pressure': {'unit': 'Pa', 'symbol': 'Pc'},
            'acentric_factor': {'unit': 'None', 'symbol': 'AcFa'},
        },
        'DEPENDENT_EQUATIONS': {
            'vapor_pressure': {'unit': 'Pa', 'symbol': 'VaPr'}
        }
    },
    "RK": {
        'DEPENDENT_DATA': {
            'critical_temperature': {'unit': 'K', 'symbol': 'Tc'},
            'critical_pressure': {'unit': 'Pa', 'symbol': 'Pc'},
            'acentric_factor': {'unit': 'None', 'symbol': 'AcFa'},
        },
        'DEPENDENT_EQUATIONS': {
            'vapor_pressure': {'unit': 'Pa', 'symbol': 'VaPr'}
        }
    },
    "vdW": {
        'DEPENDENT_DATA': {
            'critical_temperature': {'unit': 'K', 'symbol': 'Tc'},
            'critical_pressure': {'unit': 'Pa', 'symbol': 'Pc'},
        },
        'DEPENDENT_EQUATIONS': {
            'vapor_pressure': {'unit': 'Pa', 'symbol': 'VaPr'}
        }
    }
}

# activity models
ACTIVITY_MODELS = {
    "NRTL": {
        'DEPENDENT_DATA': {
            'non_randomness_parameter': {'unit': 'None', 'symbol': 'alpha'},
            'binary_interaction_parameters': {'unit': 'None', 'symbol': 'tau'},
            'interaction_energy_parameters': {'unit': 'J/mol', 'symbol': 'dg'},
            'a': {'unit': 'None', 'symbol': 'a'},
            'b': {'unit': 'None', 'symbol': 'b'},
            'c': {'unit': 'None', 'symbol': 'c'},
            'd': {'unit': 'None', 'symbol': 'd'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': "The binary-interaction-parameters (tau) can be calculated with two different methods: 1) using dg 2) using a,b,c, and d parameters."

    },
    "UNIQUAC": {
        'DEPENDENT_DATA': {
            'volume_parameter': {'unit': 'None', 'symbol': 'r'},
            'surface_area_parameter': {'unit': 'None', 'symbol': 'q'},
            'binary_interaction_parameters': {'unit': 'None', 'symbol': 'tau'},
            'interaction_energy_parameters': {'unit': 'J/mol', 'symbol': 'dU'},
            'a': {'unit': 'None', 'symbol': 'a'},
            'b': {'unit': 'None', 'symbol': 'b'},
            'c': {'unit': 'None', 'symbol': 'c'},
            'd': {'unit': 'None', 'symbol': 'd'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': "The binary-interaction-parameters (tau) can be calculated with two different methods: 1) using dU 2) using a,b,c, and d parameters."
    }
}
