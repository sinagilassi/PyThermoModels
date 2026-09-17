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
    },
    "ENRTL": {
        'DEPENDENT_DATA': {
            'non_randomness_parameter': {'unit': 'None', 'symbol': 'alpha'},
            'binary_interaction_parameters': {'unit': 'None', 'symbol': 'tau'},
            'interaction_energy_parameters': {'unit': 'J/mol', 'symbol': 'dg'},
            'ion_size': {'unit': 'documented by long-range basis', 'symbol': 'a_i'},
            'pitzer_debye_huckel_A_phi': {'unit': 'basis-dependent', 'symbol': 'A_phi'},
            'pitzer_debye_huckel_b': {'unit': 'basis-dependent', 'symbol': 'b'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Electrolyte NRTL activity model for true-species electrolyte "
            "states. The first formulation identifier is chen_evans_1986; "
            "charge is consumed from component metadata and speciation is "
            "kept outside ENRTL."
        )
    },
    "PITZER": {
        'DEPENDENT_DATA': {
            'beta0': {'unit': 'kg/mol', 'symbol': 'beta0'},
            'beta1': {'unit': 'kg/mol', 'symbol': 'beta1'},
            'c_phi': {'unit': 'kg^2/mol^2', 'symbol': 'C_phi'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Pitzer v1: one fully dissociated binary electrolyte on a molality "
            "basis, using the single-alpha formulation and reporting gamma_pm only."
        )
    },
    "WILSON": {
        'DEPENDENT_DATA': {
            'lambda_parameters': {'unit': 'None', 'symbol': 'lambda'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Wilson activity model using positive dimensionless Lambda "
            "parameters. The standard Wilson model cannot represent "
            "liquid-liquid splitting."
        )
    },
    "MARGULES": {
        'DEPENDENT_DATA': {
            'A': {'unit': 'None', 'symbol': 'A'},
            'A12': {'unit': 'None', 'symbol': 'A12'},
            'A21': {'unit': 'None', 'symbol': 'A21'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Binary two-suffix and three-suffix Margules activity model "
            "with dimensionless excess-Gibbs parameters."
        )
    },
    "VAN_LAAR": {
        'DEPENDENT_DATA': {
            'a1': {'unit': 'None', 'symbol': 'a1'},
            'a2': {'unit': 'None', 'symbol': 'a2'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Binary van Laar activity model with dimensionless a1 and a2 "
            "parameters."
        )
    },
    "REDLICH_KISTER": {
        'DEPENDENT_DATA': {
            'a': {'unit': 'None', 'symbol': 'a'},
        },
        'DEPENDENT_EQUATIONS': {},
        'DESCRIPTION': (
            "Binary Redlich-Kister excess-Gibbs expansion using ordered "
            "dimensionless coefficients a[0..n]."
        )
    }
}
