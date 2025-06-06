# APP CONSTANTS
# --------------

# import packages/modules
import math

# ! eos names
PENG_ROBINSON = "PR"
SOAVE_REDLICH_KWONG = "SRK"
REDLICH_KWONG = "RK"
VAN_DER_WAALS = "vdW"

# ! vle models
RAOULT_MODEL = 'raoult'
MODIFIED_RAOULT_MODEL = 'modified-raoult'

# ! activity coefficient model
VAN_LAAR_ACTIVITY_MODEL = 'van-laar'
WILSON_ACTIVITY_MODEL = 'wilson'


# universal gas constant [J/mol.K]
R_CONST = 8.314472

# epsilon
EPS_CONST = 1e-30

# pi
PI_CONST = math.pi

#  STP condition
#  pressure [Pa]
Pstp = 101325
#  temperature [K]
Tstp = 273.15

# reference temperature
Tref = Tstp + 25.00

# predefined parameters
PREDEFINED_PARAMETERS = {
    "kij": []
}
