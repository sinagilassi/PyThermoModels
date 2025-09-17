# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# check version
print(ptm.__version__)
# check version
print(ptdb.__version__)
# check version
print(ptdblink.__version__)

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: thermodb directory
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# ! N2
# thermodb file name
N2_thermodb_file = os.path.join(thermodb_dir, 'nitrogen-1.pkl')
# load
N2_thermodb = ptdb.load_thermodb(N2_thermodb_file)

# ! methane
# thermodb file name
CH4_thermodb_file = os.path.join(thermodb_dir, 'methane-1.pkl')
# load
CH4_thermodb = ptdb.load_thermodb(CH4_thermodb_file)

# ! ethane
# thermodb file name
C2H6_thermodb_file = os.path.join(thermodb_dir, 'ethane-1.pkl')
# load
C2H6_thermodb = ptdb.load_thermodb(C2H6_thermodb_file)

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-1.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-1.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)

# ========================================
# NOTE: INITIALIZE PYTHERMOMODELS
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# NOTE: CHECK REFERENCES
eos_model = 'PR'  # example: 'PR', 'SRK', 'RK'
# check reference
res_ = tm.check_fugacity_reference(eos_model)
print(res_)

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('N2', N2_thermodb)
thub1.add_thermodb('CH4', CH4_thermodb)
thub1.add_thermodb('C2H6', C2H6_thermodb)
thub1.add_thermodb('n-butane', n_butane_thermodb)
thub1.add_thermodb('CO2', CO2_thermodb)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE FUGACITY FOR MULTI COMPONENT
# =======================================
# NOTE: Reference
# Example 10.7 (page 378) in Introduction to Chemical Engineering Thermodynamics
# by J.M. Smith, H.C. Van Ness, M.M. Abbott
# Example 15.2 (page 587) in Introductory Chemical Engineering Thermodynamics

# init
eos = ptm.eos()
print("eos: ", eos)

# model input
# eos model
eos_model = 'PR'

# component phase
phase = "VAPOR"

# NOTE: example 1
# feed spec
N0s = {
    'N2': 0.40,
    'CH4': 0.60,
}
# temperature [K]
T = 200
# pressure [bar]
P = 30

# NOTE: Example 9.2. (page 252), The Thermodynamics of Phase and Reaction Equilibria
# feed spec
N0s = {
    'CH4': 0.35,
    'C2H6': 0.65,
}
# # temperature [K]
T = 373.15
# # pressure [bar]
P = 30

# NOTE: Example 9.5 (page 257), The Thermodynamics of Phase and Reaction Equilibria
# eos model
eos_model = 'RK'
# feed spec
N0s = {
    'CO2': 0.15,
    'n-butane': 0.85,
}
# temperature [K]
T = 444
# pressure [bar]
P = 10
# binary interaction parameter
k_ij = [[0, 0.18],
        [0.18, 0]]

# SECTION: model input content
model_input_content = """
feed-specification: {CO2: 0.15, n-butane: 0.85}
pressure: [10, 'bar']
temperature: [444, 'K']
"""

# parse model input content
model_inputs_parsed = eos.parse_model_inputs(model_input_content)

# SECTION: model input
model_input = {
    "feed-specification": N0s,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# =======================================
# ! EOS ROOT ANALYSIS
# =======================================
# eos root analysis
res_ = eos.check_eos_roots_multi_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source
)
print(res_)

# =======================================
# ! FUGACITY CALCULATION
# =======================================
# calculate fugacity
res = eos.cal_fugacity_mixture(
    model_name=eos_model,
    model_input=model_inputs_parsed,
    model_source=model_source,
    k_ij=k_ij
)
print(res)
