# import packages/modules
import pythermomodels as ptm

# check version
print(ptm.__version__)

# =======================================
# CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# model input
# eos model
eos_model = 'PR'

# component phase
phase = "gas"

# component list
comp_list = ["H2O"]
# required component input
# MW,Tc,Pc,w,Zc,Vc

# mole fraction
MoFri = []

# temperature [K]
T = 200 + 273.15

# pressure [Pa]
P = 15.55*1e5

# model input
modelInput = {
    "eos-model": eos_model,
    "phase": phase,
    "components": comp_list,
    "MoFri": MoFri,
    "params": {
        "P": P,
        "T": T,
    },
    "unit": "SI",
}

# eos
res = ptm.calculate_fugacity(modelInput)
# log
print("res: ", res)
