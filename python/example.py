# First party modules
from pyoptsparse import SLSQP, Optimization, CONMIN


# rst begin objfunc
def objfunc(xdict):
    x = xdict["xvars"]
    funcs = {}
    funcs["obj"] = -x[0] * x[1] * x[2]
    conval = [0] * 2
    conval[0] = x[0] + 2.0 * x[1] + 2.0 * x[2] - 72.0
    conval[1] = -x[0] - 2.0 * x[1] - 2.0 * x[2]
    funcs["con"] = conval
    fail = False

    return funcs, fail


# rst begin optProb
# Optimization Object
optProb = Optimization("TP037 Constraint Problem", objfunc)

# rst begin addVar
# Design Variables
optProb.addVarGroup("xvars", 3, "c", lower=[0, 0, 0], upper=[42, 42, 42], value=10)

# rst begin addCon
# Constraints
optProb.addConGroup("con", 2, lower=None, upper=0.0)

# rst begin addObj
# Objective
optProb.addObj("obj")

# rst begin print
# Check optimization problem
print(optProb)

print('---------------------')
print(' try slsqp')
print('---------------------')
optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(optProb, sens="FD") # solve
print(sol)

print('---------------------')
print(' try conmin')
print('---------------------')
optOptions = {"IPRINT": 4}
opt2 = CONMIN(options=optOptions)
sol2 = opt2(optProb, sens="FD")
print(sol2)

