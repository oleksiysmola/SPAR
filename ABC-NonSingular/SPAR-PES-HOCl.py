import numpy as np
import matplotlib.pyplot as plt
from spar import basicFunction, readBasicFunctions, operatorMap
# import sys
# print(sys.version)

toRadians = np.pi/180
wavenumberConversion = 33.7152537836138732499014581824370 # conversion factor to cm-1
# Case of HOCl
masses = np.array([15.99059462, 1.00782503223, 34.965012622636])
qVector = np.array([0.964001600451, 1.69198245766, 102.693412447*toRadians])
numberOfModes = len(qVector)

# Read inputs
basicFunctionsInputFile: str = "basicFunctions-pot-HOCl.spar"
with open(basicFunctionsInputFile) as f:
    basicFunctionsInput = f.read()

operatorInputFile: str = "potential-HOCl.spar"
with open(operatorInputFile) as f:
    operatorInput = f.read()

basicFunctionsList: dict = readBasicFunctions(basicFunctionsInput)
basicFunctionsMainHeader = basicFunctionsInput.split("\n")[0]
containsMass: bool = False
if "mass" in basicFunctionsMainHeader:
    containsMass = True
    basicFunctionsList[0] = 1/masses # Reciprocal masses act as mode 0

operator = operatorMap(operatorInput)
# grid = np.linspace(0.65, 3.4, 1000)
# grid = np.linspace(1.2, 3.5, 1000)
grid = np.linspace(4.0*toRadians, 170.0*toRadians, 1000)
potential = np.zeros(1000)
for i in range(len(grid)):
    qVector[2] = grid[i]
    potential[i] = operator.evaluatePointOfComponent("potential", basicFunctionsList, qVector)
# for i in range(len(grid)):
#     for j in range(len(operator.componentCoefficients["potential"])):
#         potential[i] += operator.componentCoefficients["potential"][j]*basicFunctionsList[1][operator.functionIndices["potential"][j, 0]].evaluate(grid[i])*basicFunctionsList[2][operator.functionIndices["potential"][j, 1]].evaluate(qVector[1])*basicFunctionsList[3][operator.functionIndices["potential"][j, 2]].evaluate(qVector[2])

plt.plot(grid, potential)
plt.savefig("HOCl-Potential-alpha.png", dpi=300)