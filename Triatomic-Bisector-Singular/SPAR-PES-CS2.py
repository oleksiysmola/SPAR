import numpy as np
import matplotlib.pyplot as plt
from spar import basicFunction, readBasicFunctions, operatorMap
# import sys
# print(sys.version)

toRadians = np.pi/180
wavenumberConversion = 33.7152537836138732499014581824370 # conversion factor to cm-1
# Case of CS2
masses = np.array([12.00000000, 31.9720707])
qVector = np.array([1.5521, 1.5521, np.pi])
numberOfModes = len(qVector)

# Read inputs
basicFunctionsInputFile: str = "basicFunctions-pot-CS2-AMES.basic"
with open(basicFunctionsInputFile) as f:
    basicFunctionsInput = f.read()

operatorInputFile: str = "potential-CS2-AMES.spar"
with open(operatorInputFile) as f:
    operatorInput = f.read()

basicFunctionsList: dict = readBasicFunctions(basicFunctionsInput)
basicFunctionsMainHeader = basicFunctionsInput.split("\n")[0]
containsMass: bool = False
if "mass" in basicFunctionsMainHeader:
    containsMass = True
    basicFunctionsList[0] = 1/masses # Reciprocal masses act as mode 0

operator = operatorMap(operatorInput)

# for key in basicFunctionsList[1]:
#     print(key, basicFunctionsList[1][key].evaluate(1.55275000))

## COMPARE TO AMES IMPLEMENTATION
energiesGridFile: str = "CS2-AMES-PotEnergy.dat"
with open(energiesGridFile) as f:
    energiesGrid = f.read()

energiesGridLines: list = energiesGrid.split("\n")
numberOfGridPoints: int = len(energiesGridLines)
AMESgrid = np.zeros((numberOfGridPoints, 3))
AMESenergies = np.zeros(numberOfGridPoints)
energiesSPAR = np.zeros(numberOfGridPoints)
for i in range(numberOfGridPoints):
    energiesGridLineSplit = energiesGridLines[i].split()
    AMESenergies[i] = float(energiesGridLineSplit[numberOfModes])
    for j in range(numberOfModes):
        AMESgrid[i, j] = float(energiesGridLineSplit[j])
    AMESgrid[i, 2] *= toRadians
    energiesSPAR[i] = operator.evaluatePointOfComponent("potential", basicFunctionsList, AMESgrid[i])

SPARenergiesFile: str = "CS2-AMES-PotEnergy-SPAR.dat"
with open(SPARenergiesFile, "w") as f:
    for i in range(numberOfGridPoints):
        f.write(f"{AMESgrid[i, 0]}  {AMESgrid[i, 1]}  {AMESgrid[i, 2]/toRadians}   {AMESenergies[i]}     {energiesSPAR[i]}  {AMESenergies[i] - energiesSPAR[i]}\n")

AMESerror = AMESenergies - energiesSPAR
print(max(AMESerror))
# print("Contains mass: ", operator.containsMass)
# print(basicFunctionsList[1][0].evaluate(1.5521))
# print(basicFunctionsList[1][0].evaluate(2.0))
# testQ = 2.0
# testAlpha = 150*toRadians
# print(basicFunctionsList[1][0].evaluate(testQ))
# print(np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[1][1].evaluate(testQ))
# print((testQ - 1.5521)*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[1][4].evaluate(testQ))
# print((testQ - 1.5521)**4*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[1][12].evaluate(testQ))
# print(basicFunctionsList[2][12].evaluate(testQ))
# print((testQ - 1.5521)**12*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[3][2].evaluate(testAlpha))
# print((np.cos(testAlpha) - np.cos(np.pi))**2*np.exp(-0.5*min(-abs(np.pi-testQ) - (5*np.pi/6-np.pi), 0)**2-1.5*min(-abs(np.pi-testQ) - (5*np.pi/6-np.pi), 0)**4))
# print("Die katze ist sehr nett:")
# print(basicFunctionsList[3][0].evaluate(testAlpha))
# print(np.exp(-0.5*min(-abs(np.pi-testAlpha) - (5*np.pi/6-np.pi), 0)**2-1.5*min(-abs(np.pi-testAlpha) - (5*np.pi/6-np.pi), 0)**4))
# print(basicFunctionsList[1][8].evaluate(testQ))
# print((testQ - 1.5521)**8*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[2][2].evaluate(testQ))
# print((testQ - 1.5521)**2*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))
# print(basicFunctionsList[1][8].evaluate(testQ)*basicFunctionsList[2][2].evaluate(testQ)*basicFunctionsList[3][2].evaluate(testAlpha))
# print(((testQ - 1.5521)**8*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))*((testQ - 1.5521)**2*np.exp(-0.200*(testQ - 1.5521)**2-0.200*(testQ - 1.5521)**4))*((np.cos(testAlpha) - np.cos(np.pi))**2*np.exp(-0.5*min(-abs(np.pi-testAlpha) - (5*np.pi/6-np.pi), 0)**2-1.5*min(-abs(np.pi-testAlpha) - (5*np.pi/6-np.pi), 0)**4)))
# qVector = [testQ, testQ, testAlpha]
# print("armes wurstchen:")
# operator.evaluatePointOfComponent("potential", basicFunctionsList, qVector)
# print(1.200000000000E+05*basicFunctionsList[1][-2].evaluate(qVector[0]))
# print(1.200000000000E+05*(1-np.exp(-1.5*(qVector[0]-1.5521)))**2)
# print(1.000000000000E+05*basicFunctionsList[1][-3].evaluate(qVector[0]))
# print(1.000000000000E+05*(1-np.exp(-1.5*(qVector[0]-1.5521)))**4)
# print(1.600000000000E+05*basicFunctionsList[1][-4].evaluate(qVector[0])*basicFunctionsList[2][-4].evaluate(qVector[1])*basicFunctionsList[3][-2].evaluate(qVector[2]))
# print(1.600000000000E+05*np.exp(-0.500*(qVector[0]-1.5521)**2-0.500*(qVector[1]-1.5521)**2)*np.cos(0.5*qVector[2])**2)
# print(())
# print(basicFunctionsList[1].keys())
# print(basicFunctionsList[1][-1].evaluate(testQ))
# print(basicFunctionsList[2][-1].evaluate(testQ))
# print(basicFunctionsList[3][-1].evaluate(testQ))

# grid = np.linspace(1.2, 2.8, 1000)
# grid = np.linspace(85*toRadians, 180*toRadians, 1000)
# potential = np.zeros(1000)
# for i in range(len(grid)):
#     qVector[2] = grid[i]
#     potential[i] = operator.evaluatePointOfComponent("potential", basicFunctionsList, qVector)

# plt.plot(grid, potential)
# plt.savefig("CS2-Potential-alpha.png", dpi=300)