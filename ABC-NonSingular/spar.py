import numpy as np

wavenumberConversion = 33.7152537836138732499014581824370 # conversion factor to cm-1

class basicFunction:
    basicFunctionType: str
    numberOfFunctions: int
    primitiveFunctionList: list

    def __init__(self, basicFunctionInputLine: str, isSeries: bool = False, order: int = 0) -> None:
        '''isSeries: Is basic function part of a series expansion
           order: Order for a repeated function type
        '''
        basicFunctionLineSplit: list = basicFunctionInputLine.split()
        self.numberOfFunctions = int(basicFunctionLineSplit[1])
        primitiveFunctionList: list = []
        startIndex = 2
        p = float(basicFunctionLineSplit[startIndex])
        functionString = basicFunctionLineSplit[startIndex + 1].lower()
        b = float(basicFunctionLineSplit[startIndex + 2])
        c = float(basicFunctionLineSplit[startIndex + 3])
        match functionString:
            case "sin":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order: np.sin(order*q)]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.sin(b*q**c)**p]
            case "cos":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order: np.cos(order*q)]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.cos(b*q**c)**p]
            case "power":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order, c=c: (q-c)**order]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(q-c))**p]
            case "morse":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order, b=b, c=c, p=p: (1-np.exp(-b*(q-c)))**order]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1-np.exp(-b*(q-c)))**p]
            case "cos(q)-cos(q0)":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order, b=b, c=c: (b*(np.cos(q)-np.cos(c)))**order]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(np.cos(q)-np.cos(c)))**p]
            case "sin(q)-sin(q0)":
                if isSeries:
                    primitiveFunctionList += [lambda q, order=order, b=b, c=c: (b*(np.sin(q)-np.sin(c)))**order]
                else:
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(np.sin(q)-np.sin(c)))**p]
            case "exp":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.exp(b*(q-c)**p)]
            case "tan":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.tan(b*q**c)**p]
            case "sec":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.cos(b*q**c))**p]
            case "csc":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.sin(b*q**c))**p]
            case "cot":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.tan(b*q**c))**p]
            case "tanh":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.tanh(b*(q-c)**p)]
            case "1-tanh":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: 1-np.tanh(b*(q-c)**p)]
            case "i":
                primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*q**c)**p]

        startIndex += 4
        for i in range(self.numberOfFunctions - 1):
            p = float(basicFunctionLineSplit[startIndex])
            functionString = basicFunctionLineSplit[startIndex + 1].lower()
            b = float(basicFunctionLineSplit[startIndex + 2])
            c = float(basicFunctionLineSplit[startIndex + 3])
            match functionString:
                case "sin":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.sin(b*q**c)**p]
                case "cos":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.cos(b*q**c)**p]
                case "power":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(q-c))**p]
                case "morse":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1-np.exp(-b*(q-c)))**p]
                case "cos(q)-cos(q0)":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(np.cos(q)-np.cos(c)))**p]
                case "sin(q)-sin(q0)":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*(np.sin(q)-np.sin(c)))**p]
                case "exp":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.exp(b*(q-c)**p)]
                case "tan":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.tan(b*q**c)**p]
                case "sec":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.cos(b*q**c))**p]
                case "csc":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.sin(b*q**c))**p]
                case "cot":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (1/np.tan(b*q**c))**p]
                case "tanh":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: np.tanh(b*(q-c)**p)]
                case "1-tanh":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: 1-np.tanh(b*(q-c)**p)]
                case "i":
                    primitiveFunctionList += [lambda q, b=b, c=c, p=p: (b*q**c)**p]

        self.primitiveFunctionList = primitiveFunctionList
    
    def evaluate(self, q: float) -> float:
        functionOutput: float = 1.0
        for i in range(self.numberOfFunctions):
            functionOutput *= self.primitiveFunctionList[i](q)
        return functionOutput
    
# Distinguishes between the fourier and taylor series types
seriesTypes: dict = {
    "cos": "fourier",
    "sin": "fourier",
    "power": "taylor",
    "morse": "taylor",
    "cos(q)-cos(q0)": "taylor",
    "sin(q)-sin(q0)": "taylor"
}
        
def readBasicFunctions(basicFunctionInput: str) -> dict:
    basicFunctionsList: dict = {}
    numberOfModes: int = basicFunctionInput.lower().count("mode")
    basicFunctionInputLines: list = basicFunctionInput.split("\nEND")[0].split("\n")
    basicFunctionMainHeader: str = basicFunctionInputLines[0]
    basicFunctionInputLines = basicFunctionInputLines[1:]
    modeHeaderIndex: int = 0
    if "series" in basicFunctionMainHeader:
        for i in range(numberOfModes):
            modeFunctionList = {} # New list of functions for mode
            modeFunctionList[0] = basicFunction("0 1 0 I 1 1") # Function 0 is always 1!
            numberOfLinesForMode: int = int(basicFunctionInputLines[modeHeaderIndex].split()[-1])
            lineBeingRead: int = 1
            functionIndexCounter: int = 1
            while lineBeingRead <= numberOfLinesForMode:
                basicFunctionInputLine: str = basicFunctionInputLines[modeHeaderIndex + lineBeingRead]
                basicFunctionInputLineSplit: list = basicFunctionInputLine.split()
                functionIndex: int = int(basicFunctionInputLineSplit[0])
                if int(basicFunctionInputLineSplit[0]) < 0:
                    modeFunctionList[functionIndex] = basicFunction(basicFunctionInputLine)
                else:
                    functionLabel: str = basicFunctionInputLineSplit[3].lower()
                    seriesType: str = seriesTypes[functionLabel]
                    orderIndex: int = 2
                    if seriesType == "taylor":
                        orderIndex = 2
                    else:
                        orderIndex = 4
                    maxOrder: int = int(basicFunctionInputLineSplit[orderIndex])
                    for j in range(1, maxOrder+1):
                        modeFunctionList[functionIndexCounter] = basicFunction(basicFunctionInputLine, True, j)
                        functionIndexCounter += 1
                lineBeingRead += 1
            basicFunctionsList[i + 1] = modeFunctionList
            modeHeaderIndex += lineBeingRead
    else:
        for i in range(numberOfModes):
            modeFunctionList = {} # New list of functions for mode
            modeFunctionList[0] = basicFunction("0 1 0 I 1 1") # Function 0 is always 1!
            numberOfLinesForMode: int = int(basicFunctionInputLines[modeHeaderIndex].split()[-1])
            lineBeingRead: int = 1
            while lineBeingRead <= numberOfLinesForMode:
                modeFunctionList[lineBeingRead] = basicFunction(basicFunctionInputLines[modeHeaderIndex + lineBeingRead])
                lineBeingRead += 1
            basicFunctionsList[i] = modeFunctionList
            modeHeaderIndex += lineBeingRead
    return basicFunctionsList

# 0 for potentials, 1 for dms, 2 for kinetic
operatorRankMapping = {
    "kinetic": 2,
    "dipole": 1,
    "potential": 0
}

# Labels for kinetic energy operator components
kineticComponents = ["gvib", "grot", "gcor", "pseudo"]

class operatorMap:
    operatorType: str
    operatorRank: int
    componentIndices: dict = {}
    componentCoefficients: dict = {}
    functionIndices: dict = {}
    numberOfModes: int
    containsMass: bool = False

    def __init__(self, operatorMappingInput: str):
        operatorMappingInputLines: list = operatorMappingInput.split("\nEND")[0].split("\n")
        operatorMappingHeader: str = operatorMappingInputLines[0].lower()
        self.operatorType = operatorMappingHeader.split()[1]
        self.containsMass = "mass" in operatorMappingHeader
        self.operatorRank = operatorRankMapping[self.operatorType]
        self.numberOfModes = len(operatorMappingInputLines[2].split()) - self.operatorRank - 3 - self.containsMass
        if self.operatorType == "kinetic":
            pass
        else:
            operatorComponentInputLines: list = operatorMappingInputLines[2 :]
            numberOfComponentTerms: int = len(operatorComponentInputLines)
            self.componentCoefficients[self.operatorType] = np.zeros(numberOfComponentTerms)
            if self.operatorRank > 0:
                self.componentIndices[self.operatorType] = np.zeros((numberOfComponentTerms, self.operatorRank))
            self.functionIndices[self.operatorType] = np.zeros((numberOfComponentTerms, self.numberOfModes + self.containsMass))
            for i in range(numberOfComponentTerms):
                operatorComponentLineSplit: list = operatorComponentInputLines[i].split()
                self.componentCoefficients[self.operatorType][i] = float(operatorComponentLineSplit[self.operatorRank + 2])
                for j in range(self.operatorRank):
                    self.componentIndices[self.operatorType][i, j] = int(operatorComponentLineSplit[j])
                for j in range(self.numberOfModes + self.containsMass):
                    self.functionIndices[self.operatorType][i, j] = int(operatorComponentLineSplit[j + self.operatorRank + 3])
    
    def evaluatePointOfComponent(self, component: str, basicFunctionsList: dict, q):
        coefficients = self.componentCoefficients[component]
        functionIndices = self.functionIndices[component]
        numberOfTerms: int = len(coefficients)
        match self.operatorRank:
            case 0:
                operatorValue = 0.0
                for i in range(numberOfTerms):
                    newTerm = coefficients[i]
                    if self.containsMass:
                        newTerm *= basicFunctionsList[0][functionIndices[i, 0]]
                        for j in range(1, self.containsMass + self.numberOfModes):
                            newTerm *= basicFunctionsList[j][functionIndices[i, j]].evaluate(q[j])
                    else:
                        for j in range(self.numberOfModes):
                            newTerm *= basicFunctionsList[j + 1][functionIndices[i, j]].evaluate(q[j])
                    operatorValue += newTerm
                return operatorValue
            case 1:
                operatorValue = np.zeros(3)
            case 2:
                match component:
                    case "pseudo":
                        operatorValue = 0.0
                    case "gvib":
                        operatorValue = np.zeros((self.numberOfModes, self.numberOfModes))
                    case "grot":
                        operatorValue = np.zeros((3, 3))
                    case "gcor":
                        operatorValue = np.zeros((self.numberOfModes, 3))

# def readCheckpointFile(checkpointFileInput: str, containsMass: bool = False) -> dict:
#     checkpointFileLineSplit: list =  checkpointFileInput.split("\n")
#     operatorType: str = checkpointFileLineSplit[0].split()[-1].lower()
#     rankOfComponent: int = componentRankMapping[operatorType]
#     pass