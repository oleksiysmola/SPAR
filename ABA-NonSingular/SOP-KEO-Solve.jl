using LinearAlgebra

toRadians::Float64 = pi/180.0
wavenumberConversion::Float64 = 33.7152537836138732499014581824370 # conversion factor to cm-1
# Case of H2S
masses::Vector{Float64} = [31.9720707, 1.00782503223]
qVector::Vector{Float64} = [1.336, 1.336, 92.110]
numberOfModes::Int64 = length(qVector)

# Read inputs
basicFunctionsInputFile::String = "combination.chk"
basicFunctionsInput::String = ""
open(basicFunctionsInputFile, "r") do basicFunctionsFileContent::IOStream
    global basicFunctionsInput = read(basicFunctionsFileContent, String)
end

kineticCheckpointFile::String = "kinetic.chk"
kineticInput::String = ""
open(kineticCheckpointFile, "r") do kineticCheckpointFileContent::IOStream
    global kineticInput = read(kineticCheckpointFileContent, String)
end

basicFunctionsInputLines::Vector{String} = split(String(split(basicFunctionsInput, "\nEND")[1]), "\n")
basicFunctionsInputLines = basicFunctionsInputLines[2:end]

elementaryFunctionMap::Dict{String, Function} = Dict([
    ("I", x -> x),
    ("COS", cosd),
    ("SIN", sind),
    ("TAN", tand),
    ("SEC", secd),
    ("CSC", cscd),
    ("COT", cotd)
]
)

function basicFunctionTemplate(basicFunctionLine::String)::Function
    basicFunctionLineSplit::Vector{String} = split(basicFunctionLine, r"\s+")
    if basicFunctionLineSplit[2] == "2"
        n1 = parse(Float64, basicFunctionLineSplit[3])
        f1 = elementaryFunctionMap[uppercase(basicFunctionLineSplit[4])]
        a1 = parse(Float64, basicFunctionLineSplit[5])
        k1 = parse(Float64, basicFunctionLineSplit[6])
        n2 = parse(Float64, basicFunctionLineSplit[7])
        f2 = elementaryFunctionMap[uppercase(basicFunctionLineSplit[8])]
        a2 = parse(Float64, basicFunctionLineSplit[9])
        k2 = parse(Float64, basicFunctionLineSplit[10])
        return q -> f1(a1*q^k1)^n1*f2(a2*q^k2)^n2
    else
        n1 = parse(Float64, basicFunctionLineSplit[3])
        f1 = elementaryFunctionMap[uppercase(basicFunctionLineSplit[4])]
        a1 = parse(Float64, basicFunctionLineSplit[5])
        k1 = parse(Float64, basicFunctionLineSplit[6])
        return q -> f1(a1*q^k1)^n1
    end
end


basicFunctionsList::Vector{Vector{Any}} = []

push!(basicFunctionsList, 1 ./ masses)
headerIndex::Int64 = 1
println(basicFunctionsInputLines)
println(split(basicFunctionsInputLines[1])[end])
for mode in 1:numberOfModes
    modeFunctionList::Vector{Function} = [] # New list of functions for mode
    push!(modeFunctionList, x -> 1) # Function 0 is always 1!
    numberOfFunctionsForMode = parse(Int64, split(basicFunctionsInputLines[headerIndex])[end])
    println(headerIndex," " , numberOfFunctionsForMode)
    for functionIndex in 1:numberOfFunctionsForMode
        push!(modeFunctionList, basicFunctionTemplate(basicFunctionsInputLines[headerIndex + functionIndex]))
    end
    push!(basicFunctionsList, modeFunctionList)
    global headerIndex += numberOfFunctionsForMode + 1
end

# println(basicFunctionsList[1])
# println(basicFunctionTemplate("12 2 2 Sin 0.5 1 2 Tan 0.5 1")(1))
# println(basicFunctionTemplate("12 1 1 Sin 1 1")(90))
println(basicFunctionsList[4][12](90))

hermitePolynomials::Vector{Function} = [x -> 1, x -> 2x]
    

# # Divide kinetic input blocks
# keywords = ["Gvib\n", "Grot\n", "Gcor\n", "pseudo\n"]
# padding = "\n987654321    0    0    0    0.00000000E+00   0"
# for i in range(numberOfModes):
#     padding += "    0"
# padding += "\n"
# kineticInputSplit = kineticInput.split(padding)[:4]
# kineticVibrationalInput = kineticInputSplit[0].split(keywords[0])[1]
# kineticRotationalInput = kineticInputSplit[1].split(keywords[1])[1]
# kineticCoriolisInput = kineticInputSplit[2].split(keywords[2])[1]
# kineticPseudoPotentialInput = kineticInputSplit[3].split(keywords[3])[1]

# GMatrixVibrational = np.zeros((numberOfModes, numberOfModes))
# GMatrixRotational = np.zeros((3, 3))
# GMatrixCoriolis = np.zeros((numberOfModes, 3))
# pseudoPotential = 0.0

# kineticPseudoPotentialInputLines = kineticPseudoPotentialInput.split("\n")
# for line in kineticPseudoPotentialInputLines:
#     kineticLineSplit = line.split()
#     newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
#     for k in range(1, numberOfModes+1):
#         newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
#     pseudoPotential += newTerm


# kineticVibrationalInputLines = kineticVibrationalInput.split("\n")
# for line in kineticVibrationalInputLines:
#     kineticLineSplit = line.split()
#     newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
#     for k in range(1, numberOfModes+1):
#         newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
#     GMatrixVibrational[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

# kineticRotationalInputLines = kineticRotationalInput.split("\n")
# for line in kineticRotationalInputLines:
#     kineticLineSplit = line.split()
#     newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
#     for k in range(1, numberOfModes+1):
#         newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
#     GMatrixRotational[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

# kineticCoriolisInputLines = kineticCoriolisInput.split("\n")
# for line in kineticCoriolisInputLines:
#     kineticLineSplit = line.split()
#     newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
#     for k in range(1, numberOfModes+1):
#         newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
#     GMatrixCoriolis[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

# pseudoPotential *= wavenumberConversion
# GMatrixVibrational *= wavenumberConversion
# GMatrixRotational *= wavenumberConversion
# GMatrixCoriolis *= wavenumberConversion
# print(pseudoPotential)
# print(GMatrixVibrational)
# print(GMatrixRotational)
# print(GMatrixCoriolis)