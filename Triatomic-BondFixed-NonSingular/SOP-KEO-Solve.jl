using LinearAlgebra
using NumericalIntegration
using ForwardDiff

toRadians::Float64 = pi/180.0
wavenumberConversion::Float64 = 33.7152537836138732499014581824370 # conversion factor to cm-1
# Case of HOCl
masses::Vector{Float64} = [15.99059462, 1.00782503223, 34.965012622636]
qVector::Vector{Float64} = [0.964001600451, 1.69198245766, 102.693412447]
referenceGeometry::Vector{Float64} = [0.964001600451, 1.69198245766, 102.693412447]
morseParameters::Vector{Float64} = [2.33, 1.9]

numberOfIntegrationPoints::Vector{Int64} = [2000, 2000, 2000]
integrationBorders::Matrix{Float64} = [
    0.65 3.4;
    1.2 3.5;
    4.0 170.0]
# integrationBorders::Matrix{Float64} = [
#     -8 8;
#     1.2 3.5;
#     4.0 170.0]


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

potentialCheckpointFile::String = "potential.chk"
potentialInput::String = ""
open(potentialCheckpointFile, "r") do potentialCheckpointFileContent::IOStream
    global potentialInput = read(potentialCheckpointFileContent, String)
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

numberOfFunctionsPerMode::Vector{Int64} = zeros(numberOfModes)

push!(basicFunctionsList, 1 ./ masses)
headerIndex::Int64 = 1
for mode in 1:numberOfModes
    modeFunctionList::Vector{Function} = [] # New list of functions for mode
    push!(modeFunctionList, x -> 1) # Function 0 is always 1!
    numberOfFunctionsPerMode[mode] = parse(Int64, split(basicFunctionsInputLines[headerIndex])[end])
    for functionIndex in 1:numberOfFunctionsPerMode[mode]
        push!(modeFunctionList, basicFunctionTemplate(basicFunctionsInputLines[headerIndex + functionIndex]))
    end
    push!(basicFunctionsList, modeFunctionList)
    global headerIndex += numberOfFunctionsPerMode[mode] + 1
end


potentialBasicFunctions::Vector{Vector{Function}} = [
    [x -> 1 - exp(-morseParameters[1]*(x - referenceGeometry[1]))],
    [x -> 1 - exp(-morseParameters[2]*(x - referenceGeometry[2]))],
    [x -> cosd(x) - cosd(referenceGeometry[3])]
]

potentialInputLines::Vector{String} = split(String(split(potentialInput, "\nEnd of potential")[1]), "\n")
potentialInputLines = potentialInputLines[2:end-2]
numberOfPotentialTerms::Int64 = length(potentialInputLines)
potentialInputLinesSplit::Vector{Vector{String}} = split.(potentialInputLines)

potentialTermsCoefficients::Vector{Float64} = zeros(numberOfPotentialTerms)
potentialTermsIndices::Matrix{Int64} = zeros(numberOfPotentialTerms, numberOfModes)
for i in 1:numberOfPotentialTerms
    potentialTermsCoefficients[i] = parse(Float64, potentialInputLinesSplit[i][3])
    potentialTermsIndices[i, :] = parse.(Int64, potentialInputLinesSplit[i][4:end])
end

expansionOrderMaximum::Vector{Int64} = zeros(numberOfModes)
for i in 1:numberOfModes
    expansionOrderMaximum[i] = maximum(potentialTermsIndices[:, i])
end

# Divide kinetic input blocks
keywords::Vector{String} = ["Gvib", "Grot", "Gcor", "pseudo"]
padding::String = "\n987654321    0    0    0    0.00000000E+00   0"
for i in 1:numberOfModes
    global padding *= "    0"
end

kineticInputSplit::Vector{String} = split(kineticInput, padding)[1:4]
kineticVibrationalInput::Vector{String} = split(kineticInputSplit[1], "\n")[2:end]
kineticRotationalInput::Vector{String} = split(kineticInputSplit[2], "\n")[3:end]
kineticCoriolisInput::Vector{String} = split(kineticInputSplit[3], "\n")[3:end]
kineticPseudoPotentialInput::Vector{String} = split(kineticInputSplit[4], "\n")[3:end]

numberOfPseudoPotentialTerms::Int64 = length(kineticPseudoPotentialInput)
pseudoPotentialCoefficients::Vector{Float64} = zeros(numberOfPseudoPotentialTerms)
pseudoPotentialMapping::Matrix{Int64} = zeros(numberOfPseudoPotentialTerms, numberOfModes + 1)
pseudoPotentialInputLinesSplit::Vector{Vector{String}} = split.(kineticPseudoPotentialInput)
for i in 1:numberOfPseudoPotentialTerms
    pseudoPotentialCoefficients[i] = parse(Float64, pseudoPotentialInputLinesSplit[i][5])
    pseudoPotentialMapping[i, :] = parse.(Int64, pseudoPotentialInputLinesSplit[i][6:end])
end

numberOfVibrationalGMatrixTerms::Int64 = length(kineticVibrationalInput)
vibrationalCoefficients::Vector{Float64} = zeros(numberOfVibrationalGMatrixTerms)
vibrationalIndices::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}(undef, numberOfVibrationalGMatrixTerms)
vibrationalMapping::Matrix{Int64} = zeros(numberOfVibrationalGMatrixTerms, numberOfModes + 1)
vibrationalGMatrixLinesSplit::Vector{Vector{String}} = split.(kineticVibrationalInput)
for i in 1:numberOfVibrationalGMatrixTerms
    vibrationalCoefficients[i] = parse(Float64, vibrationalGMatrixLinesSplit[i][5])
    vibrationalMapping[i, :] = parse.(Int64, vibrationalGMatrixLinesSplit[i][6:end])
    vibrationalIndices[i] = Tuple(parse.(Int64, vibrationalGMatrixLinesSplit[i][1:2]))
end


numberOfRotationalGMatrixTerms::Int64 = length(kineticRotationalInput)
rotationalCoefficients::Vector{Float64} = zeros(numberOfRotationalGMatrixTerms)
rotationalIndices::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}(undef, numberOfRotationalGMatrixTerms)
rotationalMapping::Matrix{Int64} = zeros(numberOfRotationalGMatrixTerms, numberOfModes + 1)
rotationalGMatrixLinesSplit::Vector{Vector{String}} = split.(kineticRotationalInput)
for i in 1:numberOfRotationalGMatrixTerms
    rotationalCoefficients[i] = parse(Float64, rotationalGMatrixLinesSplit[i][5])
    rotationalMapping[i, :] = parse.(Int64, rotationalGMatrixLinesSplit[i][6:end])
    rotationalIndices[i] = Tuple(parse.(Int64, rotationalGMatrixLinesSplit[i][1:2]))
end

numberOfCoriolisGMatrixTerms::Int64 = length(kineticCoriolisInput)
coriolisCoefficients::Vector{Float64} = zeros(numberOfCoriolisGMatrixTerms)
coriolisIndices::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}(undef, numberOfCoriolisGMatrixTerms)
coriolisMapping::Matrix{Int64} = zeros(numberOfCoriolisGMatrixTerms, numberOfModes + 1)
coriolisGMatrixLinesSplit::Vector{Vector{String}} = split.(kineticCoriolisInput)
for i in 1:numberOfCoriolisGMatrixTerms
    coriolisCoefficients[i] = parse(Float64, coriolisGMatrixLinesSplit[i][5])
    coriolisMapping[i, :] = parse.(Int64, coriolisGMatrixLinesSplit[i][6:end])
    coriolisIndices[i] = Tuple(parse.(Int64, coriolisGMatrixLinesSplit[i][1:2]))
end

pseudoPotentialCoefficients .*= wavenumberConversion
vibrationalCoefficients .*= wavenumberConversion
rotationalCoefficients .*= wavenumberConversion
coriolisCoefficients .*= wavenumberConversion

# Index of this array is i = n + 1 where n denotes the polynomial degree
maxPolyad::Int64 = 48
maxPolyadForEachMode::Vector{Int64} = [48, 48, 48]

function initializeHarmonicOscillatorBasis(basisSize::Int64, x::Vector{Float64})

    ψ = zeros(length(x), basisSize+1)

    ψ[:,1] .= π^(-0.25) .* exp.(-x.^2 ./2)

    if basisSize == 0
        return ψ
    end

    ψ[:,2] .= sqrt(2) .* x .* ψ[:,1]

    for n in 1:basisSize-1
        ψ[:,n+2] .= sqrt(2/(n+1)) .* x .* ψ[:,n+1] .-
                    sqrt(n/(n+1)) .* ψ[:,n]
    end

    return ψ
end

# println(vibrationalIndices)
# println(findall(x -> x == (2, 2), vibrationalIndices))
integrationGrid::Vector{Vector{Float64}} = []
for i in 1:numberOfModes
    push!(integrationGrid, collect(LinRange(integrationBorders[i, 1]-referenceGeometry[i], integrationBorders[i, 2]-referenceGeometry[i], numberOfIntegrationPoints[i])))
end

oscillatorBasis::Matrix{Float64} = initializeHarmonicOscillatorBasis(maxPolyad, integrationGrid[1])

println(integrate(integrationGrid[1], oscillatorBasis[:, 1].*oscillatorBasis[:, 1], SimpsonEven()))
println(integrate(integrationGrid[1], oscillatorBasis[:, 1].*oscillatorBasis[:, 2], SimpsonEven()))
println(integrate(integrationGrid[1], oscillatorBasis[:, 1].*oscillatorBasis[:, 3], SimpsonEven()))
println(integrationGrid[1])
# println(integrationGrid[1])
potentialBasicFunctionIntegrals::Vector{Vector{Matrix{Float64}}} = []
for i in 1:1 #:numberOfModes # Number of modes
    # newPotentialIntegralsList::Vector{Matrix{Float64}} = []
    # for j in 0:expansionOrderMaximum[i]
    #     newIntegralMatrix::Matrix{Float64} = zeros(maxPolyadForEachMode[i] + 1, maxPolyadForEachMode[i] + 1)
    #     for k in 1:maxPolyadForEachMode[1]+1
    #         for l in 1:k
    #             newIntegralMatrix[k, l] = integrate(integrationGrid[i], oscillatorFunctions[k].(integrationGrid[i]).*potentialBasicFunctions[i][1].(integrationGrid[i]).^j.*oscillatorFunctions[l].(integrationGrid[i]), SimpsonEven())
    #         end
    #     end
    #     push!(newPotentialIntegralsList, Hermitian(newIntegralMatrix))
    # end
    # push!(potentialBasicFunctionIntegrals, newPotentialIntegralsList)
end

# println(potentialBasicFunctionIntegrals)

# println(ForwardDiff.derivative(tan, 1))
# println(sec(1)^2)
println(integrate(testRange, oscillatorFunctions[1].(testRange).*oscillatorFunctions[1].(testRange), SimpsonEven()))
