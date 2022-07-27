module lightCurveModeling
using LinearAlgebra
using attitudeFunctions
# using Infiltrator

include("types.jl")
include("utilities.jl")
include("lightCurveModel.jl")

export targetObject, targetObjectFull, spaceScenario, Fobs, dFobs, _Fobs, _Fobs_preAlloc, _Fobs_Analysis, simpleSatellite, simpleScenario, simpleScenarioGenerator, customScenarioGenerator, objectGenerator2D, scenarioGenerator2D, customSatellite, customScenario, attLMFIM, _mapp

end # module
#
