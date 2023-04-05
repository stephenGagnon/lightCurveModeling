module lightCurveModeling
using LinearAlgebra
using attitudeFunctions
using ForwardDiff
# using Infiltrator

include("types.jl")
include("utilities.jl")
include("lightCurveModel.jl")

export targetObject, targetObjectFull, spaceScenario, Fobs, dFobs, _Fobs, _Fobs_preAlloc, _Fobs_Analysis, simpleSatellite, simpleScenario, simpleScenario2D, simpleScenarioGenerator, customScenarioGenerator, objectGenerator2D, scenarioGenerator2D, customSatellite, customScenario, attLMFIM, _mapp, dFobs_2D, flatPlate3D, flatPlate3D_missing_section

end # module
#
