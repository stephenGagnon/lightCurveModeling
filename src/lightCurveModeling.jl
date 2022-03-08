module lightCurveModeling
using LinearAlgebra
using attitudeFunctions

include("types.jl")
include("utilities.jl")
include("lightCurveModel.jl")

export targetObject, targetObjectFull, spaceScenario, Fobs, dFobs, _Fobs_Analysis, simpleSatellite, simpleScenario, customScenarioGenerator, customSatellite, customScenario

end # module
