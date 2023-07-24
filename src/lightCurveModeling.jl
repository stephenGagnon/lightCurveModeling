module lightCurveModeling
using LinearAlgebra
using attitudeFunctions
using ForwardDiff
using Infiltrator
using Plots

include("types.jl")
include("utilities/utilities.jl")
include("utilities/selfShadowingUtils.jl")
include("utilities/copyUtilities.jl")
include("utilities/iterationUtilities.jl")
include("legacyCode.jl")
include("lightCurveModel.jl")
include("lightCurveGradient.jl")
include("objectGenerators.jl")
include("selfShadowing.jl")
include("plotting.jl")

export targetObject, targetObjectFull, spaceScenario, Fobs, dFobs, _Fobs, _Fobs_preAlloc, _Fobs_Analysis, simpleSatellite, simpleScenario, simpleScenario2D, simpleScenarioGenerator, customScenarioGenerator, objectGenerator2D, scenarioGenerator2D, customSatellite, customScenario, attLMFIM, _mapp, dFobs_2D, flatPlate3D, flatPlate3D_missing_section, forwardDiffWrapper

end # module
#
