const ArrayOfVecs{T<:Real} = Array{V,1} where {V <: Vector{T}}
const ArrayofMats{T<:Real} = Array{M,1} where {M <: Matrix{T}}
const VecOfArrayOfVecs{T<:Real} = Vector{Vector{Vector{T}}}
const MatOrVecs{T <: Real} = Union{Matrix{T},ArrayOfVecs{T}}
const MatOrVec{T <: Real} = Union{Matrix{T},Vector{T}}
const Vec{T <: Real} = AbstractVector{T}
const Mat{T <: Real} = AbstractMatrix{T}

struct observer{T}
    posUnitVecRel :: Vector{T}
    distance :: T
    bodyFrameUnitVec :: Vector{T}
    positionVector :: Vector{T}
    rotation :: Matrix{T} # a rotation from the global coordinate frame into a frame where the posUnitVecRel is along the Z axis
end

struct pointSource{T}
    posUnitVecRel :: Vector{T}
    irradiance :: T
    bodyFrameUnitVec :: Vector{T}
    positionVector :: Vector{T}
    rotation :: Matrix{T} # a rotation from the global coordinate frame into a frame where the posUnitVecRel is along the Z axis
end

mutable struct clippedPolygon
    vertexIndices :: Vector{Int64}
    tempVertexIndices :: Vector{Int64}
    labels :: Vector{Int64}
    inOrOutBound :: Vector{Int64}
    intersectionNumber :: Vector{Int64}
    vertCheck :: Array{Bool,1}
    vertMap :: Vector{Int64}
    vertexNo :: Int64
    tempVertexNo :: Int64
    wrapPoints :: Vector{Int64}
    tempWrapPoints :: Vector{Int64}
    edge :: Vector{Vector{Float64}}
    labelstemp :: Vector{Int64}
    usedPoints :: Vector{Bool}
end

mutable struct surface{T}
    vertexIndices :: Vector{Int64}
    u_n :: Vector{T}
    u_u :: Vector{T}
    u_v :: Vector{T}
    Area :: T
    areaClipped :: T 
    Rdiff :: T
    Rspec :: T
    nu :: T
    nv :: T
    isVisible :: Bool 
    vertexNo :: Int64
    wrapPoints :: Vector{Int64}
    clippedGeometry :: Union{clippedPolygon, Nothing}
end

mutable struct objectProjectionData
    vertices :: Vector{Vector{Float64}}
    originalVertexNo :: Int
    relevantVertexNo :: Int
    isProjected :: Vector{Bool}
    projectedBounds :: Vector{Matrix{Float64}}
    intersection :: Vector{Vector{Float64}}
    polygons :: Vector{Vector{Int64}}
    labels :: Vector{Vector{Int64}}
    wrapPoints :: Vector{Vector{Int64}}
    wrapPointsOld :: Vector{Vector{Int64}}
    inOrOutBound :: Vector{Vector{Int64}}
    edge :: Vector{Vector{Float64}}
    usedPoints :: Vector{Vector{Bool}}
end

struct object{data <: Union{objectProjectionData, Nothing}, T, S <: surface{T}}
    vertices :: Vector{Vector{Float64}}
    surfaceNo :: Int64
    surfaces :: Vector{S}
    J :: Matrix{T}
    bodyFrame :: MatOrVecs{T}
    projectionData :: data
end

struct scenario{data <: Union{objectProjectionData, Nothing}, T, Obs <: observer{T}, obj <: object{data, T}}
    observerNo :: Int64
    observers :: Vector{Obs} 
    sun :: pointSource{T}
    objPosition :: Vector{T}
    object :: obj
    # only used to temporarily store the half angle vector between observers and the sun during calculations
    halfAngleVector :: Vector{T}
    unitScaling :: T
end

# legacy code
struct targetObject
    facetNo :: Int64
    Areas :: MatOrVec
    nvecs :: MatOrVecs
    vvecs :: MatOrVecs
    uvecs :: MatOrVecs
    nu :: MatOrVec
    nv :: MatOrVec
    Rdiff :: Union{MatOrVec, MatOrVecs}
    Rspec :: Union{MatOrVec, MatOrVecs}
    J :: Union{Matrix, Float64}
    targetObject() = new()
    targetObject(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J) = new(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J)
end

struct targetObjectFull
    facetNo :: Int64
    vertices :: Matrix
    vertList
    Areas :: MatOrVec
    nvecs :: MatOrVecs
    vvecs :: MatOrVecs
    uvecs :: MatOrVecs
    nu :: MatOrVec
    nv :: MatOrVec
    Rdiff :: Union{MatOrVec, MatOrVecs}
    Rspec :: Union{MatOrVec, MatOrVecs}
    J :: Matrix
    bodyFrame :: MatOrVecs
    targetObjectFull() = new()
    targetObjectFull(facetNo, vertices, vertList, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame) = new(facetNo, vertices, vertList, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
end

struct spaceScenario
    obsNo :: Int64
    C :: Number
    d :: MatOrVec
    sunVec :: Vector
    obsVecs :: MatOrVecs
    spaceScenario() = new()
    spaceScenario(obsNo, C, d, sunVec, obsVecs) = new(obsNo, C, d, sunVec, obsVecs)
end


