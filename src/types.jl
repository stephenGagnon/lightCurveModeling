const Vec{T<:Number} = AbstractArray{T,1}
const Mat{T<:Number} = AbstractArray{T,2}
const ArrayOfVecs{T<:Number} = Array{V,1} where V <: Vec
const ArrayofMats{T<:Number} = Array{M,1} where M <: Mat
const VecOfArrayOfVecs{T<:Number} = Vector{Vector{Vector{T}}}
const MatOrVecs = Union{Mat,ArrayOfVecs}
const MatOrVec = Union{Mat,Vec}
const anyAttitude = Union{Mat,Vec,DCM,MRP,GRP,quaternion}
const Num = N where N <: Number

struct targetObject
    facetNo :: Int64
    Areas :: MatOrVec
    nvecs :: MatOrVecs
    vvecs :: MatOrVecs
    uvecs :: MatOrVecs
    nu :: MatOrVec
    nv :: MatOrVec
    Rdiff :: MatOrVec
    Rspec :: MatOrVec
    J :: Mat
    bodyFrame :: MatOrVecs
    targetObject() = new()
    targetObject(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame) = new(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
end

struct targetObjectFull
    facetNo :: Int64
    vertices :: Mat
    vertList
    Areas :: MatOrVec
    nvecs :: MatOrVecs
    vvecs :: MatOrVecs
    uvecs :: MatOrVecs
    nu :: MatOrVec
    nv :: MatOrVec
    Rdiff :: MatOrVec
    Rspec :: MatOrVec
    J :: Mat
    bodyFrame :: MatOrVecs
    targetObjectFull() = new()
    targetObjectFull(facetNo, vertices, vertList, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame) = new(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
end

struct spaceScenario
    obsNo :: Int64
    C :: N where {N <: Number}
    d :: MatOrVec
    sunVec :: Vec
    obsVecs :: MatOrVecs
    spaceScenario() = new()
    spaceScenario(obsNo, C, d, sunVec, obsVecs) = new(obsNo, C, d, sunVec, obsVecs)
end
