const ArrayOfVecs{T<:Real} = Array{V,1} where {V <: Vector{T}}
const ArrayofMats{T<:Real} = Array{M,1} where {M <: Matrix{T}}
const VecOfArrayOfVecs{T<:Real} = Vector{Vector{Vector{T}}}
const MatOrVecs{T <: Real} = Union{Matrix{T},ArrayOfVecs{T}}
const MatOrVec{T <: Real} = Union{Matrix{T},Vector{T}}
const Vec{T <: Real} = AbstractVector{T}
const Mat{T <: Real} = AbstractMatrix{T}

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
    J :: Matrix
    bodyFrame :: MatOrVecs
    targetObject() = new()
    targetObject(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame) = new(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
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
    targetObjectFull(facetNo, vertices, vertList, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame) = new(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
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
