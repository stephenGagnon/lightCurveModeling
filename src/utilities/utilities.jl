function attLMFIM(att, sat, scen, R, parameterization = quaternion)

        if any(parameterization .== [MRP GRP quaternion])
             if typeof(att) <: Union{MRP,GRP}
                 rotFunc = ((A,v) -> any2A(A).A*v)
             elseif typeof(att) == quaternion
                 rotFunc = qRotate
             elseif length(att) == 3
                 rotFunc = ((A,v) -> p2A(A,a,f)*v)
             elseif length(att) == 4
                 rotFunc = qRotate
             end
        elseif parameterization == DCM
            if typeof(att) == Matrix
                rotFunc = ((A,v) -> A*v)
            elseif typeof(att) == DCM
                rotFunc = ((A,v) -> A.A*v)
            end
        elseif parameterization == att2D
            if typeof(att) <: Real
                rotFunc = ((tht,v) -> [cos(tht) sin(tht);-sin(tht) cos(tht)]*v)
            elseif typeof(att) == att2D
                rotFunc = ((tht,v) -> [cos(tht.tht) sin(tht.tht);-sin(tht.tht) cos(tht.tht)]*v)
            end
        else
            error("Please provide a valid attitude. Attitudes must be represented as a single 3x1 or 4x1 float array, a 3x3 float array, or any of the custom attitude types defined in the attitudeFunctions package. For a 2D object, attitude can be representated as a real number or as a att2D struct.")
        end

    dhdx = ForwardDiff.jacobian(A -> _Fobs( A, sat.nvecs, sat.uvecs, sat.vvecs, sat.Areas, sat.nu, sat.nv, sat.Rdiff, sat.Rspec, scen.sunVec, scen.obsVecs, scen.d, scen.C, rotFunc), att)

    # dhdx = dFobs(att, sat.nvecs, sat.uvecs, sat.vvecs, sat.Areas, sat.nu, sat.nv, sat.Rdiff, sat.Rspec, scen.sunVec, scen.obsVecs, scen.d, scen.C, dDotFunc, rotFunc, quaternion)

    return dhdx'*inv(R)*dhdx
end

function forwardDiffWrapper(func, dim)
    result = DiffResults.GradientResult(Array{Float64,1}(undef,dim))
    func_out = (x, grad :: Vector) ->
    begin
        if length(grad) > 0
            # @infiltrate
            result = ForwardDiff.gradient!(result, func, x)
            fval = DiffResults.value(result)
            grad[:] = DiffResults.gradient(result)
        else
            fval = func(x)
        end
        return fval
    end
    return func_out
end

function forwardDiffWrapper(func, dim, storeCost)
    result = DiffResults.GradientResult(Array{Float64,1}(undef, dim))
    func_out = (x, grad::Vector) ->
        begin
            if length(grad) > 0
                # @infiltrate
                result = ForwardDiff.gradient!(result, func, x)
                fval = DiffResults.value(result)
                push!(storeCost,fval)
                grad[:] = DiffResults.gradient(result)
            else
                fval = func(x)
            end
            return fval
        end
    return func_out
end

function dot3(a :: Vec, b :: Vec)
    
    if length(a) != 3 || length(b) != 3
        throw(DimensionMismatch("vectors have incorect length"))
    end
    dot3_unsafe(a, b)
end

function dot3_unsafe(x, y)
    s = 0.
    for k in 1:3
        @inbounds s += x[k]*y[k]
    end
    s
end

function dot2(a :: Vec, b :: Vec)
    
    if length(a) != 2 || length(b) != 2
        throw(DimensionMismatch("vectors have incorect length"))
    end
    dot2_unsafe(a, b)
end

function dot2_unsafe(x, y)
    s = 0.
    for k in 1:2
        @inbounds s += x[k]*y[k]
    end
    s
end

function dot2_unsafe(x)
    dot2_unsafe(x, x)
end

import LinearAlgebra.dot
function dot(a)
    return dot(a,a)
end

"""
    Creates a rotation matrix which transforms into a coordinate system where the input vector is the X axis
"""
function Xbasis(vec :: Vector{N}) where {N <: Number}
    rotMat = Array{N,2}(undef,3,3)
    rotMat[1,:] = vec
    if all(vec .== [1.0;0;0])
        rotMat = [1.0 0 0;0 1 0;0 0 1]
    elseif !all(vec .== [0;0;-1.0])
        rotMat[2,:] = cross(rotMat[1,:], [0;0;-1.0])
        rotMat[2,:] = rotMat[2,:]./norm(rotMat[2,:])
        rotMat[3,:] = cross(rotMat[1,:], rotMat[2,:])
        rotMat[3,:] = rotMat[3,:]./norm(rotMat[3,:])
    else
        rotMat[2,:] = cross(rotMat[1,:], [0;1.0;0])
        rotMat[2,:] = rotMat[2,:]./norm(rotMat[2,:])
        rotMat[3,:] = cross(rotMat[1,:], rotMat[2,:])
        rotMat[3,:] = rotMat[3,:]./norm(rotMat[3,:])
    end

    return rotMat
end

"""
    Creates a rotation matrix which transforms into a coordinate system where the input vector is the X axis
"""
function Zbasis(vec :: Vector{N}) where {N <: Number}
    rotMat = Array{N,2}(undef,3,3)
    rotMat[3,:] = vec
    if all(vec .== [0;0;1.0])
        rotMat = [1.0 0 0;0 1 0;0 0 1]
    elseif !all(vec .== [0;0;-1.0])
        rotMat[2,:] = cross(rotMat[3,:], [0;0;-1.0])
        rotMat[2,:] = rotMat[2,:]./norm(rotMat[2,:])
        rotMat[1,:] = cross(rotMat[2,:], rotMat[3,:])
        rotMat[1,:] = rotMat[1,:]./norm(rotMat[1,:])
    else
        rotMat[2,:] = cross(rotMat[3,:], [0;1.0;0])
        rotMat[2,:] = rotMat[2,:]./norm(rotMat[2,:])
        rotMat[1,:] = cross(rotMat[2,:], rotMat[3,:])
        rotMat[1,:] = rotMat[1,:]./norm(rotMat[1,:])
    end

    return rotMat
end

function XYextrema!(point, minMax :: Matrix)

    # assumes 2D
    # minMax array is structured as:
    #[maximum x,    minimum x;
    # maximum y,    minumum y]

    # max x
    if isnan(minMax[1,1])
        minMax[1,1] = point[1]
    elseif minMax[1,1] < point[1]
        minMax[1,1] = point[1]
    end

    # min x
    if isnan(minMax[1,2])
        minMax[1,2] = point[1]
    elseif minMax[1,2] > point[1]
        minMax[1,2] = point[1]
    end

    # max y
    if isnan(minMax[2,1])
        minMax[2,1] = point[2]
    elseif minMax[2,1] < point[2]
        minMax[2,1] = point[2]
    end

    # min y
    if isnan(minMax[2,2])
        minMax[2,2] = point[2]
    elseif minMax[2,2] > point[2]
        minMax[2,2] = point[2]
    end

end

function boundsCheck(M1 :: Matrix, M2 :: Matrix)
    # checks to see if two boxes with upper/lower bounds given by M1 and M2 overlap
    # assumes 2D
    # M1 and M2 are structured as:
    #[maximum x,    minimum x;
    # maximum y,    minumum y;
    # maximum z,    minimum z]

    if M1[1,1] > M2[1,2] && M2[1,1] > M1[1,2] && M1[2,1] > M2[2,2] && M2[2,1] > M1[2,2]
        return true
    else
        return false
    end
end

function orderIntersections(jj, vertList, vertListFull, vertMap, vertices, intersection)

    # check if there is already an intersection on this edge and find the appropriate index to insert the new intersection

    nextInd = nextIndex(jj, vertList)

    if vertList[jj + 1] == -1
        basePoint = vertices[vertList[jj]]
        d = dot2_unsafe(intersection - basePoint)

        for i = (vertMap[jj] + 1):(vertMap[jj + 1]-1)
            d2 = dot2_unsafe(vertices[vertListFull[i]] - basePoint)
            if d < d2
                return i
            end
        end

        d2 = dot2_unsafe(vertices[vertList[nextInd]] - basePoint)
        if d < d2
            return vertMap[jj + 1]
        else
            error("intersection does not lie on line segement between vertMap[jj] and the next point")
        end

    elseif vertMap[nextInd] - vertMap[jj] > 1
        # if jj < vertNo - 1 && vertMap[jj + 1] - vertMap[jj] > 1
            basePoint = vertices[vertList[jj]]
            d = dot2_unsafe(intersection - basePoint)
    
            for i = (vertMap[jj] + 1):vertMap[nextInd]
                d2 = dot2_unsafe(vertices[vertListFull[i]] - basePoint)
                if d < d2
                    return i
                end
            end
    
            error("intersection does not lie on line segement between vertMap[jj] and the next point")
    else
        return vertMap[jj] + 1
    end

    # if function has not returned yet, no intersection was found
    error("insertion point not found")
end

"""
    fullPoint!(point, adjPoints)

    linearly interpolates the first element of the point from the adjactent points  
"""
function fullPoint!(point, adjPoints)
    point[3] = normDiff2(point, adjPoints[1])/normDiff2(adjPoints[2], adjPoints[1]) * (adjPoints[2][3] - adjPoints[1][3]) + adjPoints[1][3]
    # if !isapprox(sum(cross(point - adjPoints[1], point - adjPoints[2]) - zeros(3)), 0)
    #     return false
    # else 
    #     return true
    # end
end

"""
    takes norm of the difference between two 2D vectors 
    accepts vector of any length > 2 but only considers the first two elements
"""
function normDiff2(point1,point2)
    return sqrt((point1[1]- point2[1])^2 + (point1[2]- point2[2])^2)
end

function fullPointPlane(point, Pl)
    # simple code that is readable
    # v1 = Plane[1] - Plane[2]
    # v2 = Plane[1] - Plane[3]
    # n = cross(v1,v2)
    # return (1/n[3])*(dot(n,Plane[1]) - n[1]*point[1] - n[2]*point[2]) 

    # expressions are expanded to avoid allocations in linear algebra functions
    n1 = (Pl[1][2] - Pl[2][2])*(Pl[1][3] - Pl[3][3]) - (Pl[1][3] - Pl[2][3])*(Pl[1][2] - Pl[3][2])
    n2 = (Pl[1][3] - Pl[2][3])*(Pl[1][1] - Pl[3][1]) - (Pl[1][1] - Pl[2][1])*(Pl[1][3] - Pl[3][3])
    n3 = (Pl[1][1] - Pl[2][1])*(Pl[1][2] - Pl[3][2]) - (Pl[1][2] - Pl[2][2])*(Pl[1][1] - Pl[3][1])

    return n1/n3 * (Pl[1][1] - point[1]) + n2/n3 * (Pl[1][2] - point[2]) + Pl[1][3]
end

function other12(v)
    if v == 1
        return 2
    elseif v == 2
        return 1
    end
end

function ShoelaceAlgorithm(points)
    sum = points[end][1]*points[1][2] - points[1][1]*points[end][2]
    for i = 1:length(points) - 1
        sum += points[i][1]*points[i+1][2] - points[i+1][1]*points[i][2]
    end
    return sum
end

function ShoelaceAlgorithm(indices, points, wrapPoints, lastIndex)

    area = points[indices[1]][1]*points[indices[2]][2] - points[indices[2]][1]*points[indices[1]][2]

    for i = 2:lastIndex-1
        ind = indices[i]
        if indices[i + 1] != -1
            indp = indices[i + 1]
        else
            indp = indices[wrapPoints[i + 1]]
        end

        if ind == -1
            continue
        else
            area += points[ind][1]*points[indp][2] - points[indp][1]*points[ind][2]
        end

    end
    return area
end

function ShoelaceAlgorithm(indices, points, lastIndex)

    area = points[indices[lastIndex-1]][1]*points[indices[1]][2] - points[indices[1]][1]*points[indices[lastIndex-1]][2]

    for i = 1:lastIndex-2
        ind = indices[i]
        if indices[i + 1] != -1
            indp = indices[i + 1]
        else
            for j = i:-1:1
                if indices[j] == -1
                    indp = indices[j + 1]
                elseif indices[j] == 1
                    indp = indices[1]
                end
            end
            # indp = indices[wrapPoints[i + 1]]
        end

        if ind == -1
            continue
        else
            area += points[ind][1]*points[indp][2] - points[indp][1]*points[ind][2]
        end

    end
    return area
end

function elementsExist(collection, elements)
    val = true
    for e in elements
        if !any(collection .== e)
            val = false
            break
        end
    end
    return val
end

function allInterior(labels, range)
    for i in range
        if labels[i] != -1 || labels[i] != -2
            return false
        end
    end
    return true
end

function allExterior(labels, range)
    for i in range
        if labels[i] != 0 || labels[i] != -2
            return false
        end
    end
    return true
end

function allInterior(labels, range, inOrOutBound)
    # if any point is exterior or outbound, then return false
    for i in range
        if labels[i] == 0 || inOrOutBound[i] == -1
            return false
        end
    end
    return true
end

function allExterior(labels, range, inOrOutBound)
    # if any point is interior or inbound, then return false
    for i in range
        if labels[i] == -1 || inOrOutBound[i] == 1
            return false
        end
    end
    return true
end

function clippingEdgeCaseHandling(surface1, surface2, clipNo)

    # set up to traverse the polygons (clipping polygon with be traversed counter clockwise, subject polygon will be traversed clockwise)
    if clipNo == 1
        clippingPoly = surface1.clippedGeometry
        subjectPoly =  surface2.clippedGeometry
    elseif clipNo == 2
        clippingPoly = surface2.clippedGeometry
        subjectPoly =  surface1.clippedGeometry
    end

    polygon_clip = clippingPoly.tempVertexIndices
    polygon_sub = subjectPoly.tempVertexIndices

    labels_clip = clippingPoly.labels
    labels_sub = subjectPoly.labels

    subjectPolyInds = subjectPoly.vertexIndices
    # wrapPoints = (clippingPoly.tempWrapPoints, subjectPoly.tempWrapPoints)
    wrapPoints_clip = clippingPoly.tempWrapPoints
    wrapPoints_sub = subjectPoly.tempWrapPoints
    # inOrOutBound = (clippingPoly.inOrOutBound, subjectPoly.inOrOutBound)
    inOrOutBound_clip = clippingPoly.inOrOutBound
    inOrOutBound_sub = subjectPoly.inOrOutBound

    vertNo_clip = clippingPoly.tempVertexNo
    vertNo_sub = subjectPoly.tempVertexNo
    wrapVal = 1

    # for each set of connected points, check if all points are interior or exterior. If all are interior, remove the points, if all are exterior, all points should be added to the clipped polygon
    # can definitely be optimized (broadcasting conditional, slicing, and findall are all allocating)
    wrapIndexes1 = findall(wrapPoints_clip[1:vertNo_clip] .> 0)
    wrapIndexes2 = findall(wrapPoints_sub[1:vertNo_sub] .> 0)

    index = 0
    # could be preallocated
    deletionIndexes = [Int64[],Int64[]]

    for i in eachindex(wrapIndexes1), j in eachindex(wrapIndexes2)

        if i == 1
            r1 = 1:wrapIndexes1[i]
        else
            r1 = wrapIndexes1[i-1]+1:wrapIndexes1[i]
        end

        if j == 1
            r2 = 1:wrapIndexes2[j]
        else
            r2 = wrapIndexes2[j-1]+1:wrapIndexes2[j]
        end

        # if all of the points in the subject region are interior points, and all of the points in the clipping region are exterior points, then the region is shadowed and should be removed
        # if all of the points of the clipping region are interior points, and all of the points of the subject region are exterior points, then the clipping creates a hole in the subject polygon. Label the clipping points as already used and add them to the subject polygon in reverse order

        if allInterior(labels_sub, r2, inOrOutBound_sub) && allExterior(labels_clip, r1, inOrOutBound_clip)
            # subject region is entirely inside clipping region
            # get indexes of intersections in clipping ploy that need to be removed and add them to the list
            append!(deletionIndexes[1], filter(x -> x > 0, labels_sub[r2]))
            # add the range of indexes in the subject polygon that need to be removed
            append!(deletionIndexes[2],r2)
        elseif allInterior(labels_clip, r1, inOrOutBound_clip) && allExterior(labels_sub, r2, inOrOutBound_sub)
            # clipping region is entirely inside subject polygon
            # add clipping polygon vertices in counter-clockwise order
            for i = r1[end-1]:-1:r1[1]
                index += 1
                replace!(subjectPolyInds, index, polygon_clip[r1[i]])
                labels_clip[r1[i]] = -3
                replace!(subjectPoly.wrapPoints, index, 0)
            end
            index += 1
            replace!(subjectPolyInds, index, -1)
            replace!(subjectPoly.wrapPoints, index, wrapVal)
            wrapVal = index + 1   
        end

    end

    # remove all the relevant intersections 
    deleteat!(polygon_sub, deletionIndexes[2], vertNo_sub, -2)
    deleteat!(inOrOutBound_sub, deletionIndexes[2], vertNo_sub, -3)
    deleteWrapPoints!(wrapPoints_sub, deletionIndexes[2], vertNo_sub)
    vertNo_sub = deleteat!(labels_sub, deletionIndexes[2], vertNo_sub, -2)
    subjectPoly.tempVertexNo = vertNo_sub

    deleteat!(polygon_clip, deletionIndexes[1], vertNo_clip, -2)
    deleteat!(inOrOutBound_clip, deletionIndexes[1], vertNo_clip, -3)
    deleteWrapPoints!(wrapPoints_clip, deletionIndexes[1], vertNo_clip)
    vertNo_clip = deleteat!(labels_clip, deletionIndexes[1], vertNo_clip, -2)
    clippingPoly.tempVertexNo = vertNo_clip

    # if all regions of the subject polygon are shadowed, mark the facet as non-visible
    if vertNo_sub == 0
        if clipNo == 1
            surface2.isVisible = false
        elseif clipNo == 2
            surface1.isVisible = false
        end
        return true
    end

    if index != 0 
        index += 1 # add one for the final spacing value which only exists if points were added
    end

    # update the number of relevant vertices
    if clipNo == 1
        surface2.clippedGeometry.vertexNo = index
    elseif clipNo == 2
        surface1.clippedGeometry.vertexNo = index
    end

    return false
end

"""
    pointLineBoundsCheck(point, xmax, xmin, ymax, ymin, tol)
    Checks if point is within the bounds of a line segement
"""
function pointLineBoundsCheck(point, xmax, xmin, ymax, ymin, tol = 1e-10)
    return point[1] - xmin > tol && xmax - point[1] > tol && point[2] - ymin > tol && ymax - point[2] > tol 
end

"""
    Assumes lines are both edges of a counter clockwise oriented polygons
    Determines if line 2 is inbound or outbound wrt the polygon containing line 1
    returns 1 if line 2 is inbound
    returns -1 if line 2 is outbound
    returns 0 if lines are parallel
"""
function inOrOutBoundCheck(line1,line2, tol = 1e-10)

    # # shared base points are indeterminate
    # if abs(line1[1][1] - line2[1][1]) < 1e-12 && abs(line1[1][2] == line2[1][2]) < 1e-12
    #     return 0
    # end

    v1 = line1[2] - line1[1]
    v2 = line2[2] - line2[1]
    temp = cross2D(v1,v2)
    if abs(temp) < tol
        return 0
    else
        return Int64(sign(temp))
    end
end

"""
    deals with the case where the intersection is at a shared vertex
    returns 1 if line is inbound
    returns -1 if line is outbound
    returns 0 if lines are parallel
"""
function inOrOutBoundCheckSharedVertex(line1, line2, pp1, pp2, tol = 1e-10)

    v1a = line1[2]-line1[1]
    v1b = pp1-line1[1]
    v2a = line2[2]-line2[1]
    v2b = pp2-line2[1]

    c1 = cross2D(v1a,v1b)
    c2 = cross2D(v2a,v2b)

    c1a2a = cross2D(v1a,v2a)
    c1a2b = cross2D(v1a,v2b)
    c2a1b = cross2D(v2a,v1b)

    # c1b2b = cross(v1b,v2b)

    if c1a2a > 0
        if c1a2b > 0 && c2 < 0
            inBound1 = 1
        elseif abs(c1a2b) < tol
            inBound1 = 0
        else
            inBound1 = -1
        end

        if c1 > 0 && c2a1b < 0
            inBound2 = -1
        elseif abs(c2a1b) < tol
            inBound2 = 0
        else
            inBound2 = 1
        end
    elseif c1a2a < 0 
        if c2 > 0 && c1a2b < 0
            inBound1 = -1
        elseif abs(c1a2b) < tol
            inBound1 = 0
        else
            inBound1 = 1
        end

        if c2a1b > 0 && c1 < 0
            inBound2 = 1
        elseif abs(c1) < tol
            inBound2 = 0
        else
            inBound2 = -1
        end
    end

    return inBound1, inBound2

end

"""
    deals with the case where the base point of the first line is on the second line
    returns 1 if line is inbound
    returns -1 if line is outbound
    returns 0 if lines are parallel
"""
function inOrOutBoundCheckVertInt(line1, line2, pp1, tol = 1e-10)

    v1a = line1[2] - line1[1]
    v1b = pp1 - line1[1]
    v2 = line2[2] - line2[1]

    c1 = cross2D(v1a,v1b)

    c1a2 = cross2D(v1a,v2)
    c21b = cross2D(v2,v1b)

    inBound1 = Int64(-sign(c1a2))
    
    if c1a2 > 0
        if c1 > 0 && c21b < 0
            inBound2 = -1
        elseif abs(c21b) < tol
            inBound2 = 0
        else
            inBound2 = 1
        end
    elseif c1a2 < 0
        if c21b > 0 && c1 < 0
            inBound2 = 1
        elseif abs(c1) < tol
            inBound2 = 0
        else
            inBound2 = -1
        end
    end

    return inBound1, inBound2
end

function cross2D(v1,v2)
    return v1[1]*v2[2] - v1[2]*v2[1]
end

"""
    checks if the first two elements of a vector are apprxoimately equal without any slicing or copying
"""
function vectorEquality2D(v1,v2, tol = 1e-10)
    return abs(v1[1] - v2[1]) < tol && abs(v1[2] - v2[2]) < tol
end

"""
    checks if point lies in the direction of the second point of the line
    assumes the point is on the line
"""
function lineDirectionCheck(line, point)
    temp = 0.0
    for i = 1:2
        temp += (line[2][i] - line[1][i])*(point[i] - line[1][i])
    end
    return sign(temp) == 1
end

## unused
    # function previousPoint(vertList, index)
    #     if vertList[index - 1]
    # end



    # function unsafe_Array_Insert_2D!(array, insertArray, rowRange, colRange)
    #     # assumes insertArray has the same size as [length(rowRange), length(colRange)]
    #     for (x, r) in enumerate(rowRange), (y,c) in enumerate(colRange)
    #         array[r,c] = insertArray[x,y]
    #     end
    # end

    # function appendVertex(vertices, newVertex, vertNoFull)
    #     # add intersection to list of vertices
    #     vertNoFull += 1
    #     if size(vertices,2) > vertNoFull
    #         # vertices[:, vertNoFull] = copy(newVertex)
    #         unsafe_Array_Insert_2D!(vertices, newVertex, 1:3, vertNoFull)
    #     else
    #         # append!(vertices, intersection)
    #         verticestemp = Array{eltype(vertices),2}(undef,2*size(vertices,2))
    #         # verticestemp[size(vertices)] = vertices
    #         unsafe_Array_Insert_2D!(verticestemp, vertices, 1:3, 1:vertNoFull-1)
    #         # verticestemp[:,vertNoFull] = intersection
    #         unsafe_Array_Insert_2D!(verticestemp, newVertex, 1:3, vertNoFull)
    #         vertices = verticestemp
    #     end
    #     return vertices
    # end

    # function view!(collection, outputArray, range1, range2)
    #     for (i, r) in enumerate(range1)
    #         outputArray[i][:] = view(collection[r], range2)
    #     end
    # end

