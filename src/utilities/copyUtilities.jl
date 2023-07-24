import Base.insert!
"""
    overload the insert! function to handle a case where the relevant items in the collection do not take up the entire collecction. In this case, we can shift data up one element to make space for a new item without allocating new memory
    the upper bound input is the highest index of relevant data. If the upper bound is greater than or equal to the length of the collection, the method just calls the base intsert! function
"""
function insert!(collection, index, item :: AbstractVector, upperBound)
    
    if index == upperBound + 1 && length(collection) > upperBound
        if isassigned(collection, index)
            collection[index][:] = item
        else
            collection[index] = copy(item)
        end
    elseif length(collection) > upperBound && index <= upperBound # if there are more allocated elements than relevant elements (as defined by the upperBound which is the highest index of relevant elements), use the allocated memory instead of inserting (which allocated new memory)

        #shift elements of collection beyond the index of insertion up by 1 to make space for the new element
        for i = upperBound+1:-1:index+1 #iterate through elements of collection in descending order
            collection[i] = collection[i-1]
        end
        collection[index][:] = item

    else
        temp = similar(collection)
        append!(collection,temp)
        insert!(collection, index, item, upperBound)
        # Base.insert!(collection, index, item)
    end

    return upperBound += 1
end

function insert!(collection, index, item :: Number, upperBound)
    
    if length(collection) > upperBound # if there are more allocated elements than relevant elements (as defined by the upperBound which is the highest index of relevant elements), use the allocated memory instead of inserting (which allocated new memory)

        #shift elements of collection beyond the index of insertion up by 1 to make space for the new element
        for i = upperBound+1:-1:index+1 #iterate through elements of collection in descending order
            collection[i] = collection[i-1]
        end
        collection[index] = item

    else
        temp = similar(collection)
        append!(collection,temp)
        insert!(collection, index, item, upperBound)
        # Base.insert!(collection, index, item)
    end

    return upperBound += 1
end

import Base.replace!
"""
    overload the replace! function
    this function sets the value of a collection at a specified index to the given value. It will append as necessary if the selected index is outised the range of the array
"""
function replace!(collection :: Vector{T}, index :: Int, value :: T) where T
    if index <= lastindex(collection)
        collection[index] = value
    # elseif index == lastindex(collection) + 1
    #     append!(collection,value)
    else
        temp = similar(collection)
        append!(collection, temp)
        collection[index] = value
    end
end

function replace!(collection :: Vector{T}, rangeVal, value :: Vector{T}) where T
    if rangeVal[end] <= lastindex(collection)
        for i in rangeVal
            collection[i] = value[i]
        end
    else
        # temp = Array{typeof(value),1}(undef,rangeVal[end] - lastindex(collection))
        temp = similar(collection)
        append!(collection, temp)
        # collection[rangeVal] = value
        replace!(collection :: Vector{T}, rangeVal :: AbstractVector, value :: Vector{T})
    end
end

function replace!(collection :: Vector{T}, rangeVal :: AbstractVector, value :: T) where T
    if rangeVal[end] <= lastindex(collection)
        for i in rangeVal
            collection[i] = value
        end
    else
        # temp = Array{typeof(value),1}(undef,rangeVal[end] - lastindex(collection))
        temp = similar(collection)
        append!(collection, temp)
        for i in rangeVal
            collection[i] = value
        end
    end
end

import Base.deleteat!
"""
    overload the deleteat! function
    this function shifts elements of the collection over without reducing the size of the array, creating "unused" elements at the end. It also updates an integer storing the last relevant element of the array.
    This is so that we avoid deleting some memory only to reallocate it later.
"""
function deleteat!(collection, index :: Int64, maxindex)
    for i = index:maxindex-1
        collection[i] = collection[i+1]
    end
    return maxindex - 1
end

function deleteat!(collection, range :: AbstractVector, maxindex)
    l = length(range)
    if l>0
        for i in range
            if i + l <= maxindex
                collection[i] = collection[i+l]
            else
                collection[i] = -2
            end
        end
        return maxindex - l
    else
        return maxindex
    end
end

# function deleteat!(collection, range :: AbstractVector, maxindex, fillValue)
#     l = length(range)
#     if l>0
#         for i in range[1]:maxindex
#             if i + l <= maxindex
#                 collection[i] = collection[i+l]
#             else
#                 collection[i] = fillValue
#             end
#         end
#         return maxindex - l
#     else
#         return maxindex
#     end
# end

# assumes range is continuous
function deleteat!(collection, range :: AbstractVector, maxindex, fillValue)
    l = length(range)
    if l>0
        chunkSize = 1
        for i = l:-1:1
            if i == 1
                removeChunk!(collection, range[i], chunkSize, maxindex, fillValue)
                maxindex += -chunkSize
            elseif range[i-1] == range[i]-1
                chunkSize += 1
            else
                removeChunk!(collection, range[i], chunkSize, maxindex, fillValue)
                maxindex += -chunkSize
                chunkSize = 1
            end
        end
    end
    return maxindex
end

function removeChunk!(collection, startIndex, chunkSize, maxindex ,fillValue)
    for j = startIndex:maxindex
        if j + chunkSize <= maxindex
            collection[j] = collection[j + chunkSize]
        else
            collection[j] = fillValue
        end
    end
end

function deleteWrapPoints!(collection, range :: AbstractVector, maxindex)
    l = length(range)
    if l>0
        chunkSize = 1
        for i = l:-1:1
            if i == 1
                for j = range[i]:maxindex
                    if j + chunkSize <= maxindex
                        if collection[j + chunkSize] > 0
                            collection[j] = 1
                        else
                            collection[j] = 0
                        end
                    else
                        collection[j] = 0
                    end
                end
                maxindex += -chunkSize
            elseif range[i-1] == range[i]-1
                chunkSize += 1
            else
                for j = range[i]:maxindex
                    if j + chunkSize <= maxindex
                        if collection[j + chunkSize] > 0
                            collection[j] = 1
                        else
                            collection[j] = 0
                        end
                    else
                        collection[j] = 0
                    end
                end
                maxindex += -chunkSize
                chunkSize = 1
            end
        end
    end

    prevPoint = 1
    for i = 1:maxindex
        if collection[i] == 1
            collection[i] = prevPoint
            prevPoint = i + 1
        end
    end
    
    return maxindex
end

function copyVertices!(vertListTemp, vertList, vertNo)
    if length(vertListTemp) >= vertNo
        vertListTemp[1:vertNo] = copy(vertList[1:vertNo])
    else
        temp = similar(vertListTemp)
        append!(vertListTemp, temp)
        copyVertices!(vertListTemp, vertList, vertNo)
        # vertListTemp[:] = vertList[1:length(vertListTemp)]
        # append!(vertListTemp, vertList[length(vertListTemp)+1:vertNo])
    end
end

function updateVertMap(vertMap, index)
    for i in eachindex(vertMap)
        if vertMap[i] >= index
            vertMap[i] += 1
        end
    end
end

function resetVertMap!(vertMap)
    for i in eachindex(vertMap)
        vertMap[i] = i
    end
end 

function updateClippingData!(poly :: clippedPolygon, maxIndex)
    # update the number of relevant vertices
    poly.vertexNo = maxIndex

    # these are all vectors of ints, so we can do them together in a loop without incuring type instability
    symbols = [:tempVertexIndices, :labels, :intersectionNumber, :vertMap]

    for i in eachindex(symbols)
        var = getproperty(poly,symbols[i]) :: Vector{Int64}
        if length(var) < maxIndex
            append!(var, Array{Int64,1}(undef, maxIndex))
        end
    end

    # do vertCheck seperately becaue it is a vector of booleans
    vertCheck = getproperty(poly, :vertCheck) :: Vector{Bool}
    if length(vertCheck) < maxIndex
        append!(vertCheck, Array{Bool,1}(undef, maxIndex))
    end


    # if length(poly.labels) < maxIndex
    #     temp = Array{Int64,1}(undef,length(poly.labels))
    #     append!(poly.labels, temp)
    # end

end

function setIntersections(value, check = true)
    if check
        return value
    else
        return 0
    end
end

"""
    add a point to a polygon
    adds the actual vertex coordinates to the vertex list
    adds the vertex index to the list of vertex indices
"""
function addPoint!(point, subjectPolygon, index, wrapPoints, oldInd, data, planeVertIndices)

    vertices = data.vertices
    # if the point to be added is on the clipping polygon, a new point in the subject polygon plane needs to be generated
    if point[1] == 1
        if length(vertices) > data.relevantVertexNo 
            if !isassigned(vertices, data.relevantVertexNo + 1)
                vertices[data.relevantVertexNo + 1] = similar(data.vertices[oldInd])
            end
            vert = vertices[data.relevantVertexNo + 1]
            for i = 1:2
                vert[i] = vertices[oldInd][i]
            end
            vert[3] = fullPointPlane(vert, vertices[planeVertIndices]) # compute z coordinate
        else
            vert = similar(vertices[oldInd])
            for i = 1:2
                vert[i] = vertices[oldInd][i]
            end
            vert[3] = fullPointPlane(vert, vertices[planeVertIndices]) # compute z coordinate
            append!(vertices, [vert])
        end
        data.relevantVertexNo += 1
        vertexIndex = data.relevantVertexNo
    else
        vertexIndex = oldInd
    end
    
    replace!(subjectPolygon, index, vertexIndex)
    replace!(wrapPoints[2], index, 0)

end

function resetValues!(collection, size, fillValue)
    if length(collection)<size
        append!(collection, similar(collection))
        resetValues!(collection, size, fillValue)
    else
        collection .= fillValue
    end
end