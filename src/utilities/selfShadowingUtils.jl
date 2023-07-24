"""
    modified implementation of the weiler Atherton clipping algorithm, gets the clipped subject polygon instead of the intersecting polygon
    also determines which polygon is "on top" and which should be clipped
"""
function polygonClipping!(data, surface1, surface2, diagnostics = false)

    exitFlag = 0

    vertices = data.vertices

    intersection = data.intersection 

    # check which polygon is on top
    # clipNo = checkPolyDistance(vertList, vertices, intersection, [surface1.vertexNo, surface2.vertexNo], intNo, vertCheck)
    clipNo = checkPolyDistance(vertices, intersection, surface1, surface2, data.edge)

    if clipNo == 0
        # no common point could be found (no interior vertices or intersections)
        # this means there is no intersection of the polygons
        return -2
    end

    # if diagnostics
    #     print("")
    # end

    # generate full vertex lists, including all intersections and overlaps
    # polygonIntersections!(vertNo, vertList, vertListFull, labels, intNo, vertCheck, vertMap, vertices, vertexCount, intersection, clippingPoly)
    polygonIntersections!(data, surface1, surface2, clipNo, diagnostics)

    # edge cases
    done = clippingEdgeCaseHandling(surface1, surface2, clipNo)

    if done
        return 1
    end

    # set up to traverse the polygons (clipping polygon with be traversed counter clockwise, subject polygon will be traversed clockwise)
    if clipNo == 1
        clippingPoly = surface1.clippedGeometry
        subjectPoly =  surface2.clippedGeometry
        subPolyOriginal = surface2.vertexIndices
    elseif clipNo == 2
        clippingPoly = surface2.clippedGeometry
        subjectPoly =  surface1.clippedGeometry
        subPolyOriginal = surface1.vertexIndices
    else
        error("invalid clipNo")
    end

    # polygons = (clippingPoly.tempVertexIndices, subjectPoly.tempVertexIndices)
    # labels = (clippingPoly.labels, subjectPoly.labels)
    # wrapPoints = (clippingPoly.wrapPoints, subjectPoly.wrapPoints)
    # inOrOutBound = (clippingPoly.inOrOutBound, subjectPoly.inOrOutBound)
    # wrapPointsOld = (clippingPoly.tempWrapPoints, subjectPoly.tempWrapPoints)

    polygons = data.polygons
    polygons[1] = clippingPoly.tempVertexIndices
    polygons[2] = subjectPoly.tempVertexIndices

    labels = data.labels
    labels[1] = clippingPoly.labels
    labels[2] = subjectPoly.labels

    wrapPoints = data.wrapPoints
    wrapPoints[1] = clippingPoly.wrapPoints
    wrapPoints[2] = subjectPoly.wrapPoints

    inOrOutBound = data.inOrOutBound
    inOrOutBound[1] = clippingPoly.inOrOutBound
    inOrOutBound[2] = subjectPoly.inOrOutBound

    wrapPointsOld = data.wrapPointsOld
    wrapPointsOld[1] = clippingPoly.tempWrapPoints
    wrapPointsOld[2] = subjectPoly.tempWrapPoints

    usedPoints = data.usedPoints
    usedPoints[1] = clippingPoly.usedPoints
    usedPoints[2] = subjectPoly.usedPoints

    resetValues!(usedPoints[1], clippingPoly.tempVertexNo, false)
    resetValues!(usedPoints[2], subjectPoly.tempVertexNo, false)

    subjectPolygon = subjectPoly.vertexIndices
    clippingVertNo = clippingPoly.tempVertexNo
    subjectVertNo = subjectPoly.tempVertexNo


    startingPoint = [0;0]
    # get a starting point
    found = findStartingPoint!(startingPoint, labels, inOrOutBound, wrapPointsOld, clippingVertNo, subjectVertNo, usedPoints)

    # if no valid starting point is found
    if !found
        # this should not happen
        error("unable to find starting point")
    end

    # clipping
    finished = false
    index = 1
    oldPoint = copy(startingPoint)
    currentPoint = copy(startingPoint)
    wrapVal = 1

    # add starting point to polygon
    addPoint!(currentPoint, subjectPolygon, index, wrapPoints, polygons[currentPoint[1]][currentPoint[2]], data, subjectPoly.vertexIndices[1:3])


    # print(labels[1][1:clippingPoly.tempVertexNo], "\n")
    # print(labels[2][1:subjectPoly.tempVertexNo], "\n")
    # print(startingPoint, "\n")

    while !finished

        oldPoint[:] = currentPoint

        # get the next point
        flag, output = nextPoint!(currentPoint, oldPoint, labels, wrapPointsOld, startingPoint, inOrOutBound, usedPoints)

        if flag == -1 || flag == -3
            # error("invalid point reached")
            return -1
        end
        
        # mark previous point as used (in both polygons if appropriate)
        # rework 
        if labels[oldPoint[1]][oldPoint[2]] > 0 && flag != 1
            usedPoints[other12(oldPoint[1])][labels[oldPoint[1]][oldPoint[2]]] = true
        end
        usedPoints[oldPoint[1]][oldPoint[2]] = true   

        # if we reach the original point so we have formed a closed shape
        if flag == -2
            # print(currentPoint, " ", -3, "\n")  
            # put a space in the vertex list indicating that this is the end of a close shape
            replace!(subjectPolygon, index + 1, -1)
            replace!(wrapPoints[2], index  + 1, wrapVal)

            index += 2
            wrapVal = index 

            # try to find another starting point
            found = findStartingPoint!(startingPoint, labels, inOrOutBound, wrapPointsOld, clippingVertNo, subjectVertNo, usedPoints)

            # if no starting point exists, end the procedure
            if !found
                finished = true
                continue
            end

            currentPoint[:] = startingPoint

        else# otherwise continue adding points
            # print(currentPoint, " ", labels[currentPoint[1]][currentPoint[2]] , "\n")   
            index += 1
        end

        # add current point to subject polygon
        addPoint!(currentPoint, subjectPolygon, index, wrapPoints, polygons[currentPoint[1]][currentPoint[2]], data, subPolyOriginal[1:3])

    end

    # update sizes of arrays to account for new clipped polygon
    updateClippingData!(subjectPoly, index - 1)
    return 0

end

"""
    Code to find all intersections of two polygons
    Modifies the vertListFull arrays in place to be clockwise ordered lists of points on each polygon including all intersections
    Modifies the labels arrays in place to be lists determining the status of each point on the full polygon
    -1 -> interior vertex
    0 -> exterior vertex
    n > 0 -> intersection
    where n is the index of corresponding point on the other polygon
"""
function polygonIntersections!(data, surface1 :: surface, surface2 :: surface, clippingPoly, diagnostics = false) 
    # polygonIntersections!(vertListFull, vertList, vertNoOriginal, labels, intNo, vertCheck, vertMap, vertices, vertexCount, intersection)

    intersection = data.intersection[1]
    vertices = data.vertices
    vertexCount = data.relevantVertexNo
    edge = data.edge

    if clippingPoly == 0
        return 0
    end

    vertNoOrig1, vertNoOrig2, vertNo1, vertNo2, vertList1, vertList2, vertListFull1, vertListFull2, labels1, labels2, labelstemp1, labelstemp2, inOrOutBound1, inOrOutBound2, intNo1, intNo2, vertCheck1, vertCheck2, vertMap1, vertMap2, edge1, edge2, wrapPoints1, wrapPoints2 = unpacking(surface1, surface2, clippingPoly)

    # if diagnostics
    #     if clippingPoly == 1
    #         s1 =surface1
    #         s2 = surface2
    #     else
    #         s1 =surface2
    #         s2 = surface1
    #     end
    #     plotShape2D(s1.vertexIndices,vertices,:green)
    #     plotShape2D!(s2.clippedGeometry.vertexIndices[1: s2.clippedGeometry.vertexNo],vertices,:red)
    #     # plotFacets(surface1.clippedGeometry[1:vertNo1], surface2.clippedGeometry, data.vertices)
    #     print("")
    # end
   
    # used for mapping intersections between points 
    # each intersection in a pair is labeled with mapNo and then it is incremented so that each pair has a unique number
    # at the end, these numbers are replaced with the indicies of the corresponding intersection in the opposite polygon
    mapNo = 1

    for ii = 1:vertNoOrig1 - 1
        # get the next edge from each polygon
        if ii == 1
            edge1[1] = vertices[vertList1[1]]
            edge1[2] = vertices[vertList1[2]]
            prevPoint1 = vertices[vertList1[previousIndex(1, vertList1, vertNoOrig1)]] 
        elseif vertList1[ii] != -1
            prevPoint1 = nextEdge!(edge1, vertices, vertList1, ii)
        else
            break # has to be break to exit the inner loop entirely and increment ii by 1
        end

        for jj = 1:vertNoOrig2 - 1
            if jj == 1
                edge2[1] = vertices[vertList2[1]]
                edge2[2] = vertices[vertList2[2]]
                prevPoint2 = vertices[vertList2[previousIndex(1, vertList2, vertNoOrig2)]]
            elseif vertList2[jj] != -1
                prevPoint2 = nextEdge!(edge2, vertices, vertList2, jj)
            else
                continue # continue increments the inner loop by 1
            end


            # if diagnostics && ii == 60 && (jj==8||jj==9)
            #     print("")
            # end
            
            # if diagnostics && ii == 4 && jj == 1
            #     print("")
            # end
            # get intersection for the current pair of edges
            intersectionExists, flag1, flag2, inbound1, inbound2 = lineIntersection!(edge1, edge2, intersection, vertCheck1[ii], vertCheck2[jj], prevPoint1, prevPoint2)

            # edge case handling
            # shared vertex
            if flag1 == -2 && flag2 == -2
                # mark the intersecting vertices as inbound or outbound
                replace!(inOrOutBound2, vertMap2[jj], inbound2)
                replace!(inOrOutBound1, vertMap1[ii], inbound1)

                # set label to current mapNo (increment when the same point on the other object is set)
                labels1[vertMap1[ii]] = mapNo
                labels2[vertMap2[jj]] = mapNo
                mapNo += 1
                
                vertCheck1[ii] = false
                vertCheck2[jj] = false
                continue
            elseif flag1 == -2 || flag2 == -2
                error("inconsistent return: if there is a shared vertex, both flags should be -2")
            end

            # edge one 
            if flag1 == 1
                # increment number of ray intersections associated with each point (used later to check interior/exterior status)
                intNo1[ii] += flag1
            elseif flag1 == -1
                # unknown number of intersections due to edges being co-linear
                # directly check if the point is interior or exterior and flag that it no longer needs to be checked
                labels1[vertMap1[ii]] = interiorPointCheck(edge, vertList1[ii], vertList2, vertNoOrig2, vertices)
                vertCheck1[ii] = false
            elseif flag1 == -3

                # vertex is on an edge from the other shape
                # check if there is already an intersection on this edge and find the appropriate index to insert the new intersection
                insertionIndex = orderIntersections(jj, vertList2, vertListFull2, vertMap2, vertices, vertices[vertList1[ii]])

                vertexCount = insert!(vertices, vertexCount + 1, copy(vertices[vertList1[ii]]), vertexCount) 

                # mark the intersecting vertices as inbound or outbound
                # temp = inOrOutBoundCheck(edge1,edge2)
                insert!(inOrOutBound2, insertionIndex, inbound2, vertNo2)
                replace!(inOrOutBound1, vertMap1[ii], inbound1)

                # label vertex in the first shape as being an edge vertex and insert a label for it into the second shape labels
                labels1[vertMap1[ii]] = mapNo
                insert!(labels2, insertionIndex, mapNo, vertNo2)
                insert!(wrapPoints2, insertionIndex, 0, vertNo2)
                mapNo += 1
                
                # # need to add the vertex as a separate point with a new set of 3D coordinates that place it on the other shape when de-projected
                # intersection[1][:] = vertices[vertList1[ii]]
                # fullPoint!(intersection[1], edge2)

                # insert vertex index into the other shapes vertex list at the appropriate location
                vertNo2 = insert!(vertListFull2, insertionIndex, vertexCount, vertNo2)
                updateVertMap(vertMap2, insertionIndex)
                vertCheck1[ii] = false

            elseif flag1 == 0
            else
                error("undefined flag")
            end

            # edge 2
            if flag2 == 1
                # increment number of ray intersections associated with each point (used later to check interior/exterior status)
                intNo2[jj] += flag2
            elseif flag2 == -1
                # unknown number of intersections due to edges being parallel
                # directly check if the point is interior or exterior 
                labels2[vertMap2[jj]] = interiorPointCheck(edge, vertList2[jj], vertList1, vertNoOrig1, vertices)
                vertCheck2[jj] = false
            elseif flag2 == -3
                # vertex is on edges from the other shape
                # check if there is already an intersection on this edge and find the appropriate index to insert the new intersection
                insertionIndex = orderIntersections(ii, vertList1, vertListFull1, vertMap1, vertices, vertices[vertList2[jj]])

                # # add the vertex as a separate point with a new set of 3D coordinates that place it on the other shape when de-projected
                # intersection[2][:] = vertices[vertList2[jj]]
                # fullPoint!(intersection[2], edge1)

                vertexCount = insert!(vertices, vertexCount + 1, copy(vertices[vertList2[jj]]), vertexCount) 

                # mark the intersecting vertices as inbound or outbound
                # temp = inOrOutBoundCheck(edge1,edge2)
                replace!(inOrOutBound2, vertMap2[jj], inbound2)
                insert!(inOrOutBound1, insertionIndex, inbound1, vertNo1)

                # label vertex as being an edge vertex
                labels2[vertMap2[jj]] = mapNo
                insert!(labels1, insertionIndex, mapNo, vertNo1)
                insert!(wrapPoints1, insertionIndex, 0, vertNo1)
                mapNo += 1
                
                # insert vertex into the other shape at the appropriate location
                vertNo1 = insert!(vertListFull1, insertionIndex, vertexCount, vertNo1)
                updateVertMap(vertMap1, insertionIndex)
                vertCheck2[jj] = false

            elseif flag2 == 0
            else
                error("??????")
            end

            if intersectionExists  
                # for each shape   
                # find the full point including the z coordinate for de-projection
                # insert the full point to the list of vertices
                # find where the point belongs in the list of vertices
                # insert the vertex label into the list of labels
                # insert the vertex index into the list of indices
                # update the vertexMap 

                fullPoint!(intersection, edge2)
  
                vertexCount = insert!(vertices, vertexCount + 1, intersection, vertexCount) 
                insertionIndex1 = orderIntersections(ii, vertList1, vertListFull1, vertMap1, vertices, intersection)
                insert!(labels1, insertionIndex1, mapNo, vertNo1)
                insert!(wrapPoints1, insertionIndex1, 0, vertNo1)
                vertNo1 = insert!(vertListFull1, insertionIndex1, vertexCount, vertNo1)
                updateVertMap(vertMap1, insertionIndex1)

                # fullPoint!(intersection[2], edge2)
                # vertexCount = insert!(vertices, vertexCount + 1, intersection, vertexCount)
                insertionIndex2 = orderIntersections(jj, vertList2, vertListFull2, vertMap2, vertices, intersection)
                insert!(labels2, insertionIndex2, mapNo, vertNo2)
                insert!(wrapPoints2, insertionIndex2, 0, vertNo2)
                vertNo2 = insert!(vertListFull2, insertionIndex2, vertexCount, vertNo2)
                updateVertMap(vertMap2, insertionIndex2)

                # temp = inOrOutBoundCheck(edge1,edge2)
                insert!(inOrOutBound2, insertionIndex2, inbound2, vertNo2)
                insert!(inOrOutBound1, insertionIndex1, inbound1, vertNo1)

                mapNo += 1
            end
        end
    end

    # set interior/exterior status for vertices which have not already been set
    # an even/odd number of intersections of ray originating at a point with edges of a polygon means the point is outside/inside the polygon respectively
    for i = 1:vertNoOrig1
        if vertList1[i] != -1 && vertCheck1[i]
            if mod(intNo1[i],2) == 0 #exterior
                labels1[vertMap1[i]] = 0
            else #interior
                labels1[vertMap1[i]] = -1
            end
        end
    end

    for i = 1:vertNoOrig2
        if vertList2[i] != -1 && vertCheck2[i]
            if mod(intNo2[i],2) == 0 #exterior
                labels2[vertMap2[i]] = 0
            else #interior
                labels2[vertMap2[i]] = -1
            end
        end
    end

    replace!(labelstemp1, 1:vertNo1, labels1)
    replace!(labelstemp2, 1:vertNo2, labels2)

    # for each intersection, set the value of the label to be the index of corresponding point in the other polygon
    for i = 1:mapNo-1

        # p1 = findall(labels1 .== i)
        # p2 = findall(labels2 .== i)
        index1 = 0
        index2 = 0

        for j in eachindex(labels1)
            if labels1[j] == i
                index1 = j
                break
            end
        end

        for j in eachindex(labels2)
            if labels2[j] == i
                index2 = j
                break
            end
        end

        if index1 > 0 && index2 > 0
            labelstemp1[index1] = index2
            labelstemp2[index2] = index1
        else
            error("intersection is pissing matching pair")
        end

    end


    inds = findall(wrapPoints1.>0)
    wrapPoints1[inds[1]] = 1
    for j = 2:length(inds)
        wrapPoints1[inds[j]] = inds[j-1] + 1
    end

    inds = findall(wrapPoints2.>0)
    wrapPoints2[inds[1]] = 1
    for j = 2:length(inds)
        wrapPoints2[inds[j]] = inds[j-1] + 1
    end

    data.relevantVertexNo = vertexCount

    # update clipped/subject polygons
    if clippingPoly == 1
        surface1.clippedGeometry.tempVertexNo = vertNo1
        surface2.clippedGeometry.tempVertexNo = vertNo2
    elseif clippingPoly == 2
        surface2.clippedGeometry.tempVertexNo = vertNo1
        surface1.clippedGeometry.tempVertexNo = vertNo2
    end

    replace!(labels1, 1:vertNo1, labelstemp1)
    replace!(labels2, 1:vertNo2, labelstemp2)



end

"""
    Checks which polygon is on top
    returns 1 if the first is on top, 2 if the second, and 0 if no intersection is found
"""
function checkPolyDistance(vertices, intersection, surface1, surface2, edge)

    intNo1 = surface1.clippedGeometry.intersectionNumber
    intNo2 = surface2.clippedGeometry.intersectionNumber

    vertCheck1 = surface1.clippedGeometry.vertCheck
    vertCheck2 = surface2.clippedGeometry.vertCheck

    vertList1 = surface1.vertexIndices
    vertList2 = surface2.vertexIndices

    vertNo1 = surface1.vertexNo
    vertNo2 = surface2.vertexNo

    edge1 = surface1.clippedGeometry.edge
    edge2 = surface2.clippedGeometry.edge

    
    # reset variables
    vertCheck1[1:vertNo1] .= true
    vertCheck2[1:vertNo2] .= true

    intNo1[1:vertNo1] .= 0
    intNo2[1:vertNo2] .= 0
    
    startPoint = [1;1]

    for ii = 1:vertNo1 - 1
        for jj = 1:vertNo2 - 1

            # get the next edge from each polygon
            if vertList1[ii] != -1
                nextEdge!(edge1, vertices, vertList1, ii, startPoint[1])
            else
                # wrapPoints[1][ii-1] = 1
                startPoint[1] = ii + 1
                break # has to be break to exit the inner loop entirely and increment ii by 1
            end

            if vertList2[jj] != -1
                nextEdge!(edge2, vertices, vertList2, jj, startPoint[2])
            else
                # wrapPoints[2][jj-1] = 1
                startPoint[2] = jj + 1
                continue # continue increments the inner loop by 1
            end

            # check intersection for the current pair of edges
            intersectionExists, flag1, flag2 = lineIntersection!(edge1, edge2, intersection[1], vertCheck1[ii], vertCheck2[jj])
                    
            # set X,Y coordinates of 2nd intersection to be the same as the first
            intersection[2][:] = intersection[1]

            if intersectionExists 

                fullPoint!(intersection[1], edge1)

                fullPoint!(intersection[2], edge2)

                if intersection[1][3] > intersection[2][3]
                    return 1
                elseif intersection[1][3] < intersection[2][3]
                    return 2
                else
                    # print("facets appear to intersect")
                end

            elseif flag2 == -3

                fullPoint!(intersection[1], edge1)

                if intersection[1][3] > vertices[vertList2[jj]][3]
                    return 1
                elseif intersection[1][3] < vertices[vertList2[jj]][3]
                    return 2
                else
                    vertCheck2[jj] = false
                    # error("facets appear to intersect")
                end

            elseif flag1 == -3

                fullPoint!(intersection[2], edge2)

                if intersection[2][3] < vertices[vertList1[ii]][3]
                    return 1
                elseif intersection[2][3] > vertices[vertList1[ii]][3]
                    return 2
                else
                    vertCheck1[ii] = false
                    # error("facets appear to intersect")
                end

            elseif flag1 == -2 && flag2 == -2

                if vertices[vertList1[ii]][3] > vertices[vertList2[jj]][3]
                    return 1
                elseif vertices[vertList1[ii]][3] < vertices[vertList2[jj]][3]
                    return 2
                else
                    vertCheck1[ii] = false
                    vertCheck2[jj] = false
                    # error("facets appear to intersect")
                end

            end

            # edge one 
            if flag1 == 1
                # increment number of ray intersections associated with each point (used later to check interior/exterior status)
                intNo1[ii] += flag1
            elseif flag1 == -1
                # unknown number of intersections due to edges being parallel
                # directly check if the point is interior or exterior 
                check = interiorPointCheck(edge, vertList1[ii], vertList2, vertNo2, vertices)
                if check == -1
                    z2 = fullPointPlane(vertices[vertList1[ii]], @view vertices[vertList2[1:3]])
                    if z2 > vertices[vertList1[ii]][3]
                        return 2
                    elseif z2 < vertices[vertList1[ii]][3]
                        return 1
                    else
                        # print("facets appear to intersect")
                    end
                else
                    vertCheck1[ii] = false
                end
            end

            # edge two
            if flag2 == 1
                # increment number of ray intersections associated with each point (used later to check interior/exterior status)
                intNo2[jj] += flag2
            elseif flag2 == -1
                # unknown number of intersections due to edges being parallel
                # directly check if the point is interior or exterior 
                check = interiorPointCheck(edge, vertList2[jj], vertList1, vertNo1, vertices)
                if check == -1
                    z1 = fullPointPlane(vertices[vertList2[jj]], @view vertices[vertList1[1:3]])
                    if z1 > vertices[vertList2[jj]][3]
                        return 1
                    elseif z1 < vertices[vertList2[jj]][3]
                        return 2
                    else
                        # print("facets appear to intersect")
                    end
                else
                    vertCheck2[jj] = false
                end
            end

        end
    end

    # set interior/exterior status for vertices which have not already been set
    for i = 1:vertNo1
        if vertCheck1[i]
            if !(mod(intNo1[i],2) == 0)
                z = fullPointPlane(vertices[vertList1[i]], @view vertices[vertList2[1:3]])
                if z > vertices[vertList1[i]][3]
                    return 1
                elseif z < vertices[vertList1[i]][3]
                    return 2
                else
                    # print("facets appear to intersect")
                end
            end
        end
    end

    for i = 1:vertNo2
        if vertCheck2[i]
            if !(mod(intNo2[i],2) == 0)
                z = fullPointPlane(vertices[vertList2[i]], @view vertices[vertList1[1:3]])
                if z > vertices[vertList2[i]][3]
                    return 1
                elseif z < vertices[vertList2[i]][3]
                    return 2
                else
                    # print("facets appear to intersect")
                end
            end
        end
    end

    return 0
end

"""
    lineIntersection!(line1, line2, intersection, check1, check2)
    
    This function finds the intersection of two 2D lines and checks if the intersection is on a specified line segements and ray

    Coordinates are stored as vectors. This function only considers the first and second elements of a vector as the X and Y coordinates. Coordinate vectors can generaly longer than 2 elements without issue, and the function is intended to be compatible with 3 vectors where the z coordinate (3rd element) will be determined later

    line1 and line2 should be vectors of coordinate vectors corresponding to the first and second points of each line segment

    intersection is a coordinate vector that the function writes to in-place

    check1 and check2 determine if the function should check for ray/line intersections in addition to line/line intersections check 1 determines if the ray1/line2 intersection should be considered where ray1 is a hypothetical ray eminating from the first point of line1 in the direction of the second point of line1


"""
function lineIntersection!(line1, line2, intersection, check1, check2, prevP1 = nothing, prevP2 = nothing, tol = 0)

    tol1 = 1e-6
    tol2 = 1e-8
    # bounds of lines
    xmax1 = max(line1[1][1], line1[2][1])
    xmax2 = max(line2[1][1], line2[2][1])
    ymax1 = max(line1[1][2], line1[2][2])
    ymax2 = max(line2[1][2], line2[2][2])

    xmin1 = min(line1[1][1], line1[2][1])
    xmin2 = min(line2[1][1], line2[2][1])
    ymin1 = min(line1[1][2], line1[2][2])
    ymin2 = min(line2[1][2], line2[2][2])

    # alculate slope
    m1 = (line1[2][2] - line1[1][2])/(line1[2][1] - line1[1][1])
    m2 = (line2[2][2] - line2[1][2])/(line2[2][1] - line2[1][1])

    # check if points are on the other lines
    if !isnan(m2)
        L1P1onL2 = (abs( (line1[1][2] - line2[1][2]) - m2*(line1[1][1] - line2[1][1]) ) < tol1)
        L1P2onL2 = (abs( (line1[2][2] - line2[1][2]) - m2*(line1[2][1] - line2[1][1]) ) < tol1)
    else
        L1P1onL2 = (abs(line1[1][1] - line2[1][1]) < tol1)
        L1P2onL2 = (abs(line1[2][1] - line2[1][1]) < tol1)
    end

    if !isnan(m2)
        L2P1onL1 = (abs( (line2[1][2] - line1[1][2]) - m1*(line2[1][1] - line1[1][1]) ) < tol1)
        L2P2onL1 = (abs( (line2[2][2] - line1[1][2]) - m1*(line2[2][1] - line1[1][1]) ) < tol1)
    else
        L2P1onL1 = (abs(line1[1][1] - line2[1][1]) < tol1)
        L2P2onL1 = (abs(line1[1][1] - line2[2][1]) < tol1)
    end

    if (L1P1onL2 && L1P2onL2) || (L2P1onL1 && L2P2onL1)
        # check if lines are colinear 

        # there is not a unique intersection point
        segmentIntersection = false

        # shared base point
        if abs(line1[1][1] - line2[1][1]) < tol1
            segmentIntersection = false
            ray1Intersections = -2
            ray2Intersections = -2
            return segmentIntersection, ray1Intersections, ray2Intersections, 0, 0
        end

        # set both lines to do interior/exterior check if not already done
        segmentIntersection = false
        ray1Intersections = setIntersections(-1, check1)
        ray2Intersections = setIntersections(-1, check2)

        # check if base points are on the lines (not including end points)
        if  line1[1][1] - xmin2 > tol1 && xmax2 - line1[1][1] > tol1
            ray1Intersections = -3
        end

        if line2[1][1] - xmin1 > tol1 && xmax1 - line2[1][1] > tol1
            ray2Intersections = -3
        end

        return segmentIntersection, ray1Intersections, ray2Intersections, 0, 0

    elseif abs(m1 - m2) < 1e-5 || (isequal(m1,m2) && abs(line1[1][1] - line2[1][1]) > 1e-5)

        # non colinear parallel lines do not intersect
        ray2Intersections = 0 
        ray1Intersections = 0
        segmentIntersection = false
        return segmentIntersection, ray1Intersections, ray2Intersections, -3, -3
    end

    # check for shared base points
    if L1P1onL2 && L2P1onL1 #vectorEquality2D(line1[1], line2[1], tol1) #line1[1] ≈ line2[1]

        segmentIntersection = false
        ray1Intersections = -2
        ray2Intersections = -2
        if !isnothing(prevP1) && !isnothing(prevP1)
            OutBound1, OutBound2 = inOrOutBoundCheckSharedVertex(line1,line2, prevP1, prevP2, tol1)
        else
            OutBound1 = -3
            OutBound2 = -3
        end
        return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
    end

    # compute the intersection of the lines
    if isnan(m1)
        # if only line 1 is vertical
        intersection[1] = line1[1][1]
        intersection[2] = m2*(intersection[1] - line2[1][1]) + line2[1][2]
    elseif isnan(m2)
        # if only line 2 is vertical
        intersection[1] = line2[1][1]
        intersection[2] = m1*(intersection[1] - line1[1][1]) + line1[1][2]
    else
        # coordinates of intersection
        intersection[1] = (-line2[1][1]*m2 + line1[1][1]*m1 + line2[1][2] - line1[1][2])/(m1 - m2)
        intersection[2] = m1*(intersection[1] - line1[1][1]) + line1[1][2]
    end

    # check if intersection is on line segement 1
    line1Check = (!L1P1onL2 && !L1P2onL2 && pointLineBoundsCheck(intersection, xmax1, xmin1, ymax1, ymin1, -tol2))
    # check if intersection is on line segement 2
    line2Check = (!L2P1onL1 && !L2P2onL1 && pointLineBoundsCheck(intersection, xmax2, xmin2, ymax2, ymin2, -tol2))

    dx1 = abs(line1[1][1] - intersection[1])
    dy1 = abs(line1[1][2] - intersection[2])
    # check if intersection is on ray 1
    if dx1 > tol2 && dy1 > tol2
        # if intersection is sufficiently far from the base point that the dot product check will be stable
        ray1check = lineDirectionCheck(line1, intersection)
    elseif dx1 < tol2 && abs(line1[1][1] - line1[2][1]) < tol2 && dy1 > tol2
        # deal with the case of a (near) veritcal line
        ray1check = lineDirectionCheck(line1, intersection) 
    elseif dy1 < tol2 && abs(line1[1][2] - line1[2][2])< 1e-5 && dx1 > tol2
        # deal with the case of a (near) horizontal line
        ray1check = lineDirectionCheck(line1, intersection) 
    else
        ray1check = false
    end

    dx2 = abs(line2[1][1] - intersection[1])
    dy2 = abs(line2[1][2] - intersection[2])
    # check if intersection is on ray 2
    if dx2 > tol2 && dy2 > tol2
        # if intersection is sufficiently far from the base point that the dot product check will be stable
        ray2check = lineDirectionCheck(line2, intersection) 
    elseif dx2 < tol2 &&  abs(line2[1][1] - line2[2][1]) < tol2 && dy2 > tol2
        # deal with the case of a (near) veritcal line
        ray2check = lineDirectionCheck(line2, intersection) 
    elseif dy2 < tol2 && abs(line2[1][2] - line2[2][2])< tol2 && dx2 > tol2
        # deal with the case of a (near) horizontal line
        ray2check = lineDirectionCheck(line2, intersection) 
    else
        ray2check = false
    end

    if line1Check && line2Check 
        # intersection is in the range of both line segments (excluding end/base points)
        ray1Intersections = setIntersections(1, check1)
        ray2Intersections = setIntersections(1, check2)
        segmentIntersection = true
        OutBound2 = inOrOutBoundCheck(line1,line2, 1e-4)
        OutBound1 = -OutBound2
        return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
    elseif !line1Check && !line2Check && !ray1check && !ray2check && !L1P1onL2 && !L1P2onL2 && !L2P1onL1 && !L2P2onL1
        # no relevant intersection
        ray2Intersections = 0 
        ray1Intersections = 0
        segmentIntersection = false
        OutBound1 = -3
        OutBound2 = -3
        return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
    end


    if L1P1onL2
        # base point of line1 is on line 2 in the range of the line segment
        segmentIntersection = false
        ray1Intersections = setIntersections(-3, line2Check)
        ray2Intersections = setIntersections(-1, check2 && ray2check)
        if !isnothing(prevP1) && line2Check
            OutBound1, OutBound2 = inOrOutBoundCheckVertInt(line1, line2, prevP1, 1e-4) 
        else
            OutBound1 = -3
            OutBound2 = -3
        end
        return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
    elseif check1 && ray1check
        if line2Check
            # end point of line1 is on line 2 in the range of the line segment
            ray1Intersections = 1
        elseif L2P1onL1 || L2P2onL1
            # if the intersection is shared with base or end point of line2
            ray1Intersections = -1
        else
            ray1Intersections = 0
        end
    else
        ray1Intersections = 0 
    end

    if L2P1onL1
        # base point of line1 is on line 2 in the range of the line segment
        segmentIntersection = false
        ray2Intersections = setIntersections(-3, line1Check)
        ray1Intersections = setIntersections(-1, check1 && ray1check)
        if !isnothing(prevP2) && line1Check
            OutBound2, OutBound1 = inOrOutBoundCheckVertInt(line2, line1, prevP2, 1e-4)
        else
            OutBound1 = -3
            OutBound2 = -3
        end
        return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
    elseif check2 && ray2check
        if line1Check
            # end point of line1 is on line 2 in the range of the line segment
            ray2Intersections = 1
        elseif L1P1onL2 || L1P2onL2
            # if the intersection is shared with base or end point of line2
            ray2Intersections = -1
        else
            ray2Intersections = 0
        end
    else
        ray2Intersections = 0
    end

    segmentIntersection = false
    OutBound1 = -3
    OutBound2 = -3

    return segmentIntersection, ray1Intersections, ray2Intersections, OutBound1, OutBound2
end

"""
    check if a ray intersects a line segment
    ray is given by a base point and a second point which can be anywhere along the ray and just specifies the direction
    it is assumed that the ray is not vertical (slope exists), although the line segment can be vertical
"""
function rayLineIntersection(ray, line, tol = 0.0)

    # ray slope
    if ray[2][1] - ray[1][1] == 0
        return -1
    else
        m1 = (ray[2][2] - ray[1][2])/(ray[2][1] - ray[1][1])
    end

    if line[1][1] == line[2][1] # if the line is vertical
        int1 = line[1][1]
        int2 = m1*(int1 - ray[1][1]) + ray[1][2]
    else
        m2 = (line[2][2] - line[1][2])/(line[2][1] - line[1][1])

        if isequal(m1, m2)
            return -1
        else
            int1 = (-line[1][1]*m2 + ray[1][1]*m1 + line[1][2] - ray[1][2])/(m1 - m2)
            int2 = m1*(int1 - ray[1][1]) + ray[1][2]
        end
    end

    xmax = max(line[1][1], line[2][1])
    xmin = min(line[1][1], line[2][1])
    ymax = max(line[1][2], line[2][2])
    ymin = min(line[1][2], line[2][2])

    
    if int1 - xmin > tol && int1 - xmax < -tol && int2 - ymin > tol && int2 - ymax < -tol && sign(dot2_unsafe(ray[2] - ray[1], [int1;int2;0] - ray[1])) == 1 
        return 1
    else
        return 0
    end
end

"""
    check if a 2D point is inside a polygon (shape)
    point is a 2x1 vector
    shape is a 1xn vector, where each element is the index of the corresponding vertex (ordered clockwise)
    vertices is a 1xn vector of 2x1 vectors which are the x,y coordinates of each vertex
"""
function interiorPointCheck(edge, pointIndex, shape, maxIndex, vertices)

    # generate ray that isn't vertical and isn't parallel to any of the edges of the shape
    invalid = true

    j = 0

    point = @view vertices[pointIndex]

    # while loop to retry procedure until a ray that is not vertical or parallel to side is generated
    while invalid

        invalid = false
        # preventing infinite loop
        j += 1
        if j > 1000
            error("maximum iterations reached without finding valid ray")
        end

        # make second point very large to ensure it will be outside of the object
        ray = [point; [rand(3)*10^4]]

        # counting number of intersections with shape edges to determine if the point is interior or exterior
        intersectionCount = 0

        for ii = 1:maxIndex-1
            # get the next edge from each polygon
            if shape[ii] != -1
                nextEdge!(edge, vertices, shape, ii)
            else
                continue
            end

            if point == edge[1] || point == edge[2]
                return -2
            end
            
            # check if the ray and current edge intersect
            check = rayLineIntersection(ray, edge)

            if check == -1
                # restart if ray is vertical or parallel to a side
                invalid = true
                break
            elseif check == 1
                intersectionCount += check
            else

            end

        end

        # if the ray was not valid then restart the loop
        if invalid
            continue
        end

        # if loop successfully completes, then the ray is not vertical or parallel to any side and is valid
        # check if the intersection count is odd or even to determine if the point is interior or exterior
        if mod(intersectionCount,2) == 0
            # exterior
            return 0
        else
            #interior
            return -1
        end
    end
end

function unpacking(surface1in, surface2in, clipNo)

    # select appropriate variables based on which polygon is the clipping polygon (the other is the subject polygon)
    if clipNo == 1
        surface1 = surface1in
        surface2 = surface2in
    elseif clipNo == 2
        surface2 = surface1in
        surface1 = surface2in
    else
        error("???")    
    end

    # number of relevant vertices in the vertices lists 
    vertNoOrig1 = surface1.vertexNo
    vertNoOrig2 = surface2.clippedGeometry.vertexNo

    vertNo1 = surface1.vertexNo
    vertNo2 = surface2.clippedGeometry.vertexNo

    # list of original vertices
    vertList1 = surface1.vertexIndices
    vertList2 = surface2.clippedGeometry.vertexIndices
 
    # list of all vertices including intersection points
    vertListFull1 = surface1.clippedGeometry.tempVertexIndices
    vertListFull2 = surface2.clippedGeometry.tempVertexIndices

    # label of each point in the Full vertex list
    labels1 = surface1.clippedGeometry.labels
    labels2 = surface2.clippedGeometry.labels
    labelstemp1 = surface1.clippedGeometry.labelstemp
    labelstemp2 = surface2.clippedGeometry.labelstemp

    # contains flags indiciating if intersections are inbound or outbound
    inOrOutBound1 = surface1.clippedGeometry.inOrOutBound
    inOrOutBound2 = surface2.clippedGeometry.inOrOutBound

    # number of intersections of the nth point in the original vertex list
    intNo1 = surface1.clippedGeometry.intersectionNumber
    intNo2 = surface2.clippedGeometry.intersectionNumber

    # determines whether the line intersection function should also check for ray intersections for the first shape
    vertCheck1 = surface1.clippedGeometry.vertCheck
    vertCheck2 = surface2.clippedGeometry.vertCheck

    # maps indexs in the original vertex list to indexes in the full vertex list
    vertMap1 = surface1.clippedGeometry.vertMap
    vertMap2 = surface2.clippedGeometry.vertMap

    # arrays to store the current edge
    edge1 = surface1.clippedGeometry.edge
    edge2 = surface2.clippedGeometry.edge

    # keep track of the points where the vertexLists wrap around (only the clipped polygon can have wrapping) 
    wrapPoints1 = surface1.clippedGeometry.tempWrapPoints
    wrapPoints2 = surface2.clippedGeometry.tempWrapPoints

    # reset variables
    replace!(vertCheck1, 1:vertNoOrig1, true)
    replace!(vertCheck2, 1:vertNoOrig2, true)

    replace!(intNo1, 1:vertNoOrig1, 0)
    replace!(intNo2, 1:vertNoOrig2, 0)

    #reset wrapPoints
    replace!(wrapPoints1, 1:vertNo1, surface1.wrapPoints[1:vertNo1])
    replace!(wrapPoints2, 1:vertNo2, surface2.clippedGeometry.wrapPoints[1:vertNo2])

    # reset temporary vertex list that will ahve intersections added to it
    copyVertices!(vertListFull1, vertList1, vertNo1)
    copyVertices!(vertListFull2, vertList2, vertNo2)

    # reset labels
    resetValues!(labels1, vertNo1, -2)
    resetValues!(labels2, vertNo2, -2)

    resetValues!(inOrOutBound1, vertNo1, -2)
    resetValues!(inOrOutBound2, vertNo2, -2)

    resetVertMap!(vertMap1) 
    resetVertMap!(vertMap2) 

    return vertNoOrig1, vertNoOrig2, vertNo1, vertNo2, vertList1, vertList2, vertListFull1, vertListFull2, labels1, labels2, labelstemp1, labelstemp2, inOrOutBound1, inOrOutBound2, intNo1, intNo2, vertCheck1, vertCheck2, vertMap1, vertMap2, edge1, edge2, wrapPoints1, wrapPoints2
end     