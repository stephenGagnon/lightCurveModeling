"""
if the scenario type is Nothing, then self shadoing should not be done, and only the line of sight checks are done
"""
function selfShadowing!(att :: anyAttitude, scen :: scenario{Nothing}, obsNo)
    for i = 1:scen.object.surfaceNo
        check1 = dot(scen.sun.bodyFrameUnitVec, scen.object.surfaces[i].u_n) > 0
        check2 = dot(scen.observers[obsNo].bodyFrameUnitVec, scen.object.surfaces[i].u_n) > 0
        scen.object.surfaces[i].isVisible = check1 && check2
        scen.object.surfaces[i].areaClipped = scen.object.surfaces[i].Area
    end
end

"""
    computes the area of each facet which is visible to both a specified observer and the sun
    calculations done in the body frame
    assumes the sun and observer body frame vectors are already calculated based on some rotation
"""
function selfShadowing!(att :: anyAttitude, scen :: scenario{objectProjectionData}, obsNo)

    # general procedure is to project the vertices of visible facets onto the plane normal to the sun and find any intersections of these 'shadows' in order to clip the 2D projected polygons. Then deproject the clipped polygons and re-project them onto the plane normal to the observer and perform the clipping process again to find the final area of the facet visibile to both the observer and sun.

    # temporary variable for readability 
    data = scen.object.projectionData

    # reset number of projected vertices
    data.relevantVertexNo = data.originalVertexNo

    # reset isProjected Flag
    data.isProjected .= false

    # number of surfaces
    surfaceNo = scen.object.surfaceNo

    # reset clipped polygons
    for s in scen.object.surfaces
        replace!(s.clippedGeometry.vertexIndices, eachindex(s.vertexIndices), s.vertexIndices)
        s.clippedGeometry.vertexNo = s.vertexNo
        s.clippedGeometry.wrapPoints .= 0
        s.clippedGeometry.wrapPoints[s.vertexNo] = 1
    end

    R = scen.sun.rotation * q2A(att)' 

    # for each facet: check if it is visible and then if so, project it onto the plane normal to the sun
    for i = 1:surfaceNo

        # check line of sight to see if facet is visible and store that as a Bool
        check1 = dot(scen.sun.bodyFrameUnitVec, scen.object.surfaces[i].u_n) > 0
        # check line of sight to see if facet is visible and store that as a Bool
        check2 = dot(scen.observers[obsNo].bodyFrameUnitVec, scen.object.surfaces[i].u_n) > 0
        scen.object.surfaces[i].isVisible = check1 && check2

        # if a facet is visible project all of the vertices associated with the facet onto the plane normal to the sun direction
        if scen.object.surfaces[i].isVisible

            # reset bounds on surfaces
            data.projectedBounds[i] .= NaN

            for j = 1:scen.object.surfaces[i].vertexNo - 1 #in scen.object.surfaces[i].vertexIndices
                currentIndex = scen.object.surfaces[i].vertexIndices[j]
                #if the vertex has not already been projected 
                if !data.isProjected[currentIndex]
                    # rotate the point into a frame where the sun lies along the x axis so that the y,z coordinates are the projected coordinates in the normal plane and the x coordinate is the distance to the normal plane
                    data.vertices[currentIndex] = R * scen.object.vertices[currentIndex]
                    data.isProjected[currentIndex] = true
                end
                # update the bounding x,y coordinates for the current surface
                XYextrema!((data.vertices[currentIndex]), data.projectedBounds[i])
            end
        end
    end

    # visibleFacets = findall(getfield.(scen.object.surfaces, :isVisible))

    resetCheck = false

    # loop through all pairs of visible surfaces
    for ii = 1:surfaceNo-1, jj = ii+1:surfaceNo

        # if the bounding boxes of the two projected surfaces overlap
        if scen.object.surfaces[ii].isVisible && 
            scen.object.surfaces[jj].isVisible && 
            boundsCheck(data.projectedBounds[ii], data.projectedBounds[jj])

            # if ii == 5 && jj == 12
            #     print("")
            # end

            # find the intersection of the two projected surfaces (if it exists)
            # updates the clippedPolygon field of each surface to be a list of vertices corresponding to the unshadowed and visible area 
            if obsNo == 2 && ii == 9 && jj == 11 
                diagnostics = true
                # plotClippedSat(scen)
                # print("")
            else
                diagnostics = false
            end

            flag = polygonClipping!(data, scen.object.surfaces[ii], scen.object.surfaces[jj], diagnostics)
            if flag != -2 && flag != -1
                resetCheck = true
            elseif flag == -1
                return -1
            end
            # plotShape2D(scen.object.surfaces[ii].clippedGeometry.vertexIndices, data.vertices, :green)
            # plotShape2D!(scen.object.surfaces[jj].clippedGeometry.vertexIndices, data.vertices, :red)
            # print("")
        end
    end

    # plotClippedSat(scen)
    # print("")

    # get updated visibility list if facets are completely shaddowed from the sun
    # visibleFacets = findall(getfield.(scen.object.surfaces, :isVisible))

    if data.relevantVertexNo - length(data.isProjected) > 0
        append!(data.isProjected, Vector{Bool}(undef, data.relevantVertexNo - length(data.isProjected)) ) 
    end

    data.isProjected .= false
    R = scen.observers[obsNo].rotation * scen.sun.rotation'

    # de-project and reproject for observer
    # for each facet: check if it is visible and then if so, project it onto the plane normal to the observer
    for i = 1:surfaceNo

        # if a facet is visible project all of the vertices associated with the facet onto the plane normal to the sun direction
        if scen.object.surfaces[i].isVisible

            # reset bounds on surfaces
            data.projectedBounds[i] .= NaN

            # project any points in original polygon that are not also in the clipped polygon
            for j = 1:scen.object.surfaces[i].vertexNo - 1 #in scen.object.surfaces[i].vertexIndices
                currentIndex = scen.object.surfaces[i].vertexIndices[j]
                if !data.isProjected[currentIndex]
                    # rotate the point into a frame where the observer lies along the x axis so that the y,z coordinates are the projected coordinates in the normal plane and the x coordinate is the distance to the normal plane
                    data.vertices[currentIndex] =  R * data.vertices[currentIndex]
                    data.isProjected[currentIndex] = true
                end
                # update the bounding x,y coordinates for the current surface
                XYextrema!((data.vertices[currentIndex]), data.projectedBounds[i])
            end

            # if some addional vertices were added when accounting for light source shadowing, project those vertices
            vertNo =  scen.object.surfaces[i].clippedGeometry.vertexNo
            vertexIndices = scen.object.surfaces[i].clippedGeometry.vertexIndices[1:vertNo] 

            if resetCheck
                # project all the points in the clipped polygon
                for currentIndex in vertexIndices #vertexIndices[findall(vertexIndices[1:vertNo] .> 0)]  
                    # currentIndex = scen.object.surfaces[i].clippedGeometry.vertexIndices[j]
                    #if the vertex has not already been projected 
                    if currentIndex > 0 && !data.isProjected[currentIndex] 
                        # rotate the point into a frame where the observer lies along the x axis so that the y,z coordinates are the projected coordinates in the normal plane and the x coordinate is the distance to the normal plane
                        data.vertices[currentIndex] =  R * data.vertices[currentIndex]
                        data.isProjected[currentIndex] = true

                        # update the bounding x,y coordinates for the current surface
                        XYextrema!((data.vertices[currentIndex]), data.projectedBounds[i])
                    end
                end
            end
        end
    end

    # plotClippedSat(scen)
    # print("")
    diagnostics = false

    # clipping for observer
    # loop through all pairs of visible surfaces
    # loop through all pairs of visible surfaces
    for ii = 1:surfaceNo-1, jj = ii+1:surfaceNo

        # if the bounding boxes of the two projected surfaces overlap
        if scen.object.surfaces[ii].isVisible && 
            scen.object.surfaces[jj].isVisible && 
            boundsCheck(data.projectedBounds[ii], data.projectedBounds[jj])

            if obsNo == 2 && ii == 4 && jj == 11
                diagnostics = true
                # plotClippedSat(scen)
                print("")
            end
            # find the intersection of the two projected surfaces (if it exists)
            # updates the clippedPolygon field of each surface to be a list of vertices corresponding to the unshadowed and visible area 
            flag = polygonClipping!(data, scen.object.surfaces[ii], scen.object.surfaces[jj],diagnostics)
        end
    end

    # plotClippedSat(scen)
    # print("")

    for i = 1:scen.object.surfaceNo
        if scen.object.surfaces[i].isVisible
            # clippedVertices = @view data.vertices[scen.object.surfaces[i].clippedGeometry.vertexIndices]
            # projectedVertices = @view data.vertices[scen.object.surfaces[i].vertexIndices]
            surfc = scen.object.surfaces[i].clippedGeometry

            areaClipped = ShoelaceAlgorithm(surfc.vertexIndices, data.vertices, surfc.wrapPoints, surfc.vertexNo)
            areaProjected = ShoelaceAlgorithm(scen.object.surfaces[i].vertexIndices, data.vertices, scen.object.surfaces[i].vertexNo) 
            ratio = areaClipped/areaProjected

            if surfc.vertexIndices == scen.object.surfaces[i].vertexIndices && ratio - 1 > 1e-8
                print(ratio, "\n")
                error("????????")
            elseif ratio - 1 > 1e-8
                print(ratio, "\n")
                error("????????")
            end

            scen.object.surfaces[i].areaClipped = scen.object.surfaces[i].Area * ratio
        end
    end
    return 0

end

