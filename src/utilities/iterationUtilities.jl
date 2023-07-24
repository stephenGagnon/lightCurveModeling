function nextEdge(vertices, vertList, index, startPoint)
    if vertList[index + 1] == -1
        return vertices[vertList[[index, startPoint]]]
    else
        return vertices[vertList[index:index+1]]
    end
end

function nextEdge!(edge, vertices, vertList, index, startPoint)
    if vertList[index + 1] == -1
        edge[1] = vertices[vertList[index]]
        edge[2] = vertices[vertList[startPoint]]
        # return vertices[vertList[[index, startPoint]]]
    else
        edge[1] = vertices[vertList[index]]
        edge[2] = vertices[vertList[index + 1]]
        # return vertices[vertList[index:index+1]]
    end
end

function nextEdge(vertices, vertList, index)
    if vertList[index + 1] == -1
        startPoint = 1
        for i = index:-1:1
            if vertList[i] == -1
                startPoint = i + 1
                break
            end
        end
        return vertices[vertList[[index, startPoint]]]
    else
        return vertices[vertList[index:index+1]]
    end
end

function nextEdge!(edge, vertices, vertList, index)

    if vertList[index + 1] == -1
        nextInd = 1
        for i = index:-1:2
            if vertList[i-1] == -1
                nextInd = i
                break
            end
        end
    else
        nextInd = index + 1
    end


    if index == 1 || vertList[index - 1] == -1
        prevInd = 0
        for i = index:length(vertList)-1
            if vertList[i+1] == -1
                prevInd = i
                break
            end
        end
    else
        prevInd = index - 1
    end

    previousPoint = vertices[vertList[prevInd]]
    edge[1] = vertices[vertList[index]]
    edge[2] = vertices[vertList[nextInd]]

    return previousPoint
    # edge[1] = vertices[vertList[index]]
    # edge[2] = vertices[vertList[startPoint]]
end

function nextPoint(vertList :: Vector, index :: Int64, startPoint :: Int64)
    if vertList[index + 1] == -1
        return startPoint
    else
        return index + 1
    end
end

"""
    finds the next point while raversing the boundary of the clipped polygon region. 
    updates the point variable in place to contain the next point
    returns an integer flag:
     2: currentPoint has already been reached, but it is an intersection, and the matching intersection on the other shape has not been reached
     1: two valid next points were found
     0: next point found
    -1: no valid point found
    -2: next point is the starting point (handled in the higher level function)
    -3: current point is already used
    
"""
# function nextPoint!(point, currentPoint :: Vector, labels :: Vector{Vector{Int64}}, wrapPoints :: Vector{Vector{Int64}}, startingPoint :: Vector{Int64}, inOrOutBound, usedPoints)

    #     flag = _nextPoint!(point, currentPoint, startingPoint, labels, wrapPoints, inOrOutBound, usedPoints)

    #     # outer wrapper to check if current point has already been reached
    #     if usedPoints[point[1]][point[2]]
    #         label = labels[point[1]][point[2]]
    #         # if the current point has been reached but the other point has not been reached
    #         if label > 0 && usedPoints[other12(point[1])][label] == false
    #             return 2
    #         else
    #             return -3
    #         end
    
    #     else
    #         return flag
    #     end

    # end

function nextPoint!(point, currentPoint :: Vector, labels :: Vector{Vector{Int64}}, wrapPoints :: Vector{Vector{Int64}}, startingPoint, inOrOutBound, usedPoints)

    label = labels[currentPoint[1]][currentPoint[2]]

    if label > 0 # point is an intersection

        # checks the current polygon and gets the inbound/outBound status of the intersection for both the clipping and subject polygon
        if currentPoint[1] == 1
            # the IO status of the clipping polygon is determined by the next point since the edges are evaluated counter clockwise but we proceed clockwise
            nextClipIndex = nextIndex(currentPoint[2], wrapPoints[1], -1)
            nextSubIndex = nextIndex(label, wrapPoints[2], 1)

            IO_clip = inOrOutBound[1][nextClipIndex]
            IO_sub = inOrOutBound[2][label]
            label_clip = labels[1][nextClipIndex]
            label_sub = labels[2][nextSubIndex]

        elseif currentPoint[1] == 2 
            # the IO status of the subject polygon is determined by the current point 
            nextClipIndex = nextIndex(label, wrapPoints[1], -1)
            nextSubIndex = nextIndex(currentPoint[2], wrapPoints[2], 1)

            IO_sub = inOrOutBound[2][currentPoint[2]]
            IO_clip = inOrOutBound[1][nextClipIndex]
            label_clip = labels[1][nextClipIndex]
            label_sub = labels[2][nextSubIndex]
        end 

        # check if the next point on each polygon is valid
        valid_clip = label_clip == -1 || (label_clip > 0 && (IO_clip == 1 || IO_clip == 0))
        valid_sub = label_sub == 0 || (label_sub > 0 && (IO_sub == -1 || IO_sub == 0))

        if  startingPoint[1] == 1 && ((startingPoint[2] == nextClipIndex && valid_clip) || (labels[1][startingPoint[2]] == nextSubIndex && valid_sub))
            point[1] = 1
            point[2] = nextClipIndex
            return -2, ()
        elseif startingPoint[1] == 2 && ((startingPoint[2] == nextSubIndex && valid_sub) || (labels[2][startingPoint[2]] == nextClipIndex && valid_clip))
            point[1] = 2
            point[2] = nextSubIndex
            return -2, ()
        end

        valid_clip = !usedPoints[1][nextClipIndex] && valid_clip
        valid_sub = !usedPoints[2][nextSubIndex] && valid_sub
        
        if valid_sub && valid_clip
            # if both edges are valid

            if IO_clip == 1 && IO_sub == 0 
                # prefer to take the determinate edge over the indeterminate (IO == 0)
                point[1] = 1
                point[2] = nextClipIndex
                return 0, ()
            elseif IO_clip == 0 && IO_sub == -1
                # prefer to take the determinate edge over the indeterminate (IO == 0)
                point[1] = 2
                point[2] = nextSubIndex
                return 0,()
            elseif label_clip == -1 && label_sub > 0
                point[1] = 1
                point[2] = nextClipIndex
                return 0, ()
            elseif label_clip > 0 && label_sub == 0
                point[1] = 2
                point[2] = nextSubIndex
                return 0,()
            elseif IO_clip == 1 && IO_sub == -1 
                # neither edge is indeterminate
                # prefer to switch polygons
                point[1] = other12(currentPoint[1])

                if currentPoint[1] == 1
                    point[2] = nextSubIndex
                elseif currentPoint[1] == 2 
                    point[2] = nextClipIndex
                end
                return 0, ()
            elseif IO_sub == 0 && IO_clip == 0 && label_clip == nextSubIndex && label_sub == nextClipIndex
                # both edges are indeterminate with a shared intersection as the next point
                # prefer to move along the subject polygon  
                point[1] = 2
                point[2] = nextSubIndex
                return 1, ()
            else
                error("????")
            end
            
        elseif valid_sub 
            # if the subject polygon edge is valid proceed along that edge
            # note that in the case where both edges are valid, we choose to prefer the subject polygon edge
            point[1] = 2
            point[2] = nextSubIndex
            return 0, ()

        elseif valid_clip 
            # if the clipping polygon edge is valid proceed along that edge
            point[1] = 1
            point[2] = nextClipIndex
            return 0, ()
        else 
            return -1, (nextClipIndex, nextSubIndex)
        end

    elseif (label == 0 && currentPoint[1] == 2) || (label == -1 && currentPoint[1] == 1) 
        # current point is in exterior vertex of the subject polygon or an interior vertex of the clipping polygon
        # stay on current shape and increment index by 1
        if currentPoint[1] == 1
            c = -1 
        elseif currentPoint[1] == 2 
            c = 1 
        end 
        point[1] = currentPoint[1]
        point[2] = nextIndex(currentPoint[2], wrapPoints[currentPoint[1]],c)

        if startingPoint == point
            return -2, ()
        elseif point[2] > 0 && startingPoint[1] == other12(point[1]) && startingPoint[2] == labels[point[1]][point[2]]
            return -2, ()
        end

        return 0, ()
    else
        # unable to find a valid point
        return -1, ()
    end

    
end

function nextIndex(index, wrapPoints, c)

    # increment index
    indN = index + c*1

    if indN != 0 # if next index is not 0 
        if wrapPoints[indN] > 0 # if next index is a wrapPoint
            if c == 1 # if proceeding clockwise
                return wrapPoints[indN]
            elseif c == -1 # if proceeding counter clockwise 
                while true # iterate clockwise until you find the previous wrap point and jump back to that point
                    index += 1
                    if wrapPoints[index] > 0 
                        return index - 1
                    end
                end
            end
        else
            return indN
        end
    else # next index is 0 (only happens when iterating counter clockwise)
        # iterate clockwise until you find the previous wrap point and jump back to that point
        found = false
        while !found
            index += 1
            if wrapPoints[index] > 0 
                return index-1
            end
        end
    end

end

function nextIndex(index, indices)
    if indices[index + 1] == -1
        for i = index:-1:2
            if indices[i-1] == -1
                return i
            end
        end
        return 1
    else
        return index + 1
    end
end

function previousIndex(index, indices, maxIndex)
    if index == 1
        for i = 2:maxIndex
            if indices[i] == -1
                return i - 1
            end
        end
    else
        return index - 1
    end
end

function findStartingPoint(labels, inOrOutBound, cvn, svn)
    # check for intersections and interior vertices of the clipping polygon
    found = false
    for i = 1:cvn
        if labels[1][i] == -1 # point is an interior vertex on the clipping polygon
            # set the current point to be the starting point
            startingPoint = [1, i]
            found = true
            break 
        elseif labels[1][i] > 0 # point is an intersection 
           
            if inOrOutBound[1][i] == 1
                 # if intersection is inbound
                startingPoint = [1, i]
                found = true
                break
            elseif inOrOutBound[1][i] == -1
                # if intersection is outbound, switch to the subject polygon
                startingPoint = [2, labels[1][i]]
                found = true
                break
            end
        end
    end

    # check for exterior vertices of the subject polygon
    if !found
        for i = 1:svn
            if labels[2][i] == 0 # point is an exterior vertex on the subject polygon
                startingPoint = [2, i]
                found = true
                break 
            elseif labels[2][i] > 0 # point is an intersection     
                if inOrOutBound[2][i] == 1
                    # if intersection is inbound
                   startingPoint = [2, i]
                   found = true
                   break
               elseif inOrOutBound[2][i] == -1
                   # if intersection is outbound, switch to the clipping polygon
                   startingPoint = [1, labels[2][i]]
                   found = true
                   break
               end
            end
        end
    end

    if !found
        startingPoint = Int64[]
    end

    return found, startingPoint
end

function findStartingPoint!(startingPoint, labels, inOrOutBound, wrapPoints, cvn, svn, usedPoints)
    # check for intersections and interior vertices of the clipping polygon
    found = false
    for i = 1:cvn
        if !usedPoints[1][i]
            if labels[1][i] == -1 # point is an interior vertex on the clipping polygon
                # set the current point to be the starting point
                startingPoint[1] = 1
                startingPoint[2] = i
                found = true
                break 
            elseif labels[1][i] > 0 # point is an intersection 

                check = inOrOutBound[1][nextIndex(i,wrapPoints[1],-1)]
                label = labels[1][nextIndex(i,wrapPoints[1],-1)]
            
                if label == -1 || check == 1
                    # if next point on the clipping polygon is an interior vertex or an inbound intersection
                    startingPoint[1] = 1
                    startingPoint[2] = i
                    found = true
                    break
                end
            end
        end
    end

    # check for exterior vertices of the subject polygon
    if !found
        for i = 1:svn
            if !usedPoints[2][i]
                if labels[2][i] == 0 # point is an exterior vertex on the subject polygon
                    startingPoint[1] = 2
                    startingPoint[2] = i
                    found = true
                    break 
                elseif labels[2][i] > 0 # point is an intersection     
                    if inOrOutBound[2][i] == -1
                        # if intersection is outbound
                        startingPoint[1] = 2
                        startingPoint[2] = i
                        found = true
                        break
                    end
                end
            end
        end
    end

    return found
end

function validPoint(point, label)

    if (label == -1 && point[1] == 1) || (label == 0 && point[1] == 2) || label > 0
        return true
    else 
        return false
    end

end

# function nextEdge!(edge, vertices, vertList, index)

    #     if vertList[index + 1] == -1
    #         startPoint = 1
    #         for i = index:-1:1
    #             if vertList[i] == -1
    #                 startPoint = i + 1
    #                 break
    #             end
    #         end
    #     else
    #         startPoint = index + 1
    #     end

    #     if isdefined(edge,2)
    #         edge[1] = edge[2]
    #         edge[2][:] = vertices[vertList[startPoint]]
    #     else 
    #         edge[1] = vertices[vertList[index]]
    #         edge[2] = vertices[vertList[startPoint]]
    #     end
    #     # edge[1] = vertices[vertList[index]]
    #     # edge[2] = vertices[vertList[startPoint]]
    # end


# function nextPoint(currentPoint :: Vector, labels :: Vector{Vector{Int64}}, wrapPoints :: Vector{Vector{Int64}}, startingPoint :: Vector{Int64})

    #     label = labels[currentPoint[1]][currentPoint[2]]

    #     if currentPoint[1] == 1
    #         c = -1 
    #     elseif currentPoint[1] == 2 
    #         c = 1 
    #     end 

    #     if label > 0 # point is an intersection
    #         # swap to other shape an increment index by 1
    #         shape = other12(currentPoint[1])
    #         point = [shape, nextIndex(label, wrapPoints[shape], -c)] :: Vector{Int64}
    #         flag = 0

    #         if point == startingPoint
    #             return point , -2
    #         end

    #         # if next point is not valid
    #         if !validPoint(point, labels[point[1]][point[2]])
    #             # try staying on the same shape instead
    #             point = [currentPoint[1], nextIndex(currentPoint[2], wrapPoints[currentPoint[1]],c)] :: Vector{Int64}

    #             if point == startingPoint
    #                 return point , -2
    #             elseif !validPoint(point, labels[point[1]][point[2]])
    #                 # unable to find a valid point
    #                 return Int64[0,0] , -1
    #             else
    #                 return point, flag
    #             end
    #         end

    #     elseif label == 0 && currentPoint[1] == 2 # current point is in exterior vertex of the subject polygon
    #         # stay on current shape and increment index by 1
    #         point = [currentPoint[1], nextIndex(currentPoint[2], wrapPoints[currentPoint[1]],c)] :: Vector{Int64}
    #         flag = 0

    #     elseif label == -1 && currentPoint[1] == 1 # current point is an interior vertex of the clipping polygon
    #         # stay on current shape and increment index by 1
    #         point = [currentPoint[1], nextIndex(currentPoint[2], wrapPoints[currentPoint[1]],c)] :: Vector{Int64}
    #         flag = 0
    #     else
    #         # unable to find a valid point
    #         return Int64[0,0] , -1
    #     end

    #     if point == startingPoint
    #         return point , -2
    #     end

    #     return point, flag
    # end