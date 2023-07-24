function plotShape2D(shape, vertices, color)
    plotlyjs()
    p = plot()
    plottingShape = []

    for i = 1:length(shape)
        if shape[i] != -1
            push!(plottingShape, shape[i])
        else
            for j = i-1:-1:1
                if shape[j] == -1 
                    push!(plottingShape, shape[j+1])
                    break
                elseif j == 1
                    push!(plottingShape, shape[j])
                end
            end
            plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color)
            plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color, seriestype=:scatter)
            plottingShape = []
        end
    end

    Base.invokelatest(display, p)
    return p
end

function plotShape2D!(shape, vertices, color)
    plotlyjs()
    p = plot!()
    plottingShape = []
    plotFlag = false
    for i = 1:length(shape)
        if shape[i] != -1
            push!(plottingShape, shape[i])
            plotFlag = false
        else
            for j = i-1:-1:1
                if shape[j] == -1 
                    push!(plottingShape, shape[j+1])
                    break
                elseif j == 1
                    push!(plottingShape, shape[j])
                    break
                end
            end
            plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color)
            plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color, seriestype=:scatter)
            plottingShape = []
            plotFlag = true
        end
    end

    if !plotFlag
        plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color)
        plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), color = color, seriestype=:scatter)
    end

    Base.invokelatest(display, p)
    return p
end

function plotFacets(clippingPoly,subjectPoly,vertices)
    plotShape2D(clippingPoly.tempVertexIndices[1:clippingPoly.tempVertexNo],vertices,:green)
    plotShape2D!(subjectPoly.tempVertexIndices[1:subjectPoly.tempVertexNo],vertices,:red)
end

function plotEdge2D(edge)
    plotlyjs()
    p = plot(getindex.(edge,1),getindex.(edge,2))
    plotPoint2D!(edge[1])
    Base.invokelatest(display, p)
    return p
end

function plotEdge2D!(edge)
    plotlyjs()
    p = plot!(getindex.(edge,1),getindex.(edge,2))
    plotPoint2D!(edge[1])
    Base.invokelatest(display, p)
    return p
end

function plotEdges(edge1,edge2)
    plotEdge2D(edge1)
    plotEdge2D!(edge2)
end

function plotPoint2D!(point)
    plotlyjs()
    p = scatter!([point[1]],[point[2]])
    Base.invokelatest(display, p)
    return p
end

function plotShape3D(shape, vertices, color = nothing)
    plotlyjs()
    p = plot()
    plottingShape = []
    for i = 1:length(shape)
        if shape[i] != -1
            push!(plottingShape, shape[i])
        else
            for j = i-1:-1:1
                if shape[j] == -1 
                    push!(plottingShape, shape[j+1])
                    break
                elseif j == 1
                    push!(plottingShape, shape[j])
                    break
                end
            end
            if isnothing(color)
                p = plot(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), getindex.(vertices[plottingShape],3))
            else
                p = plot(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), getindex.(vertices[plottingShape],3), color = color)
            end
            Base.invokelatest(display, p)
            plottingShape = []
        end
    end
    return p
end

function plotShape3D!(shape, vertices, color = nothing)
    plotlyjs()
    p = plot!()
    plottingShape = []
    if !any(shape .== -1)
        error()
    end
    for i = 1:length(shape)
        if shape[i] != -1
            push!(plottingShape, shape[i])
        else
            for j = i-1:-1:1
                if shape[j] == -1 
                    push!(plottingShape, shape[j+1])
                    break
                elseif j == 1
                    push!(plottingShape, shape[j])
                end
            end
            if isnothing(color)
                p = plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), getindex.(vertices[plottingShape],3))
            else
                p = plot!(getindex.(vertices[plottingShape],1), getindex.(vertices[plottingShape],2), getindex.(vertices[plottingShape],3), color = color)
            end
            Base.invokelatest(display, p)
            plottingShape = []
        end
    end
    return p
end

function plotSat(scen)

    plotlyjs()
    verts = scen.object.projectionData.vertices

    s = scen.object.surfaces[1]
    facet = s.vertexIndices
    p = plot()
    plotShape3D(facet,verts)
    # Base.invokelatest(display, p)
    
    for i = 2:length(scen.object.surfaces)
        s = scen.object.surfaces[i]
        facet = s.vertexIndices
        # p = plot!(getindex.(verts[facet],1),getindex.(verts[facet],2),getindex.(verts[facet],3))
        p = plotShape3D!(facet,verts)
        # Base.invokelatest(display, p)
        # display(plot!(getindex.(verts[facet],1),getindex.(verts[facet],2),getindex.(verts[facet],3), xlim = (-5,5), ylim = (-5,5), zlim = (-5,5), xlabel = "x", ylabel = "y", zlabel = "z"))
    end

    xlims!(-5,5)
    ylims!(-5,5)
    zlims!(-5,5)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")

    return p

end

function plotVisibleFacets(scen)

    plotlyjs()
    verts = scen.object.projectionData.vertices

    # if scen.object.surfaces[1].isVisible
    #     s = scen.object.surfaces[1]
    #     facet = s.vertexIndices
    #     # p = plot(getindex.(verts[facet],1),getindex.(verts[facet],2),getindex.(verts[facet],3))
    #     p = plotShape3D(facet,verts)
    #     # Base.invokelatest(display, p)
    # else 
    p = plot()
    # end
    
    for s in scen.object.surfaces #i = 2:length(scen.object.surfaces)
        if s.isVisible
            facet = s.vertexIndices
            p = plotShape3D!(facet,verts)
        end
    end

    xlims!(-5,5)
    ylims!(-5,5)
    zlims!(-5,5)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")

end

function plotClippedSat(scen)
    plotlyjs()
    verts = scen.object.projectionData.vertices

    # if scen.object.surfaces[1].isVisible
    #     s = scen.object.surfaces[1]
    #     facet = s.clippedGeometry.vertexIndices
    #     # p = plot(getindex.(verts[facet],1),getindex.(verts[facet],2),getindex.(verts[facet],3))
    #     p = plotShape3D(facet,verts)
    #     # Base.invokelatest(display, p)
    # else 
    p = plot()
    # end

    for s in scen.object.surfaces #i = 2:length(scen.object.surfaces)
        if s.isVisible
            # s = scen.object.surfaces[i]
            facet = s.clippedGeometry.vertexIndices[1:s.clippedGeometry.vertexNo]
            p = plotShape3D!(facet,verts)
        end
    end

    xlims!(-5,5)
    ylims!(-5,5)
    zlims!(-5,5)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")
end