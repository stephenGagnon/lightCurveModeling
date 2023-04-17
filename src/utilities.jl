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

function dot3(a :: Vec, b :: Vec)
    
    if length(a) < 3 || length(b) < 3
        throw(DimensionMismatch("vectors ahve length less than 3"))
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

