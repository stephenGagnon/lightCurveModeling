"""
  Fraction of visible light that strikes a facet and is reflected to the
  observer

 INPUTS ---------------------------------------------------------------

  A -- the attitude matrix (inertial to body)

  geometry -- a structure containg various parameters describing the
  relative possitions and directions of the observer and sun in the
  inertial frame. The comonenets are as follows:

  usun -- vector from rso to sun (inertial)
  uobs -- vector from rso to the jth observer (inertial)
  d -- distance from rso to observer j
  C -- sun power per square meter

  facet -- a structure contining various parameters describing the facet
  being observed

  Area -- facet area
  unb -- surface normal of the ith facet (body frame)
  uub,uvn body -- in plane body vectors completing the right hand rule
  Rdiff,Rspec -- spectral and diffusion parameters of the facet
  nv,nu -- coefficients to determine the in-plane distribution of
  spectral reflection

 OUTPUTS --------------------------------------------------------------

  F -- total reflected intensity (Fobs)

  ptotal -- total reflectance fraction (rho)

 CODE -----------------------------------------------------------------
"""
# wrapper for the light curve model that handles different attitude parameterizations and unwraps the input structures
function Fobs(att :: anyAttitude, scen :: scenario, a=1.0, f=1.0)

    # attitude parameterization handling
    if (typeof(att) <: Vec) & (length(att) == 3)
        rotFunc = ((A,v) -> p2A(A,a,f)*v)
    elseif ((typeof(att) <: Vec) & (length(att) == 4)) | (typeof(att) == quaternion)
        rotFunc = qRotate!
    elseif (typeof(att) <: Mat) & (size(att) == (3,3))
        rotFunc = ((A,v) -> A*v)
    elseif typeof(att) <: Union{DCM,MRP,GRP}
        rotFunc = ((A,v) -> any2A(A).A*v)
    elseif typeof(att) <: Real
        rotFunc = ((tht,v) -> [cos(tht) sin(tht);-sin(tht) cos(tht)]*v)
        return _Fobs_2D(att,obj.nvecs,obj.uvecs,obj.vvecs,obj.Areas,obj.nu,obj.nv,
         obj.Rdiff,obj.Rspec,scen.sunVec,scen.obsVecs,scen.d,scen.C,rotFunc)
    else
        error("Please provide a valid attitude. Attitudes must be represented
        as a single 3x1 or 4x1 float array, a 3x3 float array, or any of the
        custom attitude types defined in the attitueFunctions package.")
    end

    return  _Fobs(att, scen, rotFunc)
end

# handles data input as vectors of vectors
function _Fobs(att :: anyAttitude{T1}, scen :: scenario, rotFunc :: Function) where {T1 <: Real}

    # rotate inertial vectors into the body frame for calculations
    # (usun,uobst) = _toBodyFrame(att,scen,rotFunc)
    rotFunc(att, scen.sun.posUnitVecRel, scen.sun.bodyFrameUnitVec)#scen.sun.bodyFrameUnitVec[:] = 

    for i = 1:scen.observerNo
        rotFunc(att, scen.observers[i].posUnitVecRel, scen.observers[i].bodyFrameUnitVec)#scen.observers[i].bodyFrameUnitVec[:] = 
    end
    
    # intialize intensity array where each element of the array corresponds to each observer
    Ftotal = zeros(T1, scen.observerNo)

    # loop through facets # loop through observers
    for j in eachindex(scen.observers)

        # if the scenario object contains the required data for self shadowing, this function performs the self shadowing calculations and stores the modified facet areas, otherwise this function simply performs the line of sight checks
        flag = selfShadowing!(att, scen, j)
        if flag == -1
            return Ftotal, -1
        end

        # add the contribution to the total intensity from each facet
        for i in eachindex(scen.object.surfaces)   
            if scen.object.surfaces[i].isVisible
                Ftotal[j] += _Fobs_facet(scen.object.surfaces[i], scen.observers[j], scen.sun, scen.halfAngleVector, scen.unitScaling)
            end
        end
    end

    return Ftotal, 0
end

# calcuates reflected intensity contribution to a single observer from a single facet 
function _Fobs_facet(surf :: surface, obs :: observer, sun :: pointSource, uh :: Vector, unitScaling :: Float64)

    usun = sun.bodyFrameUnitVec
    uobs = obs.bodyFrameUnitVec
    un = surf.u_n
    uu = surf.u_u
    uv = surf.u_v

    # precalculate some dot products to save time
    usdun = dot(usun,un)
    uodun = dot(uobs,un)
    usduo = dot(usun,uobs)

    # compute the half angle vector between the sun and observer
    # for i in eachindex(scen.halfAngleVector)
    #     uh[i] = (usun[i] + uobs[i]) / sqrt(2 + 2*usduo)
    # end
    uh .= (usun .+ uobs) ./ sqrt(2 + 2*usduo)

    uhdun = dot(uh,un)

    # diffuse reflection
    pdiff = ((28*surf.Rdiff)/(23*pi))*(1 - surf.Rspec)*(1 - (1 - usdun/2)^5)*(1 - (1 - uodun/2)^5)

    # spectral reflection

    # calculate numerator
    if surf.nu == surf.nv # isotropic case
        pspecnum = (surf.nu + 1)/(8*pi)*(surf.Rspec + (1 - surf.Rspec)*(1 - dot(uh,usun))^5)*(uhdun^(surf.nu))
    elseif abs(uhdun-1) < 1e-15 # (1 - dot(uh,un)) < .00001 #deal with numerical issues when uhdun is close to 1
        pspecnum = sqrt((surf.nu + 1)*(surf.nv + 1))/(8*pi)*(surf.Rspec + (1 - surf.Rspec)*(1 - dot(uh,usun))^5)
    else # general non-isotropic case
        pspecnum = sqrt((surf.nu + 1)*(surf.nv + 1))/(8*pi)*(surf.Rspec + (1 - surf.Rspec)*(1 - dot(uh,usun))^5)*(uhdun^((surf.nu*dot(uh,uu)^2 + surf.nv*dot(uh,uv)^2)/(1 - uhdun^2)))
    end

    # return total intensity 
    return (sun.irradiance * unitScaling)/(obs.distance^2)*(pspecnum/(usdun + uodun - (usdun)*(uodun)) + pdiff)*(usdun)*surf.areaClipped*(uodun)
   
end