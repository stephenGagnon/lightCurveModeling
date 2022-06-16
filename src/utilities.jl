function simpleSatellite(;vectorized = false)


    ## satellite bus
    # Bus measures 1.75 x 1.7 x 1.8 m.  Difficult to discern which dimension
    # corresponds to which direction (documentation for GEOStar-2 bus does not
    # state this), but my best guess is that the side with the solar panels
    # measures 1.8 x 1.75, while the earth-facing side measures 1.75 x 1.7.
    # Based on coordinate system used for solar panel, then s1 = 1.8, s2= 1.7, s3=1.75.
    s1 = 1.8
    s2 = 1.7
    s3 = 1.75
    l1 = s1/2
    l2 = s2/2
    l3 = s3/2

    # points corresponding to the verticies of the bus
    p_bus = [l1  l2  l3; # front top right
        l1 -l2  l3; # front top left
        l1 -l2 -l3; # front bottom left
        l1  l2 -l3; # front bottom right
        -l1  l2  l3; # back top right
        -l1  l2 -l3; # back bottom right
        -l1 -l2 -l3; # back bottom left
        -l1 -l2  l3] # back top left
    npb = size(p_bus,1)

    # the sets of verticies corresponding to each facet
    K_bus = [[1 2 3 4], # front panel
        [5 6 7 8], # back panel
        [4 3 7 6], # bottom panel
        [1 5 8 2], # top panel
        [1 4 6 5], # right panel
        [2 8 7 3]] # left panel

    # bus panel areas
    Area_bus = [s3*s2 s3*s2 s1*s2 s1*s2 s3*s1 s3*s1]

    # moment of inrtia of bus about its COM
    m_bus = 1792                 # kg
    J_bus = (m_bus/12)*diagm([s2^2 + s3^2, s1^2 + s3^2, s1^2 + s2^2])

    ## satellite solar panel
    SPw = 1.6
    SPl = 4
    SP_off = l2 + SPl/2
    SP_c = [0;SP_off;0]
    offset = .01

    p_panel1 = [offset l2      -SPw/2;
        offset l2 + SPl -SPw/2;
        offset l2 + SPl  SPw/2;
        offset l2       SPw/2
        -offset l2      -SPw/2;
        -offset l2 + SPl -SPw/2;
        -offset l2 + SPl  SPw/2;
        -offset l2       SPw/2]

    p_panel2 = copy(p_panel1)
    p_panel2[:,1:2] = -p_panel2[:,1:2]
    p_panel = [p_panel1; p_panel2]

    # moment of inertia of SP about its COM
    m_SP = 50  # kg
    J_SP = (m_SP/2)/12*diagm([(SPl^2 + SPw^2), SPw^2, SPl^2])

    # Solar Panel angle offset, as measured from +X axis to normal vector
    theta = -25*pi/180

    # Solar Panel rotates about Y-axis
    R = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]

    J_SP = R'*J_SP*R

    p_panel = (R*p_panel')'

    K_panel = [[1 2 3 4] .+ npb, # front right panel
        [8 7 6 5] .+ npb, # back right panel
        [9 10 11 12] .+ npb, # front left panel
        [16 15 14 13] .+ npb] # back left panel


    npbp = npb + size(p_panel,1)

    # add solar panel areas
    Area = [Area_bus SPw*SPl*ones(1,4)]

    ## satellite antenae
    # dish radius
    d_r = 1.872/2
    # dish offset from top of bus
    d_off = l3 + d_r
    # coordinates of center of dish
    d_c = [0 0 d_off]'

    # area of dish
    Area = [Area (pi*d_r^2)*ones(1,2)]

    # generate points around the dish
    tht = 0:pi/40:2*pi
    p_dish = [zeros(length(tht),1).+offset   sin.(tht)   cos.(tht);
                zeros(length(tht),1).-offset sin.(tht)   cos.(tht)]

    for i = 1:size(p_dish,1)
        p_dish[i,:] += d_c
    end

    temp = [npbp .+ (1:length(tht));]
    K_dish = [temp, temp[end:-1:1]]

    # moment of inertia of Dish about its COM
    m_dish = 50  # kg
    J_dish = m_dish*d_r^2*diagm([0.5, 0.25, 0.25])

    ## body frame vectors
    P = [p_bus;p_panel;p_dish]
    vertices = P
    K = [K_bus; K_panel; K_dish]
    facetVerticesList = K
    facetNo = length(K)

    nvecs = zeros(3,length(K)-2)
    uvecs = zeros(3,length(K)-2)
    vvecs = zeros(3,length(K)-2)

    for i = 1:facetNo-2
        vec1 = P[K[i][2],:]-P[K[i][1],:]
        vec2 = P[K[i][3],:]-P[K[i][2],:]
        nvecs[:,i] = cross(vec1,vec2)./norm(cross(vec1,vec2))
        vvecs[:,i] = vec1./norm(vec1)
        uvecs[:,i] = cross(nvecs[:,i],vvecs[:,i])
    end

    # store body vectors
    nvecs = [nvecs [1 0 0]' [-1 0 0]']
    uvecs = [uvecs [0 1 0]' [0 -1 0]']
    vvecs = [vvecs [0 0 1]' [0 0 1]']

    bodyFrame = Matrix(1.0I,3,3)

    # in plane parameters
    nu = 1000*ones(1,facetNo)
    nv = 1000*ones(1,facetNo)

    # spectral and diffusion parameters
    Rdiff = [.6*ones(1,6) .05*ones(1,4)  .6 .6] # bus, solar panel, dish #.6*ones(1,2) .05 .26*ones(1,2) .04
    Rspec = [.26*ones(1,6) .04*ones(1,4) .275 .26]

    ## moment of inertia calcualtions

    # find COM
    # solar panels cancel and main bus is centered at origin
    COM = m_dish/(m_dish + m_bus + 2*m_SP)*d_off

    # find moment of inertia about bus center
    J_SP_bus = J_SP + (m_SP/2).*((SP_c'*SP_c).*Matrix(1.0I,3,3) - SP_c*SP_c')
    J_dish_bus = J_dish + m_dish.*((d_c'*d_c).*Matrix(1.0I,3,3) - d_c*d_c')

    J_tot = J_bus + 2*J_SP_bus + J_dish_bus

    # moment of Intertia about the COM
    J = J_tot  - (m_dish + m_bus + 2*m_SP).*((COM'*COM).*Matrix(1.0I,3,3) .- COM*COM')

    fullStruct = targetObjectFull(facetNo,vertices,facetVerticesList,Area,nvecs,
    vvecs,uvecs,nu,nv,Rdiff,Rspec,J,bodyFrame)

    if !vectorized
        Area = Area[:]
        nu = nu[:]
        nv = nv[:]
        Rdiff = Rdiff[:]
        Rspec = Rspec[:]
        nvecstemp = nvecs
        nvecs = Array{Array{Float64,1},1}(undef,size(nvecstemp,2))
        uvecstemp = uvecs
        uvecs = Array{Array{Float64,1},1}(undef,size(nvecstemp,2))
        vvecstemp = vvecs
        vvecs = Array{Array{Float64,1},1}(undef,size(nvecstemp,2))

        for i = 1:facetNo
            nvecs[i] = nvecstemp[:,i]
            uvecs[i] = uvecstemp[:,i]
            vvecs[i] = vvecstemp[:,i]
        end
    end
    simpleStruct = targetObject(facetNo,Area,nvecs,vvecs,uvecs,nu,nv,Rdiff,Rspec,J,bodyFrame)

    return simpleStruct, fullStruct
end

function customSatellite(objParams; vectorized = false)

    (satSimple, satFullSimple) = simpleSatellite(vectorized = vectorized )

    p = fieldnames(targetObjectFull)
    objvars = Array{Any,1}(undef,length(p))

    for i = 1:length(p)
        if haskey(objParams,p[i])
            objvars[i] = objParams[p[i]]
        else
            objvars[i] = getproperty(satFullSimple,p[i])
        end
    end

    sat = spaceScenario(objvars[1],objvars[4:end]...)
    satFull = spaceScenario(objvars...)

    return sat, satFull
end

function simpleScenario(;vectorized = false)

    # C -- sun power per square meter
    C = 455.0 #W/m^2

    # number of observers
    obsNo = 4

    # distance from observer to RSO
    # obsDist = 35000*1000*ones(1,obsNo)
    obsDist = 1*1000*ones(1,obsNo)

    #body vectors from rso to observer (inertial)
    r = sqrt(2)/2
    v = sqrt(3)/3
    obsVecs = [1  r  v  v  r  0 -r  0 -r
               0  r -v  v  0  r -r  r  0
               0  0  v  v  r  r  0 -r  r]
    # obsVecs = [v  r  0 -r  0 -r
    #            v  0  r -r  r  0
    #            v  r  r  0 -r  r]
    obsVecs = obsVecs[:,1:obsNo]

    # usun -- vector from rso to sun (inertial)
    sunVec = [1.0; 0; 0]

    if !vectorized
        obsvectemp = obsVecs
        obsVecs = Array{Array{Float64,1},1}(undef,size(obsvectemp,2))

        for i = 1:obsNo
            obsVecs[i] = obsvectemp[:,i]
        end

        obsDist = obsDist[:]
    end
    return spaceScenario(obsNo,C,obsDist,sunVec,obsVecs)
end

function customScenario(scenParams; vectorized = false)

    scenarioSimple = simpleScenario(vectorized  = vectorized)

    p = fieldnames(spaceScenario)
    scenvars = Array{Any,1}(undef,length(p))

    for i = 1:length(p)
        if haskey(scenParams,p[i])
            scenvars[i] = scenParams[p[i]]
        else
            scenvars[i] = getproperty(scenarioSimple,p[i])
        end
    end

    # set observer number to match the number of observer vectors
    if typeof(scenvars[5]) == Vector{Vector{Float64}}
        scenvars[1] = length(scenvars[5])
        dtemp = zeros(length(scenvars[5]))
        if length(scenvars[3]) < length(scenvars[5])
            dtemp[1:length(scenvars[3])] = scenvars[3]
            dtemp[length(scenvars[3])+1:end] .= scenvars[3][1]
        else
            dtemp[1:length(scenvars[5])] = scenvars[3][1:length(scenvars[5])]
        end
        scenvars[3] = dtemp
    elseif typeof(scenvars[5]) == Matrix{Float64}
        scenvars[1] = size(scenvars[5],2)
        dtemp = zeros(size(scenvars[5],2))
        dtemp[1:length(scenvars[3])] = scenvars[3]
        dtemp[length(scenvars[3])+1:end] .= scenvars[3][1]
        scenvars[3] = dtemp
    end

    scenario = spaceScenario(scenvars...)

    return scenario
end

function simpleScenarioGenerator(;vectorized = false)
    (sat, satFull) = simpleSatellite(vectorized = vectorized )
    scenario = simpleScenario(vectorized  = vectorized)
    return sat, satFull, scenario
end

function objectGenerator2D()

    obsvecs = Array{Array{Float64,1},1}(undef,2)
    obsvecs[1] = [1;1]./norm([1;1])
    obsvecs[2] = [1;-1]./norm([1;-1])

    facetNo = 3

    # side lengths
    L = [3,4,5]
    # in 2D length is the analogue of area
    Area = L

    p1 = [L[1] 0]
    p2 = [0 0]
    p3 = [0 -L[2]]

    verts = [p1;p2;p3]
    vertList = [[2,1],[3,2],[1,3]]

    # in-plane geometry
    # side normals
    nvecs = Array{Array{Float64,1},1}(undef,facetNo)
    uvecs = Array{Array{Float64,1},1}(undef,facetNo)
    # create empty unused array for vvecs
    vvecs = Array{Array{Float64,1},1}(undef,0)

    for i = 1:facetNo
        temp = verts[vertList[i][2],:] - verts[vertList[i][1],:]
        temp = temp./norm(temp)
        uvecs[i] = temp
        temp = reverse(temp)
        temp[1] = - temp[1]
        nvecs[i] = temp
    end

    # in plane parameters
    # nu = [1000 1000 1000];
    # nv = [1000 1000 1000];
    nu = 1000*ones(1,facetNo)
    nv = 1000*ones(1,facetNo)

    # spectral and diffusion parameters
    Rdiff = .5*ones(1,facetNO);
    # Rdiff = [.1 .2 .3];
    # Rspec = [.4 .5 .6];
    Rspec = .5*ones(1,facetNo);

    bodyFrame = [1 0; 0 1]

    Jx = (L[1]*L[2]^3)/12
    Jy = (L[2]*L[1]^3)/12
    J = [Jx 0 ; 0 Jy]

    return targetObject(facetNo, Areas, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame), targetObjectFull(facetNo, verts, vertList, Area, nvecs, vvecs, uvecs, nu, nv, Rdiff, Rspec, J, bodyFrame)
end

function scenarioGenerator2D()
    # C -- sun power per square meter
    C = 455.0 #W/m^2

    # number of observers
    obsNo = 2

    # distance from observer to RSO
    # obsDist = 35000*1000*ones(1,obsNo)
    obsDist = 35000*1000*ones(1,obsNo)

    #body vectors from rso to observer (inertial)
    r = sqrt(2)/2
    obsVecs = [1  r  r  0 -r -r -1  0;
               0  r -r  1 -r  r  0 -1]

    obsVecs = obsVecs[:,1:obsNo]

    # usun -- vector from rso to sun (inertial)
    sunVec = [1.0; 0]

    if !vectorized
        obsvectemp = obsVecs
        obsVecs = Array{Array{Float64,1},1}(undef,size(obsvectemp,2))

        for i = 1:obsNo
            obsVecs[i] = obsvectemp[:,i]
        end

        obsDist = obsDist[:]
    end
    return spaceScenario(obsNo,C,obsDist,sunVec,obsVecs)
end

function customScenarioGenerator(;scenParams = nothing, objParams = nothing, vectorized = false)

    (satSimple, satFullSimple) = simpleSatellite(vectorized = vectorized )

    if ~isnothing(scenParams)
        scenario = customScenario(scenParams, vectorized = vectorized)
    else
        scenario = simpleScenario()
    end

    if ~isnothing(objParams)
        sat, satFull = customSatellite(objParams, vectorized = vectorized)
    else
        sat, satFull = simpleSatellite()
    end

    return sat, satFull, scenario
end

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

    FIM = dhdx'*inv(R)*dhdx
end
