"""
  calculates the sensitivity of the reflectance from one facet due to small
  changes in attitude

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

  dFobs -- the sensitivty vector of the refelctance (Fobs) wrt small changes
     in attitude

  dpt -- the sensitifity vector of the brightness (rho) wrt small changes
  in the attitude

  Ft -- the total reflectance of the facet

  Code -----------------------------------------------------------------
"""
function dFobs(att :: Vector, unm :: ArrayOfVecs, uum :: ArrayOfVecs, uvm :: ArrayOfVecs, Area :: Vector, nu :: Vector, nv :: Vector, Rdiff :: Vector, Rspec :: Vector, usun :: Vector, uobst :: ArrayOfVecs, d :: Vector, C :: Float64, dDotFunc :: Function, rotFunc :: Function, parameterization)

    Ftotal = zeros(length(uobst),)
    dFtotal = zeros(length(uobst),length(att))

    (unI,uuI,uvI) = _toInertialFrame(att,unm,uum,uvm,rotFunc,parameterization)


    # total = zeros(length(uobst),)
    # dtotal = zeros(length(uobst),3)

    uh = Array{typeof(usun[1]),1}(undef,length(usun))
    # A = p2A(att)

    for i = 1:eachindex(unm)
        un = unI[i]
        uv = uvI[i]
        uu = uuI[i]

        for j = 1:eachindex(uobst)
            uobs = uobst[j]

            check1 = dot(usun,un) <= 0
            check2 = dot(uobs,un) <= 0
            visFlag = check1 | check2

            if visFlag

            else
                # calculate the half angle vector
                # uh = (usun + uobs)./norm(usun + uobs)


                usduo = dot(usun,uobs)
                for k = 1:length(usun)
                    uh[k] = (usun[k] + uobs[k])/sqrt(2 + 2*usduo)
                end
                # uh[1] = (usun[1] + uobs[1])/sqrt(2 + 2*usduo)
                # uh[2] = (usun[2] + uobs[2])/sqrt(2 + 2*usduo)
                # uh[3] = (usun[3] + uobs[3])/sqrt(2 + 2*usduo)

                usuh = dot(usun,uh)
                # define dot products and their derivatives
                uhun = dot(uh,un)
                duhun = dDotFunc(uh,unm[i],att)

                uhuu = dot(uh,uu)
                duhuu = dDotFunc(uh,uum[i],att)

                uhuv = dot(uh,uv)
                duhuv = dDotFunc(uh,uvm[i],att)

                unus = dot(usun,un)
                dunus = dDotFunc(usun,unm[i],att)

                unuo = dot(uobs,un)
                dunuo = dDotFunc(uobs,unm[i],att)


                # diffuse reflection
                k = ((28*Rdiff[i])/(23*pi))*(1 - Rspec[i])
                pdiff = k*(1 - (1 - unus/2)^5)*(1 - (1 - unuo/2)^5)

                dpdiff = (5/2)*k*( ( ((1-.5*unus)^4)*(1-(1-.5*unuo)^5)*dunus ) +
                ( ((1-.5*unuo)^4)*(1 - (1 - .5*unus)^5)*dunuo ) )

                # spectral reflection

                # calculate numerator and account for the case where the half angle
                # vector lines up with the normal vector

                K = sqrt((nu[i] + 1)*(nv[i] + 1))/(8*pi)
                Fref = Rspec[i] + (1-Rspec[i])*(1-usuh)^5


                if uhun == 1 #(uhun≈1)
                    num = K*Fref
                    dnum = [0;0;0]
                else
                    fac = (nu[i]*uhuu^2 + nv[i]*uhuv^2)/(1 - uhun^2)
                    dfac = ((2*nu[i]*uhuu*duhuu + 2*nv[i]*uhuv*duhuv)*(1-uhun^2) +
                     2*uhun*(nu[i]*uhuu^2 + nv[i]*uhuv^2)*duhun)/(1-uhun^2)^2
                    num = K*Fref*uhun^(fac)
                    dnum = (fac*num/uhun)*duhun + num*log(uhun)*dfac
                end


                den = (unus + unuo - (unus)*(unuo))
                dden = (1-unuo)*dunus + (1-unus)*dunuo

                pspec = num/den
                dpspec = (den*dnum-num*dden)/(den^2)

                ptotal = pspec + pdiff
                dptotal = dpdiff + dpspec

                Fsun = C*ptotal*(unus)
                dFsun = C*(dptotal*unus + ptotal*dunus)

                Fobs_ = Fsun/(d[j]^2)*Area[i]*(unuo)
                dFobs_ = Area[i]/(d[j]^2)*(dFsun*unuo + Fsun*dunuo)

                # total[j] += pdiff
                # dtotal[j,:] += dpdiff

                Ftotal[j] += Fobs_
                dFtotal[j,:] += dFobs_
                if uhun == 1

                end
            end

        end
    end
    return Ftotal, dFtotal #, total, dtotal
end

function dFobs_2D(att, obj::targetObject, scen::spaceScenario)
    return dFobs_2D(att, obj.nvecs, obj.uvecs, obj.vvecs, obj.Areas, obj.nu, obj.nv, obj.Rdiff, obj.Rspec, scen.sunVec, scen.obsVecs, scen.d, scen.C)
end

function dFobs_2D(att :: Float64, unm :: ArrayOfVecs, uum :: ArrayOfVecs, uvm :: ArrayOfVecs, Area :: Vector, nu :: Vector, nv :: Vector, Rdiff :: Vector, Rspec :: Vector, usun :: Vector, uobst :: ArrayOfVecs, d :: Vector, C :: Float64)

    Ftotal = zeros(length(uobst),)
    dFtotal = zeros(length(uobst),)

    (unI,uuI,uvI) = _toInertialFrame(att,unm,uum,uvm,rotate_2D, att2D)


    # total = zeros(length(uobst),)
    # dtotal = zeros(length(uobst),3)

    uh = Array{typeof(usun[1]),1}(undef,length(usun))
    # A = p2A(att)

    for i = 1:eachindex(unm)
        un = unI[i]
        uu = uuI[i]

        for j = 1:eachindex(uobst)
            uobs = uobst[j]

            check1 = dot(usun,un) <= 0
            check2 = dot(uobs,un) <= 0
            visFlag = check1 | check2

            if visFlag

            else
                # calculate the half angle vector
                # uh = (usun + uobs)./norm(usun + uobs)


                usduo = dot(usun,uobs)
                for k = 1:length(usun)
                    uh[k] = (usun[k] + uobs[k])/sqrt(2 + 2*usduo)
                end
                # uh[1] = (usun[1] + uobs[1])/sqrt(2 + 2*usduo)
                # uh[2] = (usun[2] + uobs[2])/sqrt(2 + 2*usduo)
                # uh[3] = (usun[3] + uobs[3])/sqrt(2 + 2*usduo)

                usuh = dot(usun,uh)
                # define dot products and their derivatives
                uhun = dot(uh,un)
                # duhun = dDotFunc(uh,unm[i],att)
                duhun = uh' * dAdtht(att) * unm[i]

                uhuu = dot(uh,uu)
                # duhuu = dDotFunc(uh,uum[i],att)
                duhuu = uh' * dAdtht(att) * uum[i]

                unus = dot(usun,un)
                # dunus = dDotFunc(usun,unm[i],att)
                dunus = usun' * dAdtht(att) * unm[i]

                unuo = dot(uobs,un)
                # dunuo = dDotFunc(uobs,unm[i],att)
                dunuo = uobs' * dAdtht(att) * unm[i]


                # diffuse reflection
                k = ((28*Rdiff[i])/(23*pi))*(1 - Rspec[i])
                pdiff = k*(1 - (1 - unus/2)^5)*(1 - (1 - unuo/2)^5)

                dpdiff = (5/2)*k*( ( ((1-.5*unus)^4)*(1-(1-.5*unuo)^5)*dunus ) +
                ( ((1-.5*unuo)^4)*(1 - (1 - .5*unus)^5)*dunuo ) )

                # spectral reflection

                # calculate numerator and account for the case where the half angle
                # vector lines up with the normal vector

                K = sqrt((nu[i] + 1)*(nv[i] + 1))/(8*pi)
                Fref = Rspec[i] + (1-Rspec[i])*(1-usuh)^5


                if uhun == 1 #(uhun≈1)
                    num = K*Fref
                    dnum = 0
                else
                    fac = (nu[i]*uhuu^2)/(1 - uhun^2)
                    dfac = ((2*nu[i]*uhuu*duhuu)*(1-uhun^2) +
                     2*uhun*(nu[i]*uhuu^2)*duhun)/(1-uhun^2)^2
                    num = K*Fref*uhun^(fac)
                    dnum = (fac*num/uhun)*duhun + num*log(uhun)*dfac
                end


                den = (unus + unuo - (unus)*(unuo))
                dden = (1-unuo)*dunus + (1-unus)*dunuo

                pspec = num/den
                dpspec = (den*dnum-num*dden)/(den^2)

                ptotal = pspec + pdiff
                dptotal = dpdiff + dpspec

                Fsun = C*ptotal*(unus)
                dFsun = C*(dptotal*unus + ptotal*dunus)

                Fobs_ = Fsun/(d[j]^2)*Area[i]*(unuo)
                dFobs_ = Area[i]/(d[j]^2)*(dFsun*unuo + Fsun*dunuo)

                # total[j] += pdiff
                # dtotal[j,:] += dpdiff

                Ftotal[j] += Fobs_
                dFtotal[j] += dFobs_
                if uhun == 1

                end
            end

        end
    end
    return Ftotal, dFtotal #, total, dtotal
end