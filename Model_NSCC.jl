using Interpolations, QuadGK, Roots, DifferentialEquations

#Function to find (closest) position in array of certain value
#!!It works for increasing values on column!!
function findpos(matrix,value,column)
    for b = 1:length(matrix[:,column])
        Δx1 = abs(matrix[b,column] - value)
        Δx2 = abs(matrix[b+1,column] - value)
        if matrix[b+1,column] >= value
            if Δx1 < Δx2
                return b
                break
            else
                return b+1
                break
            end
        end
    end
end

#Function to find (closest) position in array of certain value
#!!It works for decreasing values on column!!
function findposDec(matrix,value,column)
    for j = 1:length(matrix[:,column])
        Δx1 = abs(matrix[j,column] - value)
        Δx2 = abs(matrix[j+1,column] - value)
        if matrix[j+1,column] < value
            if Δx1 < Δx2
                return j
                break
            else
                return j+1
                break
            end
        end
    end
end

#Evaluating derivative functions to find quantities at threshold
function derThr(eos, k)
    s = 0.
    for j = 1:k
        Δp = eos[j+1,1] - eos[j,1]
        Δϵ = eos[j+1,2] - eos[j,2]
        s = Δp/Δϵ
    end
    return s
end

#creating an interpolated function for speed of sound modification
function csInterp(eos, k, p1, matr_ρ, matr_cs, ds)
    #define matrix to store (ρ,cs) for segments modification
    segm_matrix = Matrix{Float64}(undef,0,2) 
    #segm_matrix = vcat(segm_matrix, [eos[k,3] p1])
    segm_matrix = vcat(segm_matrix, [eos[k,3] sqrt(p1)])
    for j = 1:7
      segm_matrix = vcat(segm_matrix, [matr_ρ[ds,j] matr_cs[ds,j]])  
    end
    return linear_interpolation(segm_matrix[:,1],segm_matrix[:,2])
end

#create high energy EOS
function heEOS(he_eos::AbstractMatrix{T}, eos, k, cs_interp) where T
    Δρ = 10.0^-5
    M = div((10.0^-3 - eos[k,3]),10.0^-5)+div((10.0^-2 - 10.0^-3),10.0^-4)+
        div((ρfin - 10.0^-2),10.0^-3)    
    ρold = eos[k,3]
    ϵold = eos[k,2]
    pold = eos[k,1]
    
    for j = 1:M-1
        ρnew = ρold + Δρ
        ϵnew = ϵold + Δρ * (ϵold + pold)/ρold
        pnew = pold + cs_interp(ρold)^2 * Δρ * (ϵold + pold)/ρold #cs^2=dp/dϵ
        if ρnew > ρfin #break cycle if we reach last mass density
            break
        end
        he_eos = vcat(he_eos, [pnew ϵnew ρnew])
        #store new data in old variable for next loop
        ρold = ρnew 
        ϵold = ϵnew
        pold = pnew
        if ρnew < 10.0^-3
            Δρ = 10.0^-5
        elseif ρnew < 10.0^-2
            Δρ = 10.0^-4
        elseif ρnew < 10^-1
            Δρ = 10.0^-3
        else 
            Δρ = 10.0^-2
        end
    end
    return he_eos
end

#merge two EOSs (low density and QCD phase)
function mergeEOS(eos_matrix::AbstractMatrix{T}, he_eos, eos, k) where T
    N = k + length(he_eos[:,1])
    for j = 1:N
        if j <= k
            eos_matrix = vcat(eos_matrix, [eos[j,1] eos[j,2] eos[j,3]])
        else
            eos_matrix = vcat(eos_matrix, [he_eos[j-k,1] he_eos[j-k,2] he_eos[j-k,3]])
        end 
    end
    return eos_matrix
end

#BUILDING EOS before Λ transition
function build(eos_matrix::AbstractMatrix{T}, eos, matr_ρ, matr_cs, ds) where T
    #find quantities at transition threshold
    k = findpos(eos,ρt,3)
    
    #Evaluating derivative functions to find quantities at threshold
    p1 = derThr(eos, k)
    
    #creating an interpolated function for speed of sound modification
    cs_interp = csInterp(eos, k, p1, matr_ρ, matr_cs, ds)
    
    #define matrix to store new part of EOS
    he_eos = Matrix{Float64}(undef,0,3)
    he_eos = heEOS(he_eos, eos, k, cs_interp)
        
    #merge the two EOSs at right threshold
    eos_matrix = mergeEOS(eos_matrix, he_eos, eos, k)
    return eos_matrix
end

#define interpolated EOS before and after transition
function eos_interp(x) 
    if x < pc
        ϵ_fluid(x)
    elseif x >= pc
        ϵ_fluid(x + Λ) + Λ
    end
end

#define derivative of interpolated EOS before and after transition
function eos_prime_interp(x) 
    if x < pc
        ϵ_prime(x)
    elseif x >= pc
        ϵ_prime(x + Λ)
    end
end

#Computation of mass density, it is needed for moment of inertia calculation
#dln(ρ)=dϵ/(p+ϵ), everything is total=fluid+Λ contribution
function intF(x)
    1/(x - Λ + p_fluid(x - Λ) )
end

#function to define the new ρ after transition
function completeEos(ρ_matr::AbstractMatrix{T}, matr) where T
    #find value of ϵ and ρ just at threshold (outside)
    f(x) = p_fluid(x) - pc
    ϵ_out = find_zero(f, (matr[1,2],matr[end,2]) )
    ρ_out = ρ_fluid(ϵ_out)
    
    #find quantities after the transition at the threshold
    fIn(x) = p_fluid(x) - Λ - pc
    ϵ_fIn = find_zero(fIn, (matr[1,2],matr[end,2]) ) #ϵ_fluid inside at junction
    ϵ_tIn = ϵ_fIn + Λ #ϵ_tot inside at junction
    #we are assuming that ρtot is not only fluid inside, but has Λ contribution
    ρ_In = ρ_out * ((ϵ_tIn + pc)/(ϵ_out + pc))
    
    k = findpos(matr,ϵ_fIn,2) #find closest value to ϵfIn in matrix
     
    ρ_matr = vcat(ρ_matr, [ϵ_tIn ρ_In] )#first element at threshold
    for i=1:length(matr[:,1])-k
        integral, err = quadgk(intF, ϵ_tIn, matr[i+k,2] + Λ)#value of integral and its error
        ρtemp = ρ_In * exp( integral ) 
        ρ_matr = vcat(ρ_matr, [matr[i+k,2]+Λ ρtemp] )
    end
    return ρ_matr
end;

#define interpolated ρ before and after transition
function ρ_interp(x)
    if x < pc
        ρ_fluid(eos_interp(x) )
    elseif x >= pc
        ρ_inside(eos_interp(x) )
    end
end;

#Compute moment of inertia #REMOVE THIS??
function InCalc(sol, scaling_function, Rs)
    #remove knots
    r_array = Interpolations.deduplicate_knots!(sol.t;move_knots = false)
    pressure = linear_interpolation(r_array,sol[3,:])
    #h_func = linear_interpolation(sol.t,sol[2,:])
    #gtt(x) = scaling_function * h_func(x)
    grrInv = linear_interpolation(sol.t,sol[1,:])#This is f, but f=grr^-1
    densityI(x) = ρ_interp(pressure(x)) * x^4 * sqrt( 1/grrInv(x) )#* sqrt( gtt(x)/grrInv(x) )
    
    integral, err = quadgk(densityI, r0, Rs)
    return 4.0 * π * integral
end

function integrator(P0)
    #condition when defining surface and stop integration
    condition(u, r, integrator) =  u[3]/P0 < 10.0^-12 #10^-10
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    
    u0 = [f0, h0, P0, H0, β0]
    
    prob = ODEProblem(tov!, u0, rspan)
    
    sol = solve(prob, RK4(), callback = cb, dt = 0.0005, adaptive = false)
    return sol
end;

#Define the system of TOV and tidal deformability eqs to solve
#Note that ϵ'[r]=dϵ/dp*p'=dϵ/dp*du[3]
function tov!(du, u, p, r) #u=f,h,P,H,β with f=grr^-1, h=gtt
    du[1] = (1.0 - u[1] - 8.0 * π * r^2 * eos_interp(u[3]) )/r
    du[2] = -(u[2] * (-1.0 + u[1] - 8.0 * π * r^2 * u[3]))/(r * u[1])
    du[3] = (-1.0 + u[1] - 8.0 * π * r^2 * u[3]) * (u[3] + 
        eos_interp(u[3]) ) /(2.0 * r * u[1])
    du[4] = u[5]
    du[5] = (u[4] * (-u[1]^3 + (1.0 +  8.0 * π * r^2 * u[3])^3 -
            u[1] * (1.0 + 8.0 * π * r^2 * u[3]) * (-3.0 + 60.0 * π * r^2 * u[3] + 
            20.0 * π * r^2 * eos_interp(u[3]) ) 
            + u[1]^2 * (-3.0 + 60.0 * π * r^2 * u[3] +
            8.0 * π * r^3 * eos_prime_interp(u[3]) * (-1.0 + u[1] - 
            8.0 * π * r^2 * u[3]) * (u[3] + eos_interp(u[3]) ) /(2.0 * r * u[1])+
            20.0 * π * r^2 * eos_interp(u[3]) )) + 
            r * u[1] * (-1 + u[1] - 8.0 * π * r^2 * u[3]) * (1.0 + u[1] + 
            4.0 * π * r^2 * u[3] -
            4.0 * π * r^2 * eos_interp(u[3]) ) * u[5])/(r^2 * u[1]^2 * (1.0 - 
            u[1] + 8.0 * π * r^2 * u[3]))
end;

#Solve a cycle of TOV to retrieve a M-R curve (and other relations)
function cycleTOV(DATA_matrix::AbstractMatrix{T}, P0, Pf) where T
    #determine how long the loop is, based on Pf
    if Pf > 10.0^-2
        N = div((10.0^-4 - P0),(2.5*10.0^-6)) + div((10.0^-3 - 10.0^-4),(2.5*10.0^-5)) +
        div((10.0^-2 - 10.0^-3),(2.5*10.0^-4))
    elseif Pf > 10.0^-3
        N = div((10.0^-4 - P0),(2.5*10.0^-6)) + div((10.0^-3 - 10.0^-4),(2.5*10.0^-5)) + 
        div((Pf - 10.0^-3),(2.5*10.0^-4))
    elseif Pf > 10.0^-4
        N = div((10.0^-4 - P0),(2.5*10.0^-6)) + div((Pf - 10.0^-4),(2.5*10.0^-5))
    elseif Pf > 10.0^-5
        N = div((Pf - P0),(2.5*10.0^-6))
    end

    for i = 1:N
        sol = integrator(P0)
    
        Rs = sol.t[end] #radius
        M = Rs/2 * (1 - sol[1,end]) #mass
        scaling_function = (1 - 2 * M / Rs) / sol[2,end] #rescale for h
        y = Rs * sol[5,end]/sol[4,end]
        C = M / (Rs / Lu * 10^-5) #compactness
        #Love number k2
        k2 = (8.0 * (1.0 - 2.0 * C)^2 * C^5 * (2.0 + 2.0 * C * (-1.0 + y) -
            y))/(5.0 * (2.0 * C * (6.0 + C^2 * (26.0 - 22.0 * y) - 3.0 * y +
            4 * C^4 * (1.0 + y) + 3.0 * C * (-8.0 + 5.0 * y) + C^3 * (-4.0 + 6.0 * y)) +
            3.0 * (1.0 - 2.0 * C)^2 * (2.0 + 2.0 * C * (-1.0 + y) - y) * log(1 - 2*C)))
        λ = 2.0/3.0 * k2 * (Rs / Lu * 10^-5)^5 #tidal deformability
        
        #compute inertia if required
        if Inertia == 1
            I = InCalc(sol, scaling_function, Rs)
        else
            I = 0
        end
        
        #store values in matrix, create a new row each M-R loop
        DATA_matrix = vcat(DATA_matrix, [P0/ρu M Rs/Lu*10.0^-5 C k2 λ I/M^3 λ/M^5])
    
        #increment on P0
        if P0 < 10.0^-4       
            P0 = P0 + 2.5*10.0^-6
        elseif P0 < 10.0^-3
            P0 = P0 + 2.5*10.0^-5
        elseif P0 < 10.0^-2
            P0 = P0 + 2.5*10.0^-4
        else 
            P0 = P0 + 2.5*10.0^-3
        end
    end
    return DATA_matrix
end

#Computation of radius and tidal deformability and fixed mass
function λRCalc(DATA_matrix, iMax)
    #remove knots
    M_array = Interpolations.deduplicate_knots!(DATA_matrix[1:iMax,2];move_knots = true)
    
    #interpolate R-M curve and find corresponding radius for Mtest
    RMcurve = linear_interpolation(M_array,DATA_matrix[1:iMax,3])
    Rt = RMcurve(Mtest)
    
    #interpolate λ-M curve and find corresponding λ for Mtest
    λMcurve = linear_interpolation(M_array,DATA_matrix[1:iMax,6])
    λt = λMcurve(Mtest) * KK^5/(Mtest)^5 #scaling factor to check with GW results
    
    return λt, Rt
end;