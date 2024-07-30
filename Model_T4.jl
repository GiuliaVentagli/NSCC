using DifferentialEquations, Interpolations

# function to find (closest) position in array of certain value
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

#Point particle contribution
function Fpp(x)
    γ = 0.577126 #Euler's constant
    a0 = 1.
    a1 = 0.
    a2 = -743.0/336.0 - 11 * ν/4.
    a3 = 4. * π
    a4 = 34103.0/18144.0 + 13661.0/2016.0 * ν + 59.0/18.0 * ν^2
    a5 = - π * (4159.0/672.0 + 189.0/8.0 * ν)
    a6 = 16447322263.0/139708800.0 - 1712.0/105.0 * γ + 16.0/3. * π^2  -
        856.0/105.0 * log(16. * x) +  ν * (451.0/48.0 * π^2 - 56198689.0/217728.0) 
        + 541.0/896.0 * ν^2 - 5605.0/2592.0 * ν^3
    a7 = - π * (4415.0/4032.0 - 358675.0/6048.0 * ν -
        91495.0/1512.0 * ν^2)
    
    return (64. * ν)/(5. * M)  * x^5 * (a0 + a1 * x^(1/2) + a2 * x + a3 * x^(3/2) +
        a4 * x^2 + a5 * x^(5/2) + a6 * x^3 + a7 * x^(7/2)) 
end

#Tidal contribution
function Ftid(x)
    return (32. * Ma * λb) / (5. * M^7) * ( 12. * (1. + 11. * Ma/M) * x^10 +
    (4421.0/28.0 - 12263.0/28.0 * Mb/M + 1893.0/2.0 * (Mb/M)^2 - 661. * (Mb/M)^3 ) *x^11 ) +
    (32. * Mb * λa) / (5. * M^7) * ( 12. * (1. + 11. * Mb/M) * x^10 +
    (4421.0/28.0 - 12263.0/28.0 * Ma/M + 1893.0/2.0 * (Ma/M)^2 - 661. * (Ma/M)^3 ) *x^11 )
end

#Define the system of differential eqs to solve
function system!(du, u, p, t) # u = x, ϕ
    du[1] = Fpp(u[1]) + Ftid(u[1])
    du[2] = u[1]^(3/2) / M
end

function integrator()
    #for the initial ω I am using third Kepler's law
    ω0 = sqrt(G * M/(d^3))
    x0 = (G * M * ω0)^(2/3)
    ϕ0 = 0.
    tspan = (0, tf)
    u0 = [x0, ϕ0]
    
    condition(u, t, integrator) = G^(1/3)* M / u[1] < dend
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    prob = ODEProblem(system!, u0, tspan)
    
    sol = solve(prob, RK4(), callback = cb, dt = δt, adaptive = false)
    return sol
end

#Function to store data after solving differential equations
function dumpData(sol, name_file)
    #Store important data in first commented line
    dumping = open(name_file, "a");
    writedlm( dumping, ["#Ma = $Ma, Mb = $Mb, λa = $λa, λb = $λb, d_i = $d, d_f = $dend, d_obs = $dobs"])
    writedlm( dumping, ["#t/M x ϕ ω*M ϕ*M h22 Ih22"])
    close(dumping)
    for i=1:length(sol.t)-1
        tresc = sol.t[i]/M
        omegaresc = sol[1,i]^(3/2)
        phaseM = sol[2,i] * M
        
        h22 = coeff * sol[1,i] * cos(2 * sol[2,i])
        ih22 = coeff * sol[1,i] * sin(2 * sol[2,i])
        
        dumping = open(name_file, "a");
        writedlm( dumping,  [tresc sol[1,i] sol[2,i] omegaresc phaseM h22 ih22], ' ')
        close(dumping)
    end
end

#Function to compute the Newmann-Penrose scalar Ψ4
function strain(data)
    #Interpolate h_22 and Ih_22
    t_array = Interpolations.deduplicate_knots!(data[:,1];move_knots = true)
    h_22_interp = linear_interpolation(t_array,data[:,6])
    ih_22_interp = linear_interpolation(t_array,data[:,7])
    
    #Possibility to improve this part of code, not straightforward
    #Evaluate first order derivatives
    h22_prime(x) = only(Interpolations.gradient(h_22_interp, x))
    ih22_prime(x) = only(Interpolations.gradient(ih_22_interp, x));
    h22primeMat = Matrix{Float64}(undef,0,2)
    ih22primeMat = Matrix{Float64}(undef,0,2)
    for i = 1:length(data[:,1])
        val = h22_prime(data[i,1])
        ival = ih22_prime(data[i,1])
        h22primeMat = vcat(h22primeMat, [data[i,1] val])
        ih22primeMat = vcat(ih22primeMat, [data[i,1] ival])
    end
    
    h22prime_array = Interpolations.deduplicate_knots!(h22primeMat[:,1];move_knots = true)
    ih22prime_array = Interpolations.deduplicate_knots!(ih22primeMat[:,1];move_knots = true)
    h_22_prime_interp = linear_interpolation(h22prime_array,h22primeMat[:,2])
    ih_22_prime_interp = linear_interpolation(ih22prime_array,ih22primeMat[:,2])
    
    Reψ4(x) = only(Interpolations.gradient(h_22_prime_interp, x))
    Imψ4(x) = only(Interpolations.gradient(ih_22_prime_interp, x))
    
    return Reψ4, Imψ4
end;