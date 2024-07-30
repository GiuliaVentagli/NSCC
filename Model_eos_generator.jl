using Interpolations

#function to find (closest) position in array of certain value 
#!!It works for increasing values on column!!
function findpos(matrix,value,column)
    K = length(matrix[:,column])
    
    #check if element is in list
    if value > matrix[end,column] || value < matrix[1,column] 
        return 0
    end
    
    for b = 1:length(matrix[:,column])
        Δx1 = abs(matrix[b,column]  - value)
        Δx2 = abs(matrix[b+1,column] - value)
        if matrix[b+1,column] >= value
            if Δx1 < Δx2
                return b
            else
                return b+1
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
    #@show(ρold)
    ϵold = eos[k,2]
    pold = eos[k,1]
    
    for j = 1:M-1
        ρnew = ρold + Δρ
        #@show(ρnew,Δρ)
        ϵnew = ϵold + Δρ * (ϵold + pold)/ρold
        pnew = pold + cs_interp(ρold)^2 * Δρ * (ϵold + pold)/ρold #cs^2=dp/dϵ
        if ρnew > ρfin #break cycle if we reach last mass density
            break
        end
        he_eos = vcat(he_eos, [pnew ϵnew ρnew])
        ρold = ρnew #store new data in old variable for next loop
        ϵold = ϵnew
        pold = pnew
        if ρnew < 10.0^-3
            Δρ = 10.0^-5
            #Δρ = 1.0 * 10.0^-5
        elseif ρnew < 2.0 * 10.0^-3
            Δρ = 10.0^-4
        elseif ρnew < 10.0^-2
            #Δρ = 10.0^-4
            Δρ = 2.0 * 10.0^-4
        elseif ρnew < 10^-1
            #Δρ = 10.0^-3
            Δρ = 2.0 * 10.0^-3
        else 
            #Δρ = 10.0^-2
            Δρ = 2.0 * 10.0^-2
        end
    end
    return he_eos
end

#merge two EOSs
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
    k = findpos(eos,ρt,3) #find closest position of ρt
    #@show eos[k,1] #shows pressure at QCD threshold

    #Evaluating derivative functions to find quantities at threshold
    p1 = derThr(eos, k)
   
    #creating an interpolated function for speed of sound modification
    cs_interp = csInterp(eos, k, p1, matr_ρ, matr_cs, ds)
    
    #define matrix to store new part of eos
    he_eos = Matrix{Float64}(undef,0,3)
    he_eos = heEOS(he_eos, eos, k, cs_interp)
    
    #merge the two eos at right threshold
    eos_matrix = mergeEOS(eos_matrix, he_eos, eos, k)
    return eos_matrix
end

#Introducing Λ transition
function Λtransition(eos_complete::AbstractMatrix{T}, matr) where T
    k = findpos(matr,pc,1) #find line in matrix corrisponding to pc
    
    #NB: p<pc total=fluid, p>pc total=fluid+Λ contribution
    for j=1:length(matr[:,1])
        if j < k
            eos_complete = vcat(eos_complete, [matr[j,1] matr[j,2]])
        elseif j >= k
            pfluid = matr[j,1] + Λ
            #if it goes outside of tables, it stops.
            if Λ >= 0 && pfluid > matr[end,1] 
                break
            end
            if Λ < 0 && pfluid < matr[1,1]
                break
            end
            #find closest position of pfluid in tables, and corresponding energy density
            l = findpos(matr,pfluid,1)
            eos_complete = vcat(eos_complete, [matr[j,1] matr[l,2]+Λ])
        end
    end
    return eos_complete
end;