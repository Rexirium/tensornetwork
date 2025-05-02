using ITensors, ITensorMPS
using LinearAlgebra, MKL
include("../MajoranaRep.jl")

function KitaevChain_BdG(Ls::Int, mu::Real, tt::Number, delta::Number)
    A = diagm(0=>fill(-mu, Ls), 1=>fill(-tt, Ls-1), -1=>fill(-tt', Ls-1))
    B = diagm(1=>fill(delta, Ls-1), -1=>fill(-delta, Ls-1))
    return A, B
end

function KitaevChain_fermion(s::Vector{<:Index}, mu::Real, tt::Number, delta::Number)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls
        os+= -mu, "N", j
        os+= mu/2, "Id", j
    end  
    for j in 1:Ls-1
        os+= -tt, "Cdag", j, "C", j+1
        os+= -tt, "Cdag", j+1, "C", j
        os+= delta, "C",j, "C", j+1
        os+= delta',"Cdag",j+1, "Cdag", j
    end      
    return MPO(os,s)
end

function KitaevChain_fermion(s::Vector{<:Index}, mu::Real, tt::Number, delta::Number, V::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls
        os += -mu, "N", j
        os += mu/2, "Id", j
    end
    for j in 1:Ls-1
        os += -tt, "Cdag", j, "C", j+1
        os += -tt, "Cdag", j+1, "C", j
        os += delta, "C", j, "C", j+1
        os += delta', "Cdag", j+1, "Cdag", j
        os += V, "N", j, "N", j+1
    end
    return MPO(os,s)
end

function KitaevChain_transf(s::Vector{<:Index}, mu::Real, tt::Real, delta::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls
        os += -im*mu/2, "Gamma1", j, "Gamma2", j
    end
    for j in 1:Ls-1
        os += im*(tt+delta)/2, "Gamma2", j, "Gamma1", j+1
        os += im*(-tt+delta)/2, "Gamma1", j, "Gamma2", j+1
    end
    return MPO(os,s)
end

function KitaevChain_origin(s::Vector{<:Index}, mu::Real, tt::Real, delta::Real)
    Ls = length(s)
    os, oso, ose = OpSum(), OpSum(), OpSum()
    a, b = delta - tt, delta + tt
    for j in 1:Ls-1
        os += -im*mu/2, "Gamma1", j, "Gamma2", j
        if isodd(j)
            oso += im*a/2, "Gamma1", j, "Gamma2", j+1
            ose += im*b/2, "Gamma2", j, "Gamma1", j+1
        else
            oso += im*b/2, "Gamma2", j, "Gamma1", j+1
            ose += im*a/2, "Gamma1", j, "Gamma2", j+1
        end
    end
    os += -im*mu/2, "Gamma1", Ls, "Gamma2", Ls
    H0 = MPO(os, s)
    Ho = MPO(oso, s)
    He = MPO(ose, s)
    return [H0, Ho, He]
end

function KitaevChain_origin2(s::Vector{<:Index}, mu::Real, tt::Real, delta::Real)
    Ls = length(s)
    oso, ose = OpSum(), OpSum()
    a, b = delta - tt, delta + tt
    for j in 1:Ls-1
        oso += -im*mu/4, "Gamma1", j, "Gamma2", j
        oso += im*b/2, "Gamma2", j, "Gamma1", j+1
        ose += im*mu/4, "Gamma2",j, "Gamma1", j
        ose += im*a/2, "Gamma1", j, "Gamma2", j+1
    end
    oso += -im*mu/4, "Gamma1", Ls, "Gamma2", Ls
    ose += im*mu/4, "Gamma2", Ls, "Gamma1", Ls
    Ho = MPO(oso, s)
    He = MPO(ose, s)
    return [Ho, He]
end

function KitaevChain_jordan(s::Vector{<:Index}, mu::Real, tt::Number, delta::Number)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls
        os += -mu, "S-", j, "S+", j
        os += mu/2, "Id", j
    end
    for j in 1:Ls-1
        os += -tt, "S-", j, "S+", j+1
        os += -tt, "S-", j+1, "S+", j
        os += -delta, "S+", j, "S+", j+1
        os += -delta', "S-", j+1, "S-", j
    end
    return MPO(os,s)
end
