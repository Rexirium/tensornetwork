using LinearAlgebra, MKL
using ITensors, ITensorMPS
include("../MajoranaRep.jl")

function majorana2BdG(Ls::Int, t1::Real, t2::Real)
    A = diagm(0=>(2t1) .* ones(Ls), 1=>(-t2).*ones(Ls-1), -1=>(-t2').*ones(Ls-1) )
    B = diagm(1=>t2 .* ones(Ls-1), -1=>(-t2).*ones(Ls-1))
    H = [A -conj(B); B -transpose(A)]
    return Hermitian(H)
end

function majorana_hamiltonian(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        os += im*t1, "Gamma1", j, "Gamma2", j
        os += im*t2, "Gamma2", j, "Gamma1", j+1
    end
    os += im*t1, "Gamma1", Ls, "Gamma2", Ls
    return MPO(os,s)
end

function jordan_hamiltonian(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        os += -t1, "Z", j
        os += -t2, "X", j, "X", j+1
    end
    os += -t1, "Z", Ls
    return MPO(os,s)
end

function fermion_hamiltonian(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls
        os += 2*t1, "Cdag", j, "C", j
        os += -t1, "Id", j
    end
    for j in 1:Ls-1
        os += -t2, "Cdag", j, "C", j+1
        os += -t2, "Cdag", j+1, "C", j
        os += t2, "C", j, "C", j+1
        os += t2, "Cdag", j+1, "Cdag", j
    end
    return MPO(os, s)
end

function spin_hamiltonian(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        if isodd(j)
            os += t1, "X", j, "Y", j+1
        else
            os += t2, "Y", j, "X", j+1
        end
    end
    return MPO(os, s)
end



