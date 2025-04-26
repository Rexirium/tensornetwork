using ITensors, ITensorMPS
using LinearAlgebra, MKL
include("../MajoranaRep.jl")

function ssh_hamiltonian(Ls::Int, t1::Number, t2::Number)
    arr = repeat([t1, t2], Ls ÷ 2)[1:Ls-1]
    return Hermitian(diagm(1=>arr))
end

function ssh_fermion(s::Vector{<:Index}, t1::Number, t2::Number)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        t = isodd(j) ? t1 : t2
        os += t, "Cdag", j, "C", j+1
        os += t, "Cdag", j+1, "C", j
    end
    return MPO(os,s)
end

function ssh_jordan(s::Vector{<:Index}, t1::Number, t2::Number)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        t = isodd(j) ? t1 : t2
        os += t, "S-", j, "S+", j+1
        os += t, "S-", j+1, "S+", j
    end
    return MPO(os, s)
end

function ssh_transf(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    os = OpSum()
    for j in 1:Ls-1
        t = isodd(j) ? t1 : t2
        os += im*t/2, "Gamma1", j, "Gamma2", j+1
        os += im*t/2, "Gamma1", j+1, "Gamma2", j
    end
    return MPO(os, s)
end

function ssh_origin(s::Vector{<:Index}, t1::Real, t2::Real)
    Ls = length(s)
    oso, ose = OpSum(), OpSum()
    for j in 1:Ls-1
        if isodd(j)
            oso += im*t1/2, "Gamma1", j, "Gamma2", j+1
            ose += -im*t1/2, "Gamma2", j, "Gamma1", j+1
        else
            oso += -im*t2/2, "Gamma2", j, "Gamma1", j+1
            ose += im*t2/2, "Gamma1", j, "Gamma2", j+1
        end
    end
    Ho = MPO(oso, s)
    He = MPO(ose, s)
    return [Ho, He]
end

function ssh_origin2(so::Vector{<:Index}, se::Vector{<:Index}, t1::Real, t2::Real)
    Lh = length(so)
    oso, ose = OpSum(), OpSum()
    for j in 1:Lh-1
        oso += im*t1/2, "Gamma1", j, "Gamma2", j
        oso += -im*t2/2, "Gamma2", j, "Gamma1", j+1
        ose += -im*t1/2, "Gamma2", j, "Gamma1", j
        ose += im*t2/2, "Gamma1", j, "Gamma2", j+1
    end
    oso += im*t1/2, "Gamma1", Lh, "Gamma2", Lh
    ose += -im*t1/2, "Gamma2", Lh, "Gamma1", Lh
    Ho = MPO(oso, so)
    He = MPO(ose, se)
    return [Ho, He]
end

