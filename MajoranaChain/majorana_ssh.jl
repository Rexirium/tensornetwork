using LinearAlgebra
using ITensors, ITensorMPS
include("../MajoranaRep.jl")

function majorana_BdG(Ns::Int, t1::Real, t2::Real, retstate::Bool = false)
    a, b = im*t1/2, im*t2/2
    arr = repeat([a,b], Ns)
    deleteat!(arr, 2*Ns)
    hamil = Hermitian(diagm(1=>arr))
    if retstate
        return eigen(hamil)
    else
        return eigvals(hamil)
    end
end

function majorana_hamiltonian(s::Vector{Index{Int}}, t1::Real, t2::Real)
    Ns = length(s)
    os = OpSum()
    for j in 1:Ns-1
        os += im*t1, "Gamma1", j, "Gamma2", j
        os += im*t2, "Gamma2", j, "Gamma1", j+1
    end
    os += im*t1, "Gamma1", Ns, "Gamma2", Ns
    return MPO(os,s)
end

function jordan_hamiltonian(s::Vector{Index{Int}}, t1::Real, t2::Real)
    Ns = length(s)
    os = OpSum()
    for j in 1:Ns-1
        os += -t1, "Z", j
        os += -t2, "X", j, "X", j+1
    end
    os += -t1, "Z", Ns
    return MPO(os,s)
end

function fermion_hamiltonian(s::Vector{Index{Int}}, t1::Real, t2::Real)
    Ns = length(s)
    os = OpSum()
    for j in 1:Ns
        os += 2*t1, "Cdag", j, "C", j
        os += -t1, "Id", j
    end
    for j in 1:Ns-1
        os += -t2, "Cdag", j, "C", j+1
        os += -t2, "Cdag", j+1, "C", j
        os += t2, "C", j, "C", j+1
        os += t2, "Cdag", j+1, "Cdag", j
    end
    return MPO(os, s)
end

function spin_hamiltonian(s::Vector{Index{Int}}, t1::Real, t2::Real)
    Ns = length(s)
    os = OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += t1, "X", j, "Y", j+1
        else
            os += t2, "Y", j, "X", j+1
        end
    end
    return MPO(os, s)
end

L, D = 40, 4
v, w = 1.0, 2.0
sw = Sweeps(15)
setmaxdim!(sw, 200)
setcutoff!(sw, 1E-14)
krydim=4

let 
    sites = siteinds("MF", L)
    psi0 = random_mps(sites; linkdims=D)
    H = majorana_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in Majorana basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    H=majorana_hamiltonian(sites, v, w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in Fermion basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    H=fermion_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites simplified Kitaev chain for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    sites=siteinds("S=1/2",L)
    psi0=random_mps(sites; linkdims=D)
    H=jordan_hamiltonian(sites, v, w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in spin basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end




