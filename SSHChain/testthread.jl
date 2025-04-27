using Base.Threads
using LinearAlgebra, MKL
include("../ExactDiagonal.jl")
include("sshmodel.jl")
include("../entanglement.jl")


let 
    L, D = 40, 4
    Lhalf = L ÷ 2
    b = Lhalf
    v, w = 0.9, 1.0

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1.0E-12)
    krydim = 4

    
    fcsites = siteinds("Fermion", L)
    fcpsi0 = random_mps(fcsites; linkdims=D)

    H_fermion = ssh_fermion(fcsites, v, w)
    
    fcenergy, fcpsi = dmrg(H_fermion, fcpsi0, sw; eigsolve_krylovdim=krydim)
    linkdims(fcpsi)

end

