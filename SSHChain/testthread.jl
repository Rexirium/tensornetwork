using Base.Threads
using LinearAlgebra, MKL 
using BenchmarkTools
include("sshmodel.jl")
#include("../ExactDiagonal.jl")
#include("../entanglement.jl")

nthreads()=8

let 
    L, D = 40, 4
    Lhalf = L ÷ 2
    b = Lhalf
    v, w = 0.5, 1.0

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1.0E-12)
    krydim = 4

    fcsites = siteinds("Fermion", L)
    mosites = siteinds("MF", Lhalf)
    mesites = siteinds("MF", Lhalf)

    fcpsi0 = random_mps(fcsites; linkdims=D)
    mpsi0 = random_mps(mosites, mesites; linkdims=D)

    H_fermion = ssh_fermion(fcsites, v, w)
    H_majorana = ssh_origin2(mosites, mesites, v, w)
    
    #println("simplified majorana fermion representation:")
    #@benchmark menergy, mpsi = dmrg($H_majorana, $mpsi0, $sw; eigsolve_krylovdim=$krydim, outputlevel=0)
    #println("==========================================================")
    println("complex fermion representation")
    @benchmark fcenergy, fcpsi = dmrg($H_fermion, $fcpsi0, $sw; eigsolve_krylovdim=$krydim, outputlevel=0)
    #println("==========================================================")
end

