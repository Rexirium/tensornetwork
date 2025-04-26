using Random
using HDF5
using Base.Threads
using BenchmarkTools
include("sshmodel.jl")
include("../ExactDiagonal.jl")
include("../entanglement.jl")

nthreads() = 4
nthreadpools() = 2
BLAS.set_num_threads(4)

seeds = rand(RandomDevice(), UInt128)
rngs = [MersenneTwister(seeds + 1<<tid) for tid in 1:nthreads()] 


let 
    L, D = 40, 4
    Lhalf = L ÷ 2
    nalg = 5
    v, w = 0.5, 1.0
    vs = [0.0, 0.5, 1.0, 1.5,  2.0]
    ns = length(vs)

    sw = Sweeps(10)
    setmaxdim!(sw, 300)
    setcutoff!(sw, 1.0E-14)
    krydim = 4

    fcsites = siteinds("Fermion", L)
    fmsites = siteinds("Fermion", L)
    mmsites = siteinds("MF", L)
    mosites = siteinds("MF", Lhalf)
    mesites = siteinds("MF", Lhalf)
    jwsites = siteinds("S=1/2", L)

    entangles = zeros(ns, Lhalf+1, nalg)
    runtime = zeros(ns, nalg)
    usespace = zeros(ns, nalg) 

    @threads for i = 1:ns
        v = vs[i]
        rng = rngs[threadid()]
        #initial states
        fcpsi0 = random_mps(rng, fcsites; linkdims=D)
        fmpsi0 = random_mps(rng, fmsites; linkdims=D)
        mmpsi0 = random_mps(rng, mmsites; linkdims=D)
        mpsi0 = random_mps(rng, mosites, mesites; linkdims=D)
        jwpsi0 = random_mps(rng, jwsites; linkdims=D)
        # make hamiltonian MPO
        H_fermion = ssh_fermion(fcsites, v, w)
        H_transf = ssh_transf(fmsites, v, w)
        H_origin = ssh_origin(mmsites, v, w)
        H_origin2 = ssh_origin2(mosites, mesites, v, w)
        H_jordan = ssh_jordan(jwsites, v, w)


        entropies = zeros(Lhalf+1, nalg)
        fcbench = @btimed dmrg($H_fermion, $fcpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        entropies[:,1] .= [entangle_entropy(fcbench.value[2], 2c) for c = 0:Lhalf]
        fmbench = @btimed dmrg($H_transf, $fmpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        entropies[:,2] .= [entangle_entropy(fmbench.value[2], 2c) for c = 0:Lhalf]
        mmbench = @btimed dmrg($H_origin, $mmpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        entropies[:,3] .= [entangle_entropy(mmbench.value[2], 2c) for c = 0:Lhalf]
        mbench = @btimed dmrg($H_origin2, $mpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        entropies[:,4] .= [entangle_entropy(mbench.value[2], c) for c = 0:Lhalf]
        jwbench = @btimed dmrg($H_jordan, $jwpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0)
        entropies[:,5] .= [entangle_entropy(jwbench.value[2], 2c) for c = 0:Lhalf]

        entangles[i, :, :] .= entropies
        runtime[i, :] .= [fcbench.time, fmbench.time, mmbench.time, mbench.time, jwbench.time]
        usespace[i, :] .= [fcbench.bytes, fmbench.bytes, mmbench.bytes, mbench.bytes, jwbench.bytes] ./1<<20
    end

    h5open("SSHChain/sshbenchdata.h5", "w") do file
        write(file, "vs", vs)
        write(file, "entangles", entangles)
        write(file, "runtime", runtime)
        write(file, "usespace", usespace)
    end


end