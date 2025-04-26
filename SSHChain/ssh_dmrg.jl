using Random
using HDF5
using Base.Threads
include("sshmodel.jl")
include("../ExactDiagonal.jl")
include("../entanglement.jl")

nthreads() = 8
nthreadpools() = 2
BLAS.set_num_threads(4)

seeds = rand(RandomDevice(), UInt128)
rngs = [MersenneTwister(seeds + 1<<tid) for tid in 1:nthreads()] 


let 
    L, D = 40, 4
    Lhalf = L ÷ 2
    b = Lhalf
    num, nalg = 50, 6
    v, w = 0.5, 1.0
    vs = LinRange(0.0,2.5, num+1)

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1.0E-12)

    fcsites = siteinds("Fermion", L)
    fmsites = siteinds("Fermion", L)
    mmsites = siteinds("MF", L)
    mosites = siteinds("MF", Lhalf)
    mesites = siteinds("MF", Lhalf)
    jwsites = siteinds("S=1/2", L)

    energies = zeros(num+1, nalg)
    entropies = zeros(num+1, nalg-1)
    maxbonds = zeros(Int, (num+1, nalg-1))

    @threads for i = 1:num+1
        v = vs[i]
        rng = rngs[threadid()]

        fcpsi0 = random_mps(rng, fcsites; linkdims=D)
        fmpsi0 = random_mps(rng, fmsites; linkdims=D)
        mmpsi0 = random_mps(rng, mmsites; linkdims=D)
        m2psi0 = random_mps(rng, mosites, mesites; linkdims=D)
        jwpsi0 = random_mps(rng, jwsites; linkdims=D)

        fmobs = BondsObserver()
        fcobs = BondsObserver()
        mmobs = BondsObserver()
        mobs = BondsObserver()
        jwobs = BondsObserver()

        H_ssh = ssh_hamiltonian(L, v, w)
        H_fermion = ssh_fermion(fcsites, v, w)
        H_transf = ssh_transf(fmsites, v, w)
        H_origin = ssh_origin(mmsites, v, w)
        H_origin2 = ssh_origin2(mosites, mesites, v, w)
        H_jordan = ssh_jordan(jwsites, v, w)

        denergy = groundstate_energy(H_ssh)

        fcenergy, fcpsi = dmrg(H_fermion, fcpsi0, sw; observer=fcobs, outputlevel=0)
        fcentropy = entangle_entropy(fcpsi, b)
        fcbond = maximum(fcobs.bonds)

        fmenergy, fmpsi = dmrg(H_transf, fmpsi0, sw; observer=fmobs, outputlevel=0)
        fmentropy = entangle_entropy(fmpsi, b)
        fmbond = maximum(fmobs.bonds)
        
        mmenergy, mmpsi = dmrg(H_origin, mmpsi0, sw; observer=mmobs, outputlevel=0)
        mmentropy = entangle_entropy(mmpsi, b)
        mmbond = maximum(mmobs.bonds)

        menegy, mpsis = dmrg(H_origin2, m2psi0, sw; observer=mobs, outputlevel=0)
        mentropy = entangle_entropy(mpsis, b ÷ 2)
        mbond = maximum(mobs.bonds)
        
        jwenergy, jwpsi = dmrg(H_jordan, jwpsi0, sw; observer=jwobs, outputlevel=0)
        jwentropy = entangle_entropy(jwpsi, b)
        jwbond = maximum(jwobs.bonds)
        
        energies[i,:] .= [fcenergy, fmenergy, mmenergy, menegy, jwenergy, denergy]
        entropies[i,:] .= [fcentropy, fmentropy, mmentropy, mentropy, jwentropy]
        maxbonds[i,:] .= [fcbond, fmbond, mmbond, mbond, jwbond]
    end

    h5open("SSHChain/sshplotdata.h5", "w") do file
        write(file, "vs", collect(vs))
        write(file, "energies", energies)
        write(file, "entropies", entropies)
        write(file, "maxbonds", maxbonds)
    end

end
