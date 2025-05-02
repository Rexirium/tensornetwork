using HDF5
using Base.Threads
include("kitaevchain.jl")
include("../ExactDiagonal.jl")
include("../entanglement.jl")

nthreads() = 8
nthreadpools() = 2
BLAS.set_num_threads(4)

let
    L, D = 40, 4
    num, nalg = 50, 4
    t, Δ = 1.0, 1.0
    mus = LinRange(0.0, 4.0, num+1)
    b = L ÷ 2

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1E-10)
    krydim = 4
    
    fcsites = siteinds("Fermion", L)
    fmsites = siteinds("Fermion", L)
    mmsites = siteinds("MF", L)
    jwsites = siteinds("S=1/2", L)

    energies = zeros(num+1, nalg+1)
    entropies = zeros(num+1, nalg)
    maxbonds = zeros(Int, (num+1, nalg))
    @threads for i in 1:(num+1)
        μ = mus[i]

        fcpsi0 = random_mps(fcsites; linkdims=D)
        fmpsi0 = random_mps(fmsites; linkdims=D)
        mmpsi0 = random_mps(mmsites; linkdims=D)
        jwpsi0 = random_mps(jwsites; linkdims=D) 

        fmobs = BondsObserver()
        fcobs = BondsObserver()
        mmobs = BondsObserver()
        jwobs = BondsObserver()

        A, B = KitaevChain_BdG(L, μ, t, Δ)
        H_fermion = KitaevChain_fermion(fcsites, μ, t, Δ)
        H_transf = KitaevChain_transf(fmsites, μ, t, Δ)
        H_origin2 = KitaevChain_origin2(mmsites, μ, t, Δ)
        H_jordan = KitaevChain_jordan(jwsites, μ, t, Δ)

        denergy = groundstate_energy(A, B)

        fcenergy, fcpsi = dmrg(H_fermion, fcpsi0, sw; observer = fcobs, outputlevel=0)
        fcentropy = entangle_entropy(fcpsi, b)
        fcbond = maximum(fcobs.bonds)
        
        fmenergy, fmpsi = dmrg(H_transf, fmpsi0, sw; observer = fmobs, outputlevel=0)
        fmentropy = entangle_entropy(fmpsi, b)
        fmbond = maximum(fmobs.bonds)

        mmenergy, mmpsi = dmrg(H_origin2, mmpsi0, sw; observer = mmobs, outputlevel=0)
        mmentropy = entangle_entropy(mmpsi, b)
        mmbond = maximum(mmobs.bonds)

        jwenergy, jwpsi = dmrg(H_jordan, jwpsi0, sw; observer = jwobs, outputlevel=0)
        jwentropy = entangle_entropy(jwpsi, b)
        jwbond = maximum(jwobs.bonds)

        energies[i,:] .= [fcenergy, fmenergy, mmenergy, jwenergy, denergy]
        entropies[i, :] .= [fcentropy, fmentropy, mmentropy, jwentropy]
        maxbonds[i,:] .= [fcbond, fmbond, mmbond, jwbond]
    end

    h5open("KitaevChain/kcplotdata.h5", "w") do file
        write(file, "mus", collect(mus))
        write(file, "energies", energies)
        write(file, "entropies", entropies)
        write(file, "maxbonds", maxbonds)
    end
end  

