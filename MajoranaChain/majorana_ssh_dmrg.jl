using Base.Threads
using HDF5
include("majorana_ssh.jl")
include("../entanglement.jl")

nthreads() = 8
nthreadpools() = 1
BLAS.set_num_threads(4)

sw = Sweeps(10)
setmaxdim!(sw, 200)
setcutoff!(sw, 1E-12)
krydim = 4

let 
    L, D = 40, 4
    w, b = 1.0, L÷2
    num = 50

    vs = LinRange(0.0, 2.5, num + 1)
    energies = zeros(num + 1, 5)
    entropies = zeros(num+1, 4)
    maxbonds = zeros(Int, (num+1, 4))
    # sites and initial state
    fsites = siteinds("Fermion", L; conserve_qns = false)
    fpsi0 = random_mps(fsites; linkdims = D)
    msites = siteinds("MF", L ; conserve_qns = false)
    mpsi0 = random_mps(msites; linkdims = D)
    jsites = siteinds("S=1/2", L)
    jpsi0 = random_mps(jsites; linkdims = D)

    # different parameters
    @threads for i = 1:(num+1)
        v = vs[i]
        H_majorana = majorana_hamiltonian(msites, v, w)
        H_mfermion = majorana_hamiltonian(fsites, v, w)
        H_cfermion = fermion_hamiltonian(fsites, v, w)
        H_jordan = jordan_hamiltonian(jsites, v, w)
        HBdG = majorana2BdG(L, v, w)

        mobs = BondsObserver()
        fmobs = BondsObserver()
        fobs = BondsObserver()
        jobs = BondsObserver()

        denergy = sum(eigvals(HBdG, 1:L))./2
        menergy ,mpsi = dmrg(H_majorana, mpsi0, sw; observer=mobs, eigsolve_krylovdim = krydim, outputlevel = 0 )
        fmenergy, fmpsi = dmrg(H_mfermion, fpsi0, sw; observer=fmobs, eigsolve_krylovdim = krydim, outputlevel = 0 )
        fcenergy, fcpsi = dmrg(H_cfermion, fpsi0, sw; observer=fobs, eigsolve_krylovdim = krydim, outputlevel = 0 )
        jenergy, jpsi = dmrg(H_jordan, jpsi0, sw; observer=jobs, eigsolve_krylovdim = krydim, outputlevel = 0)

        mentropy = entangle_entropy(mpsi, b)
        fmentropy = entangle_entropy(fmpsi, b)
        fcentropy = entangle_entropy(fcpsi, b)
        jentropy = entangle_entropy(jpsi, b)
        
        energies[i, :] .= [menergy, fmenergy, fcenergy, jenergy, denergy]
        maxbonds[i, :] .= [maximum(mobs.bonds), maximum(fmobs.bonds), maximum(fobs.bonds), maximum(jobs.bonds)]
        entropies[i, :] .= [mentropy, fmentropy, fcentropy, jentropy]
    end

    h5open("MajoranaChain/majoranasshplotdata.h5", "w") do file
        write(file, "vs", collect(vs))
        write(file, "energies", energies)
        write(file, "entropies", entropies)
        write(file, "maxbonds", maxbonds)
    end
end