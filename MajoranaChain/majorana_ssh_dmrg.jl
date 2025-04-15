using BenchmarkTools
using HDF5
include("majorana_ssh.jl")
include("../entanglement.jl")

sw=Sweeps(10)
setmaxdim!(sw,200)
setcutoff!(sw,1E-14)
krydim = 4

let 
    L, D = 40, 4
    w, b = 1.0, Int(L/2)
    num = 50

    vs = LinRange(0.0, 2.5, num + 1)
    energies = zeros((num + 1, 4))
    entropies = zeros((num+1, 4))
    maxbonds = zeros(Int, (num+1, 4))
    # sites and initial state
    fsites = siteinds("Fermion", L; conserve_qns = false)
    fpsi0 = random_mps(fsites; linkdims = D)
    msites = siteinds("MF", L ; conserve_qns = false)
    mpsi0 = random_mps(msites; linkdims = D)
    jsites = siteinds("S=1/2", L)
    jpsi0 = random_mps(jsites; linkdims = D)
    # different parameters
    for (i, v) in enumerate(vs)
        H_majorana = majorana_hamiltonian(msites, v, w)
        H_mfermion = majorana_hamiltonian(fsites, v, w)
        H_cfermion = fermion_hamiltonian(fsites, v, w)
        H_jordan = jordan_hamiltonian(jsites, v, w)

        menergy ,mpsi = dmrg(H_majorana, mpsi0, sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        fmenergy, fmpsi = dmrg(H_mfermion, fpsi0, sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        fcenergy, fcpsi = dmrg(H_cfermion, fpsi0, sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        jenergy, jpsi = dmrg(H_jordan, jpsi0, sw; eigsolve_krylovdim = krydim, outputlevel = 0)

        mentropy = entangle_entropy(mpsi, b)
        fmentropy = entangle_entropy(fmpsi, b)
        fcentropy = entangle_entropy(fcpsi, b)
        jentropy = entangle_entropy(jpsi, b)
        
        energies[i, :] .= [menergy, fmenergy, fcenergy, jenergy ]
        maxbonds[i, :] .= [maxlinkdim(mpsi), maxlinkdim(fmpsi), maxlinkdim(fcpsi), maxlinkdim(jpsi)]
        entropies[i, :] .= [mentropy, fmentropy, fcenergy, jentropy]
    end

    h5open("MajoranaChain/majoranasshplotdata.h5", "cw") do file
        write(file, "vs", collect(vs))
        write(file, "energies", energies)
        write(file, "entropies", entropies)
        write(file, "maxbonds", maxbonds)
    end
end