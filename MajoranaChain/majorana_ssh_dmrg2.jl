using BenchmarkTools
using Base.Threads
using HDF5
include("majorana_ssh.jl")
include("../entanglement.jl")

nthreads() = 4
nthreadpools() = 1
BLAS.set_num_threads(4)

sw = Sweeps(10)
setmaxdim!(sw, 200)
setcutoff!(sw, 1E-14)
krydim = 4

let
    L, D = 40, 4
    w, b = 1.0, L ÷ 2
    vs = [0.5, 1.0, 2.0]
    ns, nalg = length(vs), 4

    correlations = zeros(ns, L, nalg)
    entangles = zeros(ns, L+1, nalg)
    runtime = zeros(ns, nalg)
    usespace = zeros(ns, nalg) 

    # sites
    fsites = siteinds("Fermion", L; conserve_qns = false)
    msites = siteinds("MF", L ; conserve_qns = false)
    jsites = siteinds("S=1/2", L)

    @threads for i = 1:ns
        #initial states
        mpsi0 = random_mps(msites; linkdims = D)
        fpsi0 = random_mps(fsites; linkdims = D)
        jpsi0 = random_mps(jsites; linkdims = D)
        # make hamiltonian MPO
        v = vs[i]
        H_majorana = majorana_hamiltonian(msites, v, w)
        H_mfermion = majorana_hamiltonian(fsites, v, w)
        H_cfermion = fermion_hamiltonian(fsites, v, w)
        H_jordan = jordan_hamiltonian(jsites, v, w)

        fcbench = @btimed dmrg($H_cfermion, $fpsi0, $sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        fmbench = @btimed dmrg($H_mfermion, $fpsi0, $sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        mbench = @btimed dmrg($H_majorana, $mpsi0, $sw; eigsolve_krylovdim = krydim, outputlevel = 0 )
        jbench = @btimed dmrg($H_jordan, $jpsi0, $sw; eigsolve_krylovdim = krydim, outputlevel = 0)

        psis = [mbench.value[2], fmbench.value[2], fcbench.value[2],jbench.value[2]]
        corrs = zeros(L, nalg)
        entropies = zeros(L+1, nalg)
        for a in eachindex(psis)
            corrs[:, a] .= correlation_matrix(psis[a], "N", "N")[b, :]
            entropies[:, a] .= [entangle_entropy(psis[a], j) for j=0:L]
        end 
        correlations[i, :, :] .= corrs
        entangles[i, :, :] .= entropies
        runtime[i, :] .= [mbench.time, fmbench.time, fcbench.time, jbench.time]
        usespace[i, :] .= [mbench.bytes, fmbench.bytes, fcbench.bytes, jbench.bytes] ./1<<20
    end

    h5open("MajoranaChain/majoranasshbenchdata.h5", "w") do file
        write(file, "vs", vs)
        write(file, "correlations", correlations)
        write(file, "entangles", entangles)
        write(file, "runtime", runtime)
        write(file, "usespace", usespace)
    end
end