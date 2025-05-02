using HDF5
using Base.Threads
using BenchmarkTools
include("kitaevchain.jl")
include("../ExactDiagonal.jl")
include("../entanglement.jl")

nthreads() = 4
nthreadpools() = 2
BLAS.set_num_threads(4)

let 
    L, D = 40, 4
    num, nalg = 50, 4
    t, Δ = 1.0, 1.0
    muls = [0.0, 1.0, 2.0, 3.0, 4.0]
    ns = length(muls)
    b = L÷2

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1E-12)
    krydim = 4
    
    fcsites = siteinds("Fermion", L)
    fmsites = siteinds("Fermion", L)
    mmsites = siteinds("MF", L)
    jwsites = siteinds("S=1/2", L)

    entangles = zeros(ns, L+1, nalg)
    correlations = zeros(ns, L, nalg)
    runtime = zeros(ns, nalg)
    usespace = zeros(ns, nalg) 

    @threads for  i in 1:ns
        μ = muls[i]

        fcpsi0 = random_mps(fcsites; linkdims=D)
        fmpsi0 = random_mps(fmsites; linkdims=D)
        mmpsi0 = random_mps(mmsites; linkdims=D)
        jwpsi0 = random_mps(jwsites; linkdims=D) 

        H_fermion = KitaevChain_fermion(fcsites, μ, t, Δ)
        H_transf = KitaevChain_transf(fmsites, μ, t, Δ)
        H_origin2 = KitaevChain_origin2(mmsites, μ, t, Δ)
        H_jordan = KitaevChain_jordan(jwsites, μ, t, Δ)

        fcbench = @btimed dmrg($H_fermion, $fcpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        fmbench = @btimed dmrg($H_transf, $fmpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        mmbench = @btimed dmrg($H_origin2, $mmpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0 )
        jwbench = @btimed dmrg($H_jordan, $jwpsi0, $sw; eigsolve_krylovdim = $krydim, outputlevel = 0)

        psis = [fcbench.value[2], fmbench.value[2], mmbench.value[2],jwbench.value[2]]
        corrs = zeros(L, nalg)
        entropies = zeros(L+1, nalg)
        for a in 1:nalg
            corrs[:, a] .= correlation_matrix(psis[a], "N", "N")[b, :]
            entropies[:, a] .= [entangle_entropy(psis[a], j) for j=0:L]
        end 
        correlations[i, :, :] .= corrs
        entangles[i, :, :] .= entropies
        runtime[i, :] .= [fcbench.time, fmbench.time, mmbench.time, jwbench.time]
        usespace[i, :] .= [fcbench.bytes, fmbench.bytes, mmbench.bytes, jwbench.bytes] ./(1<<20)
    end

    h5open("KitaevChain/kcbenchdata.h5", "w") do file
        write(file, "muls", muls)
        write(file, "correlations", correlations)
        write(file, "entangles", entangles)
        write(file, "runtime", runtime)
        write(file, "usespace", usespace)
    end
end