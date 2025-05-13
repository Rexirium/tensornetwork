using Plots, LaTeXStrings
using Printf
using Base.Threads
nthreadpools() = 2
nthreads() = 8
BLAS.set_num_threads(4)

include("kitaevchain.jl")
include("../ExactDiagonal.jl")
include("../entanglement.jl")
default(
    grid=true, 
    framestyle=:box,
    guidefontsize=14,
    legendfontsize=10,
    tickfontsize=10
)

let 
    L, D = 40, 4
    μ, t, Δ = 3.0, -1.0 , 1.0
    num = 100
    mus = LinRange(-4.0, 4.0, num+1)
    b = L ÷ 2
    height = (abs(μ) + abs(t) + abs(Δ))*1.5

    energies = zeros(num+1, 2L)
    @threads for i in 1:(num+1)
        mu = mus[i]
        A, B = KitaevChain_BdG(L, mu, t, Δ)
        energies[i, :] .= spectrum_BdG(A, B; retstate=false)
    end

    A0, B0 = KitaevChain_BdG(L, μ, t, Δ)
    denergy = groundstate_energy(A0, B0)

    sw = Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1E-10)
    krydim = 4

    sites = siteinds("Fermion", L)
    psi0 = random_mps(sites; linkdims=D)

    H_mpo = KitaevChain_fermion(sites, μ, t, Δ)
    energy, psi = dmrg(H_mpo, psi0, sw; eigsolve_krylovdim=krydim)

    density = expect(psi, "N")
    corrs = correlation_matrix(psi, "N", "N")
    entangles = [entangle_entropy(psi, j) for j=0:L]
    
    pl = plot(mus, energies, xlabel=L"\mu", ylabel=L"E", leg=false, c=:black, title="spectrum")
    vline!([μ], line=(1.5,:dash))
    annotate!([(0.0, -4, @sprintf("ED E = %.4f", denergy)), 
        (0.0,4, @sprintf("DMRG E = %.4f", energy))])
    bm = bar(density, xlabel="site "*L"j", ylabel=L"n_j", ylim=(0.0,1.0), leg=false, title="density")
    hm = heatmap(transpose(corrs), xlabel=L"i", ylabel=L"j", title="correlation")
    pe = plot(0:L, entangles, lw=2, xlabel="site "*L"j", ylabel=L"S(j)", title="entropy", ylim=(0.0,1.0))
    plot(pl, bm, hm, pe, layout=@layout([a b; c d]), size=(1000, 800))
end

