using ITensors, ITensorMPS
using LinearAlgebra
using Plots, LaTeXStrings

include("ExactDiagonal.jl")
include("entanglement.jl")

function hamiltonian_mpo(s::Vector{Index{Int}}, hamil::AbstractMatrix)
    length(s) == size(hamil, 1) || error("Wrong length of MPS indices")
    Ns = length(s)
    os = OpSum()
    for i in 1:Ns, j in i:Ns
        if hamil[i, j] == 0.0
            continue
        end
        if i == j
            os += hamil[i, i], "N", i
        else
            os += hamil[i, j], "Cdag", i, "C", j
            os += hamil[i, j]', "Cdag", j, "C", i
        end
    end
    return MPO(os, s)
end
#=
let 
    L, D = 10, 4
    v, w = 0.5, 1.0

    sw =Sweeps(15)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1.0E-14)
    krydim = 4
    obs = BondsObserver()

    sites = siteinds("Fermion", L)
    psi0 = random_mps(sites; linkdims = D)
    H = Hermitian(rand(L, L))
    H_MPO = hamiltonian_mpo(sites, H)

    energy_DMRG, psi_DMRG = dmrg(H_MPO, psi0, sw, observer = obs; eigsolve_krylovdim = krydim)
    energy_ED, psi_ED = groundstate(H)
    maxbond = max((obs.bonds)...)
    spectrum = spectrum_BdG(H, retstate = false)

    density_DMRG = expect(psi_DMRG, "N")
    density_ED = density_vec(psi_ED[1])
    
    
    println("DMRG energy is $energy_DMRG, ED energy is $energy_ED")
    println("maxbond during DMRG is $maxbond")
    p = scatter(spectrum, xlabel = L"n", ylabel="energy", framestyle = :box, legend = false, title="spectrum")
    b1 = bar(density_DMRG, ylabel = "density", 
        legend = false, framestyle = :box, color = :red, title = "DMRG density")
    b2 = bar(density_ED, xlabel = "site "*L"n", ylabel = "density", 
        legend = false, framestyle = :box, color = :blue, title = "ED density")
    plot(p, b1, b2, layout = @layout([a  [b; c]]), size = (1000, 600))
end
=#

function makeKitaevChain(Ls::Int, mu::Real, tt::Number, delta::Number)
    arrt = -tt .* ones(Ls-1)
    arrtc = -tt' .* ones(Ls-1)
    arrd = delta .* ones(Ls-1)
    A = diagm(0=>(-mu .* ones(Ls)), 1=>arrt, -1=>arrtc)
    B = diagm(1=>arrd, -1=>(-arrd))
    return A, B
end


let 
    L, num = 40, 50
    μ, t, Δ = 0.05, 1.0, 1.0
    mus = LinRange(0.0,4.0, num+1)
    
    A, B = makeKitaevChain(L, μ, t, Δ)
    spec, U, V = spectrum_BdG(A, B; retstate =true)
    n = 1
    bar(abs2.(U[:, n]),framestyle=:box,lw=2, leg=false)
end