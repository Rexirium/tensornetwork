using LinearAlgebra
using Plots, LaTeXStrings
include("../ExactDiagonal.jl")
default(
    grid=true, 
    framestyle=:box,
    guidefontsize=14,
    legendfontsize=13,
    tickfontsize=10
)

let 
    L, num = 20, 50
    t, Δ = 1.0, 1.0
    mus = LinRange(0.0,4.0, num+1) .* t
    mu_list = [0.5, 2.0, 4.0] .* t
    ns = length(mu_list)

    Λ = diagm(1=>Δ .* ones(L-1), -1=>(-Δ) .* ones(L-1))
    Ξ = diagm(1=>(-t).* ones(L-1), -1=>(-t').*ones(L-1))

    energies = zeros(num+1, 2L)
    gsstates = zeros(L, ns)
    exstates = zeros(L, ns)
    for (i, μ) in enumerate(mus)
        Ξ[diagind(Ξ)] .= - μ
        energies[i, :] .= spectrum_BdG(Ξ, Λ; retstate=false)
    end
    legs_gs = String[]
    legs_ex = String[]
    for (i, μ) in enumerate(mu_list)
        Ξ[diagind(Ξ)] .= - μ
        λ, U, V = spectrum_BdG(Ξ, Λ; retstate=true)
        gsstates[:, i] .= abs2.(U[:, 1]) .+ abs2.(V[:, 1])
        exstates[:, i] .= abs2.(U[:, 2]) .+ abs2.(V[:, 2])
        push!(legs_gs, "zero mode "*L"\mu = "*"$(μ)")
        push!(legs_ex, "1st excited "*L"\mu = "*"$(μ)")
    end
    legs_gs = reshape(legs_gs, 1, :)
    legs_ex = reshape(legs_ex, 1, :)
    colors = [:red :green :blue]
    p1 = plot(mus, energies, leg=false, c=:black,
        xlabel=L"\mu", ylabel=L"E/t")
    vline!([2.0], line=(1,:dash))

    p2 = plot(gsstates, lw=1.5, ylim=(-0.02,0.4), label=legs_gs, c=colors,
    xlabel="site "*L"j", ylabel=L"n_j")
    plot!(exstates, line=(1.5, :dash), label=legs_ex, c=colors, leg=:top)

    p = plot(p1, p2, layout=@layout([a b]), size=(800, 400), leftmargin=4Plots.mm, bottommargin=3Plots.mm)
    savefig(p, "KitaevChain/kitaevfigs/kitaev_spectrum.svg")
end

#=
let
    L, num = 20, 50
    t, Δ = 1.0, 1.0
    mus = LinRange(0.0,4.0, num+1) .* t

    Ξ = zeros(L, L)
    Λ = diagm(1=>(Δ-t).*ones(L-1), -1=>(-t-Δ).*ones(L-1))
    X = [Ξ Λ ; -transpose(Λ) -Ξ]

    energies = zeros(num+1, 2L)
    for (i,μ) in enumerate(mus)
        X[diagind(X, L)] .= - μ
        X[diagind(X, -L)] .= μ
        F = schur(X)
        energies[i, :] .= sort(imag.(F.values))
    end
    plot(mus, energies, framestyle=:box, leg=false, c=:black,
        xlabel=L"\mu", ylabel=L"E/t")
    vline!([2.0], line=(1,:dash))
end
=#