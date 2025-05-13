using LinearAlgebra, MKL
using HDF5
using LaTeXStrings
using Plots
default(
    grid=false, 
    titlelocation=:left,
    framestyle=:box,
    guidefontsize=16,
    legendfontsize=10,
    tickfontsize=10
)

let 
    file = h5open("KitaevChain/kcplotdata2.h5", "r")
    mus = read(file, "mus")
    energies = read(file, "energies")
    entropies = read(file, "entropies")
    close(file)
    num, nalg = size(energies)
    num -= 1
    nalg -= 1

    energyerr = zeros(num+1, nalg)
    for i in 1:nalg
        energyerr[:, i] .= max.(abs.((energies[:,i] .- energies[:,nalg+1])./ energies[:,nalg+1]), 1.0E-16)
    end

    alglabel = ["CF" "MF(trans)" "MF(origin)" "JW" ]
    alllabel = hcat(alglabel, "BdG")

    p1 = plot(mus, energies, xtick=:in, ylabel=L"E", lw = 2,  label = alllabel)
    vline!([2.0], line=(1,:dash), label=false,legend_position=:bottomleft)

    p1e = plot(mus, energyerr, xlabel=L"\mu/t", ylabel="err",yscale=:log10, ylim=(1.0E-14,1.0E-4),
        lw =2, label=alglabel)
    vline!([2.0], line=(1,:dash), label=false)

    p2 = plot(mus, entropies, xlabel=L"\mu/t", ylabel=L"S_\mathrm{vN}", 
        lw=2, label=alglabel)
    vline!([2.0], line=(1,:dash), label=false)

    P = plot(p1, p1e, p2, layout = @layout([[a{0.55h};b] c]), size=(800, 450), title=["a)" "b)" "c)"], 
        leftmargin=2Plots.mm, bottommargin=2.4Plots.mm)
    savefig(P, "KitaevChain/kitaevfigs/kc_en.svg")
end