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
    file = h5open("SSHChain/sshplotdata.h5", "r")
    vs = read(file, "vs")
    num = length(vs)-1
    energies = read(file, "energies")
    entropies = read(file, "entropies")
    close(file)
    nalg = 6

    energyerr = zeros(num+1, nalg-1)
    for i in 1:(nalg-1)
        energyerr[:, i] .= abs.((energies[:,i] .- energies[:,nalg])./ energies[:,nalg])
    end

    alglabel = ["CF" "MF(trans)" "MF(origin)" "MF(simpl)" "JW" ]
    alllabel = hcat(alglabel, "ED")

    p1 = plot(vs, energies, xtick=:in, ylabel=L"E", lw = 2,  label = alllabel)
    vline!([1.0], line=(1,:dash), label=false,legend_position=:bottomleft)

    p1e = plot(vs, energyerr, xlabel=L"v", ylabel="err",yscale=:log10, 
        lw =2, label=alglabel)
    vline!([1], line=(1,:dash), label=false)

    p2 = plot(vs, entropies, xlabel=L"v", ylabel=L"S_\mathrm{vN}",ylim=(-0.05,0.8), 
        lw=2, label=alglabel)
    vline!([1], line=(1,:dash), label=false)

    P = plot(p1, p1e, p2, layout = @layout([[a;b] c]), size=(800, 450), title=["a)" "b)" "c)"])
    #savefig(P, "SSHChain/sshfigs/ssh_en.pdf")
end