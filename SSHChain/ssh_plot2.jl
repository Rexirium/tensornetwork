using LinearAlgebra, MKL
using HDF5
using LaTeXStrings
using Plots
default(
    grid=false, 
    titlelocation=:left,
    framestyle=:box,
    guidefontsize=14,
    legendfontsize=10,
    tickfontsize=10
)

let 
    file = h5open("MajoranaChain/majoranasshbenchdata.h5", "r")
    vs = read(file, "vs")
    entangles = read(file, "entangles")
    close(file)

    ns, Lc, nalg = size(entangles)
    Lc -= 1
    L = 2*Lc
    xs = range(0, Lc)

    alglabel = ["CF" "MF(trans)" "MF(origin)" "MF(simpl)" "JW" ]

    ps = []
    for i in 1:ns
        if i==4
            continue
        end
        p = plot(xs, entangles[i,:,:], xlabel="cell j", lw=1.5, label=alglabel,
            title = latexstring("v = $(vs[i]),\\, w = 1.0"))
        if i==1 || i==3
            ylabel!(L"S(j)")
        end
        push!(ps, p)
    end
    P = plot(ps... , layout= @layout([a b; c d]), size = (800, 600))

    savefig(P, "MajoranaChain/majoranafigs/mfssh_entangle.pdf")
end