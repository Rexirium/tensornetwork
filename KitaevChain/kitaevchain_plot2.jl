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
    file = h5open("KitaevChain/kcbenchdata.h5", "r")
    muls = read(file, "muls")
    entangles = read(file, "entangles")
    close(file)

    ns, L, nalg = size(entangles)
    L -= 1
    xs = range(0, L)

    alglabel = ["CF" "MF(trans)" "MF(origin)" "JW" ]
    ylims = [(0.0,0.1),(0.0,0.1),(0.0,0.6),(0.0,0.3),(0.0,0.3)]
    ps = []
    for i in 1:ns
        if i==1
            continue
        end
        p = plot(xs, entangles[i,:,:], lw=1.5, label=alglabel,
            ylim=ylims[i], title = latexstring("\\mu = $(muls[i]),\\, t= \\Delta = 1.0"),
            legend_position=(i==3 ? :bottom : :topright) )
        if i==2 || i==4
            ylabel!(L"S(x)")
        end
        if i==4 || i==5
            xlabel!("site "*(L"x"))
        end
        push!(ps, p)
    end
    P = plot(ps... , layout= @layout([a b; c d]), size = (800, 600))

    #savefig(P, "KitaevChain/kitaevfigs/kc_ent.pdf")
end