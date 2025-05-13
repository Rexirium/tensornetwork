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
    correlations = read(file, "correlations")
    close(file)

    ns, L, nalg = size(correlations)
    xs = range(1, L)

    alglabel = ["CF" "MF(trans)" "MF(origin)" "JW" ]
    ylims = [(0.2,0.6),(0.3,0.7),(0.6,1.0),(0.8,1.0),(0.8,1.0)]
    ps = []
    for i in 1:ns
        if i==4
            continue
        end
        p = plot(xs, correlations[i,:,:], lw=1.5, label=alglabel,
             ylim=ylims[i], title = latexstring("\\mu = $(muls[i]),\\, t= \\Delta = 1.0")
             )
        if i==1 || i==3
            ylabel!(L"C(20,j)")
        end
        if i==3 || i==5
            xlabel!("site "*(L"j"))
        end
        push!(ps, p)
    end
    P = plot(ps... , layout= @layout([a b; c d]), size = (800, 600), leftmargin=2Plots.mm)

    savefig(P, "KitaevChain/kitaevfigs/kc_corr.svg")
end