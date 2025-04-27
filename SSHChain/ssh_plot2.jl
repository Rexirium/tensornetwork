using LinearAlgebra, MKL
using HDF5
using LaTeXStrings
using Plots
default(
    grid=true, 
    titlelocation=:left,
    framestyle=:box,
    guidefontsize=14,
    legendfontsize=10,
    tickfontsize=10
)

let 
    file = h5open("SSHChain/sshbenchdata2.h5", "r")
    vls = read(file, "vls")
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
        p = plot(xs, entangles[i,:,:], lw=1.5, label=alglabel,
            ylim=(0.0,1.0), title = latexstring("v = $(vls[i]),\\, w = 1.0"),
            legend_position=(i==5 ? :top : :bottom))
        if i==1 || i==3
            ylabel!(L"S(j)")
        end
        if i==3 || i==5
            xlabel!("cell "*(L"j"))
        end
        push!(ps, p)
    end
    P = plot(ps... , layout= @layout([a b; c d]), size = (800, 600))

    #savefig(P, "SSHChain/sshfigs/ssh_ent.pdf")
end