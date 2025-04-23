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
    file1 = h5open("MajoranaChain/majoranasshbenchdata.h5", "r")
    vs = read(file1, "vs")
    entangles = read(file1, "entangles")
    close(file1)

    ns, L, nalg = size(entangles)
    L -= 1

    entangleerr = zeros(ns, L+1, nalg-1)
    for i=1:ns, a=1:(nalg - 1)
        if a !=3
            arr = abs.((entangles[i,:,a] .- entangles[i,:,3])./entangles[i,:,3])
        else
            arr = abs.((entangles[i,:,4] .- entangles[i,:,3])./entangles[i,:,3])
        end
        mask = arr .< 1.0E-16
        arr[mask] .= 1.0E-16
        entangleerr[i, :, a] .= arr
    end
    

    alglabel = ["MF(origin)" "MF(trans)" "CF" "JW" ]
    fewlabel = ["MF(origin)" "MF(trans)" "JW"]
    ylims =[(-0.001, 0.05), (-0.01, 0.55), (-0.01, 0.25) ]
    ps1 = []
    for i in 1:ns
        p = plot(entangles[i,:,:], xlabel="site j", ylim=ylims[i], lw=2, label=alglabel,
            title = latexstring("v = $(vs[i]),\\, w = 1.0"), leg =(i==ns))
        if i==1
            ylabel!(L"S(j)")
        end
        push!(ps1, p)
    end
    p1 = plot(ps1... , layout= @layout([a b c]), size = (800, 250))

    ps2 = []
    colors = palette(:default)[[1 2 4]]
    for i in 1:ns
        p = plot(entangleerr[i,:,:], xlabel="site j", ylim=(1.0E-16,1E-8), lw=2, label=fewlabel,
        yscale=:log10, c=colors, leg= (i==ns) )
        if i==1
            ylabel!(L"\epsilon")
        end
        push!(ps2, p)
    end
    p2 = plot(ps2... , layout= @layout([a b c]), size = (800, 250))
    P = plot(p1, p2, layout=@layout([a; b{0.45h}]), size=(800, 500))
    #savefig(P, "MajoranaChain/figs/mfssh_entangle.pdf")
end