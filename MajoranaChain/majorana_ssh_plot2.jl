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
    file1 = h5open("MajoranaChain/majoranasshbenchdata.h5", "r")
    vs = read(file1, "vs")
    correlations = read(file1, "correlations")
    entangles = read(file1, "entangles")
    runtime = read(file1, "runtime")
    usespace = read(file1, "usespace")
    close(file1)
    
    ns, L, nalg = size(correlations)
    correlationerr = zeros(ns, L, nalg-1)
    for i=1:ns, a=1:(nalg - 1)
        if a !=3
            arr = abs.((correlations[i,:,a] .- correlations[i,:,3])./correlations[i,:,3])
        else
            arr = abs.((correlations[i,:,4] .- correlations[i,:,3])./correlations[i,:,3])
        end
        mask = arr .< 1.0E-16
        arr[mask] .= 1.0E-16
        correlationerr[i,:,a] .= arr
    end

    alglabel = ["MF(origin)" "MF(trans)" "CF" "JW" ]
    fewlabel = ["MF(origin)" "MF(trans)" "JW"]

    ps1 = []
    for i in 1:ns
        p = plot(correlations[i,:,:], xlabel="site j", ylim = (-0.02,0.4), lw=2, label=alglabel,
            title = latexstring("v = $(vs[i]),\\, w = 1.0"))
        if i==1
            ylabel!(L"C(20,j)")
        end
        push!(ps1, p)
    end
    p1 = plot(ps1... , layout= @layout([a b c]), size = (800, 250))

    ps2 = []
    colors = palette(:default)[[1 2 4]]
    for i in 1:ns
        p = plot(correlationerr[i,:,:], xlabel="site j", ylim=(1.0E-15,2.5E-9), lw=2, label=fewlabel,
        yscale=:log10, c=colors, leg= (i==ns) )
        if i==1
            ylabel!(L"\epsilon")
        end
        push!(ps2, p)
    end
    p2 = plot(ps2... , layout= @layout([a b c]), size = (800, 250))
    P = plot(p1, p2, layout=@layout([a; b{0.45h}]), size=(800, 500))
    savefig(P, "MajoranaChain/majoranafigs/mfssh_corr.pdf")
end