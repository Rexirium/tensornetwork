using LinearAlgebra, MKL
using HDF5
using LaTeXStrings
using Plots
gr()
default(
    grid=true, 
    titlelocation=:left,
    framestyle=:box,
    guidefontsize=14,
    legendfontsize=10,
    tickfontsize=10
)

let 
    file1 = h5open("KitaevChain/kcbenchdata.h5", "r")
    muls = read(file1, "muls")
    runtime = read(file1, "runtime")
    usespace = read(file1, "usespace")
    close(file1)

    file2 = h5open("KitaevChain/kcplotdata.h5", "r")
    mus = read(file2, "mus")
    maxbonds = read(file2, "maxbonds")
    close(file2)

    ns, nalg = size(runtime)
    shifts = [-0.15 -0.05 0.05 0.15]
    mucoord = broadcast(+, muls, shifts)

    alglabel = ["CF" "MF(trans)" "MF(origin)" "JW" ]
    colors = palette(:default)[[1 2 3 4 5]]
    #=
    b1 = bar(mucoord, runtime, xlabel=L"\mu/t", ylabel= "runtime (s)", yscale=:log10, ylim=(0.1,100), 
        bar_width = 0.1, width=0.5, label=alglabel, legend_position=:topleft)

    plot!(twinx(), mus, maxbonds, lw=1.5, label=alglabel, yscale=:log10, ylim=(5,50),  
        yticks=([10, 20, 40],["10", "20", "40"]), ylabel=L"D", legend_position=:topright)
    #savefig(b1, "KitaevChain/kitaevfigs/kc_bench1.pdf")
   =#
    
    b2 = bar(mucoord, usespace, xlabel=L"\mu/t", ylabel= "space (MB)", yscale=:log10, ylim=(10^3, 10^5),
        bar_width = 0.1, width=0.5, label=alglabel, legend_position=:topleft)
    plot!(twinx(), mus, maxbonds, lw=1.5, label=alglabel, yscale=:log10, ylim=(5,50),
        yticks=([10, 20, 40],["10", "20", "40"]), ylabel=L"D", legend_position=:topright)
    
    savefig(b2, "KitaevChain/kitaevfigs/kc_bench2.pdf")
 
end