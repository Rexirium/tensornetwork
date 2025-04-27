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
    file1 = h5open("MajoranaChain/majoranasshbenchdata2.h5", "r")
    vs = read(file1, "vs")
    runtime = read(file1, "runtime")
    usespace = read(file1, "usespace")
    close(file1)

    file2 = h5open("MajoranaChain/majoranasshplotdata.h5", "r")
    vss = read(file2, "vs")
    maxbonds = read(file2, "maxbonds")
    close(file2)

    shifts = [-0.15 -0.05 0.05 0.15]
    vcoord = broadcast(+, vs, shifts)

    alglabel = ["MF(origin)" "MF(trans)" "CF" "JW" ]
    colors = palette(:default)[[1 2 3 4]]
    b1 = bar(vcoord, runtime, xlabel=L"v/w", ylabel= "runtime (s)", yscale=:log10, ylim=(0.1,1000),
        bar_width = 0.1, width=0.5, label=alglabel, legend_position=:topleft)
    plot!(twinx(), vss, maxbonds, ylim=(0,120), lw=1.5, label=alglabel, 
        c=colors, legend_position=:topright)
    
    b2 = bar(vcoord, usespace, xlabel=L"v/w", ylabel= "space (MB)", yscale=:log10, ylim=(100, 10^6),
        bar_width = 0.1, width=0.5, label=alglabel, legend_position=:topleft)
    plot!(twinx(), vss, maxbonds, ylim=(0,120), lw=1.5, label=alglabel, 
        c=colors, ylabel=L"D", legend_position=:topright)
    
    P = plot(b1, b2, layout=@layout([a b]), size=(700, 400))
    #savefig(P, "MajoranaChain/majoranafigs/mfssh_bench.pdf")
end