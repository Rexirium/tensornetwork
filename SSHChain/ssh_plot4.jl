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
    file1 = h5open("SSHChain/sshbenchdata2.h5", "r")
    vls = read(file1, "vls")
    runtime = read(file1, "runtime")
    usespace = read(file1, "usespace")
    close(file1)

    file2 = h5open("SSHChain/sshplotdata2.h5", "r")
    vs = read(file2, "vs")
    maxbonds = read(file2, "maxbonds")
    close(file2)

    ns, nalg = size(runtime)
    shifts = [-0.1 -0.05 0.0 0.05 0.1]
    vcoord = broadcast(+, vls, shifts)

    alglabel = ["CF" "MF(trans)" "MF(origin)" "MF(simpl)" "JW" ]
    colors = palette(:default)[[1 2 3 4 5]]
    
    b1 = bar(vcoord, runtime, xlabel=L"v/w", ylabel= "runtime (s)", yscale=:log10, ylim=(1,200),
        bar_width = 0.05, width=0.5, label=alglabel, legend_position=:topleft)

    plot!(twinx(), vs, maxbonds, ylim=(1,200), lw=1.5, label=alglabel, yscale=:log10,
        ylabel=L"D", legend_position=:topright)
    #savefig(b1, "SSHChain/sshfigs/ssh_bench1.pdf")
   


    #=
    b2 = bar(vcoord, usespace, xlabel=L"v/w", ylabel= "space (MB)", yscale=:log10, ylim=(10^3, 10^6),
        bar_width = 0.05, width=0.5, label=alglabel, legend_position=:topleft)
    plot!(twinx(), vs, maxbonds, lw=1.5, label=alglabel, yscale=:log10, ylim=(1,200),
         ylabel=L"D", legend_position=:topright)
    
    #savefig(b2, "SSHChain/sshfigs/ssh_bench2.pdf")
   =#
end