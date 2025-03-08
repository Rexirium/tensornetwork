using HDF5
using Plots, LinearAlgebra
using LaTeXStrings, Printf

let 
    file = h5open("sshplotdata.h5","r")
    vs = read(file,"vs")
    densities = read(file,"densities")
    correlations = read(file,"correlations")
    close(file)
    L = length(densities[1,:])
    b=Int(L/2)
    cg1=cgrad(:matter, rev = true, scale=:exp10,categorical=false)
    p1=heatmap(vs,1:L,transpose(densities),xlabel=L"v",ylabel="site",
    framestyle=:box,colorbar_title="density",colorbar=:top,c=cg1)
    cg2=cgrad(:thermal, rev = false, scale=:exp10,categorical=false)
    p2=heatmap(vs,(1-b):(L-b),transpose(correlations),xlabel=L"v",ylabel="distance",
    framestyle=:box,colorbar_title="correlation",c=cg2)
    lay=@layout([a b])
    plot(p1,p2,layout=lay,size=(800,400))
end