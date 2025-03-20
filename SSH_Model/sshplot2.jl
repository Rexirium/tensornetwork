using HDF5
using Plots, LinearAlgebra
using LaTeXStrings, Printf

P = let 
    file = h5open("SSH_Model/sshplotdata.h5","r")
    vs = read(file,"vs")
    densities = read(file,"densities")
    correlations = read(file,"correlations")
    close(file)
    L = length(densities[1,:])
    b=Int(L/2)
    cg1=cgrad(:matter, rev = true, scale=:exp10,categorical=false)
    p1=heatmap(vs,1:L,transpose(densities),xlabel=L"v",ylabel="site",
        framestyle=:box, title="density", c=cg1)
    cg2=cgrad(:roma, rev = false, scale=:exp10,categorical=false)
    p2=heatmap(vs,(1-b):(L-b),transpose(correlations),xlabel=L"v",ylabel="distance",
        framestyle=:box, title="correlation",c=cg2)
    lay=@layout([a b])
    plot(p1,p2,layout=lay,size=(1000,400))
end

savefig(P, "SSH_Model/sshfigs/ssh_density_corr.pdf")