using HDF5
using Plots, LinearAlgebra
using LaTeXStrings, Printf

let 
    file = h5open("sshplotdata.h5","r")
    vs = read(file,"vs")
    vs_ite = read(file, "vs_ite")
    energies = read(file, "energies")
    errors = read(file, "errors")
    ITEdata = read(file, "ITEdata")
    entropies = read(file, "entropies")
    close(file)

    p1=plot(vs,energies,line=[:solid :dash],label=["DMRG" "ED"], ylabel=L"E",w=2,framestyle=:box)
    scatter!(vs_ite,ITEdata[:,1],label="ITE")
    p2=plot(vs,errors,framestyle=:box,w=2,line=:dashdot,xlabel=L"v",ylabel=L"\epsilon",yscale=:log10,label="DMRG err")
    scatter!(vs_ite,ITEdata[:,2],label="ITE err",yscale=:log10)
    p3=plot(vs, entropies[:,1], label="DMRG entropy", color=:red,ylabel=L"S_\mathrm{vN}",w=2,framestyle=:box)
    scatter!(vs_ite,ITEdata[:,3],label="ITE entropy")
    p4=plot(vs,entropies[:,2],xlabel=L"v",ylabel=L"D",ylim=(0,220),label="DMRG bond dim",w=2,c=:green,framestyle=:box)
    scatter!(vs_ite,ITEdata[:,4],label="ITE bond dim",ylim=(0,220))
    
    lay=@layout([a b ; c d])
    plot(p1,p3,p2,p4,layout=lay,size=(800,600))
end