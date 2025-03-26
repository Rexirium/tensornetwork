using Plots
using LaTeXStrings,Printf
include("sshmodel.jl")

p1=let 
    L=40
    v=1.0
    ws=LinRange(0,2.5,101)
    spec=SSH_spectrum_obc.(L,v,ws)
    data=transpose(hcat(spec...))
    plot(ws,data,w=1.5,leg=false,xlabel=L"w",ylabel=L"E",framestyle=:box)
    vline!([1.0],line=(1,:dash))
    annotate!([(0.5,3,(L"v=1"))])
end
p2=let 
    L=40
    w=1.0
    vs=LinRange(0,2.5,101)
    spec=SSH_spectrum_obc.(L,vs,w)
    data=transpose(hcat(spec...))
    plot(vs,data,w=1.5,leg=false,xlabel=L"v",ylabel=L"E",framestyle=:box)
    vline!([1.0],line=(1,:dash))
    annotate!([(0.5,3,(L"w=1"))])
end
p3=let 
    L=40
    w=1.0
    th=LinRange(0.0,1.0,101)
    vx=1.5.-cos.(π*th)
    ux=sin.(π*th)
    spec=SSH_spectrum_obc.(L,vx,w,ux)
    data=transpose(hcat(spec...))
    plot(th,data,w=1.5,leg=false,xlabel=L"\theta",ylabel=L"E",framestyle=:box)
    vline!([1/3],line=(1,:dash))
    annotate!([(0.15,3,(L"w=1"))])
end
p4=let
    L=40
    v,w=0.5,1.0
    ls=[35,20,5]
    lay=@layout([a;b;c])
    spec,states=SSH_spectrum_obc(L,v,w;retstate=true)
    wavefunc=[]
    labels=[]
    for (i,l) in enumerate(ls)
        push!(wavefunc,states[:,l])
        push!(labels,@sprintf("ϵ = %.2f",spec[l]))
    end
    labels=reshape(labels,(1,3))
    plot(wavefunc,xlabel="site "*(L"n"),ylabel=L"\psi",ylim=(-0.4,0.6),framestyle=:box,w=1.5,label=labels)
    annotate!([(15,0.5,(L"v=0.5,w=1"))])
end

layout=@layout([a b;c d])
P = plot(p1,p2,p3,p4,layout=layout,size=(900,700))
#savefig(P, "SSH_Model/sshfigs/ssh_spectrum.pdf")

