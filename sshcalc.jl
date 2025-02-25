using Plots, LaTeXStrings
include("sshmodel.jl")
#=
let 
    L,D=40,5
    v,w=1.0, 2.0
    V1,V2=1.0,0.0
    sw=Sweeps(15)
    setmaxdim!(sw,200)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    Hssh=SSH_obc(sites,v,w,[V1,V2])
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites SSH model for v = $v, w = $w and intracell interaction V1= $V1, intercell interaction V2 = $V2")
    energy, psi=dmrg(Hssh,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end
=#

p1=let 
    L=40
    v=1.0
    ws=LinRange(0,2.5,101)
    spec=SSH_spectrum_obc.(L,v,ws)
    data=transpose(hcat(spec...))
    plot(ws,data,w=1.5,leg=false,xlabel=L"w",ylabel=L"E",framestyle=:box)
    vline!([1.0],line=(1,:dash))
end
p2=let 
    L=40
    w=1.0
    vs=LinRange(0,2.5,101)
    spec=SSH_spectrum_obc.(L,vs,w)
    data=transpose(hcat(spec...))
    plot(vs,data,w=1.5,leg=false,xlabel=L"v",ylabel=L"E",framestyle=:box)
    vline!([1.0],line=(1,:dash))
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
end
layout=@layout([a b;c])
plot(p1,p2,p3,layout=layout,size=(800,600))

