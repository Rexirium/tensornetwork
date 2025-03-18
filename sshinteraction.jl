using Plots
using Printf, LaTeXStrings
include("sshmodel.jl")

#DMRG parameters
sw=Sweeps(15)
setmaxdim!(sw, 200)
setcutoff!(sw, 1E-14)
krydim=4

let 
    L, D = 40,5
    v, w = 1.0, 2.0
    V, W = 1.0, 0.0
    b_ent=Int(L/2)
    obs=EntangleObserver(b_ent,[])
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    Hssh=SSH_obc(sites,v,w,[V,W])
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites SSH model for v = $v, w = $w and V1= $V, V2 = $W")
    energy, psi=dmrg(Hssh,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
    data= hcat((obs.data)...)
    plot(data[1,:],xlabel="sweep",ylabel=L"S_\mathrm{vN}",ylim=(0.66,0.72),label="entanglement entropy",color=:red,w=2,legend=(0.6,0.8))
    bar!(twinx(),data[2,:],ylabel=L"D",ylim=(0,100),label="maximum bond dimension",fillalpha=0.5,legend=(0.6,0.9))
end