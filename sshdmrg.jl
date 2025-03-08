using HDF5
using Printf
include("sshmodel.jl")

#DMRG parameters
sw=Sweeps(15)
setmaxdim!(sw,200)
setcutoff!(sw,1E-14)
krydim=4
#=
let 
    L,D=40,5
    v,w=1.0, 2.0
    V1,V2=1.0,0.0
    b_ent=Int(L/2)
    obs=EntangleObserver(b_ent,[])
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    Hssh=SSH_obc(sites,v,w,[V1,V2])
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites SSH model for v = $v, w = $w and V1= $V1, V2 = $V2")
    energy, psi=dmrg(Hssh,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
    data= hcat((obs.data)...)
    plot(data[1,:],xlabel="sweep",ylabel=L"S_\mathrm{vN}",ylim=(0.66,0.72),label="entanglement entropy",color=:red,w=2,legend=(0.8,0.8))
    bar!(twinx(),data[2,:],ylabel=L"D",ylim=(0,100),label="maximum bond dimension",fillalpha=0.5,legend=(0.8,0.9))
end
=#
let
    L,D=40,5
    w,b=1.0,Int(L/2)
    num,num_ite=50,5
    ratio=Int(num/num_ite)
    vs=LinRange(0.0,2.5,num+1)
    vs_ite=LinRange(0.0,2.5,num_ite+1)
    sites=siteinds("Fermion",L;conserve_qns=false)
    psi0=random_mps(sites;linkdims=D)
    energies=zeros((num+1,2))
    entropies=zeros((num+1,2))
    densities=zeros((num+1,L))
    correlations=zeros((num+1,L))
    for (i,v) in enumerate(vs)
        energy_ED=SSH_ED(L,v,w)
        Hssh=SSH_obc(sites,v,w)
        energy_DMRG,psi_DMRG=dmrg(Hssh,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=0)
        entropy=entangle_entropy(psi_DMRG,b)
        ldim=maxlinkdim(psi_DMRG)
        density=expect(psi_DMRG,"N")
        nncorr=correlation_matrix(psi_DMRG,"N","N")[:,b]
        energies[i,:]=[energy_DMRG,energy_ED]
        entropies[i,:]=[entropy,ldim]
        densities[i,:]=density
        correlations[i,:]=nncorr
    end
    errors = energies[:,1]-energies[:,2]
    # ITE parameters
    cutoff=1.0E-12
    finaltemp=0.1
    tau=0.04
    steps=Int(ceil(1/finaltemp/tau))
    data_ITE=zeros((num_ite+1,4))
    for (i,v) in enumerate(vs_ite)
        Hssh=SSH_obc(sites,v,w)
        sshgate=SSH_gate(sites,v,w,tau)
        energy_ITE,psi_ITE=ImaginaryTimeEvolve(Hssh,sshgate,psi0,steps;cutoff=cutoff,display=false)
        ldim_ITE=maxlinkdim(psi_ITE)
        entropy_ITE=entangle_entropy(psi_ITE,b)
        err_ITE=energy_ITE-energies[ratio*(i-1)+1,2]
        entropy_ITE=entangle_entropy(psi_ITE,b)
        data_ITE[i,:]=[energy_ITE,err_ITE,entropy_ITE,ldim_ITE]
    end

    h5open("sshplotdata.h5","w") do file
        write(file,"vs",collect(vs))
        write(file,"vs_ite",collect(vs_ite))
        write(file,"energies", energies)
        write(file,"entropies", entropies)
        write(file,"errors",errors)
        write(file,"densities",densities)
        write(file,"correlations",correlations)
        write(file,"ITEdata",data_ITE)
    end
end
