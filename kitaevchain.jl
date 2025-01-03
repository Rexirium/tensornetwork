using ITensors, ITensorMPS
using LinearAlgebra

mutable struct EntangleObserver <: AbstractObserver
    bond::Int
    EntangleObserver(bond=1)=new(bond)
end
# calculate entanglement entropy
function entangle_entropy(psi::MPS, b::Int)
    psi_tmp=orthogonalize(psi,b)
    llink = linkinds(psi_tmp,b-1)
    lsite = siteinds(psi_tmp,b)
    U,S,V = svd(psi_tmp[b],(llink...,lsite...))
    SvN=0.0
    for n=1:dim(S,1)
        p = abs2(S[n,n])
        SvN-= p*log(p)
    end
    return SvN
end
#inspect entanglement entropy after each sweep of DMRG
function ITensorMPS.measure!(O::EntangleObserver; psi,sweep_is_done,kwargs...)
    if sweep_is_done
        b=O.bond
        SvN=entangle_entropy(psi,b)
        println("  von Neumann SvN=$SvN")
    end
end

ITensors.op(::OpName"Beta",::SiteType"Fermion")=[0 1 ; 1 0]
ITensors.op(::OpName"Gamma",::SiteType"Fermion")=[0 -im ; im 0]
ITensors.op(::OpName"Id",::SiteType"Fermion")=[1 0; 0 1]
    
ITensors.has_fermion_string(::OpName"Beta",::SiteType"Fermion")=true
ITensors.has_fermion_string(::OpName"Gamma",::SiteType"Fermion")=true
ITensors.has_fermion_string(::OpName"Id",::SiteType"Fermion")=false

function KitaevChainSpectrum(Ls::Int,mu::Float64,tt::Float64,delta::Number)
    Ns=2*Ls
    ev1=repeat([mu,tt+delta],Ls)[1:Ns-1]
    ev3=repeat([0,-tt+delta],Ls-1)[1:Ns-3]
    Ham=diagm(1=>ev1,3=>ev3)
    Ham=Hermitian(Ham)
    return eigvals(Ham)
end

function KitaevChainED(Ls::Int,mu::Float64,tt::Float64,delta::Number)
    spec=KitaevChainSpectrum(Ls,mu,tt,delta)
    energy=sum(spec[spec[:].<=mu])
    return energy
end

function KitaevChain(s::Vector{Index{Int64}}, mu::Float64,tt::Float64,delta::Number)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns
        os+= -mu, "N", j
        os+= mu/2, "Id", j
    end  
    for j in 1:Ns-1
        os+= -tt, "Cdag", j, "C", j+1
        os+= -tt, "Cdag", j+1, "C", j
        os+= delta, "C",j, "C", j+1
        os+= delta',"Cdag",j+1, "Cdag", j
    end      
    Hamil=MPO(os,s)
    return Hamil
end
function KitaevChainMF(s::Vector{Index{Int64}}, mu::Float64, tt::Float64, delta::Float64)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns
        os+=-im*mu/2, "Beta", j, "Gamma", j
    end
    for j in 1:Ns-1
        os+=im*(tt+delta)/2, "Gamma", j, "Beta", j+1
        os+=im*(-tt+delta)/2, "Beta", j, "Gamma", j+1
    end
    Hamil=MPO(os,s)
    return Hamil
end
#=
let
    L,D=40,6
    μ, t, Δ =1.0,2.0,0.5
    b_ent=Int(L/2)

    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    obs=EntangleObserver(b_ent)
    
    sites=siteinds("Fermion",L;conserve_qns=false)
    Hkc=KitaevChainMF(sites,μ,t,Δ)
    psi0=random_mps(sites;linkdims=D)

    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Kitaev chain with μ =$μ, t=$t, Δ=$Δ")
    energy, psi=dmrg(Hkc,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------")

    H2=inner(Hkc,psi,Hkc,psi)
    varsq=H2-energy*energy
    density=expect(psi,"N")
    SvN=entangle_entropy(psi,Int(L/2))
    println("ground state energy is $energy  variance of energy is $varsq")
    println("entanglement entropy across middle bond is $SvN \n")
    for (j,nc) in enumerate(density)
        println("density on site $j is \t $nc")
    end
end  sso
=#