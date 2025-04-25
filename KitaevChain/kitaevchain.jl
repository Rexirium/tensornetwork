include("../MajoranaRep.jl")

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

function KitaevChainSpectrum(Ls::Int,mu::Real,tt::Real,delta::Number)
    Ns=2*Ls
    ev1=repeat([mu,tt+delta],Ls)[1:Ns-1]
    ev3=repeat([0,-tt+delta],Ls-1)[1:Ns-3]
    Ham=diagm(1=>ev1,3=>ev3)
    Ham=Hermitian(Ham)
    return eigvals(Ham)
end

function KitaevChainED(Ls::Int,mu::Real,tt::Real,delta::Number)
    spec = KitaevChainSpectrum(Ls,mu,tt,delta)
    energy = sum(spec[spec[:].<=0.0])
    return energy
end

function KitaevChain(s::Vector{<:Index}, mu::Real,tt::Real,delta::Number)
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
    return MPO(os,s)
end

function KitaevChain(s::Vector{<:Index}, mu::Real, tt::Real, delta::Number, V::Real)
    Ns = length(s)
    os = OpSum()
    for j in 1:Ns
        os += -mu, "N", j
        os += mu/2, "Id", j
    end
    for j in 1:Ns-1
        os += -tt, "Cdag", j, "C", j+1
        os += -tt, "Cdag", j+1, "C", j
        os += delta, "C", j, "C", j+1
        os += delta', "Cdag", j+1, "Cdag", j
        os += V, "N", j, "N", j+1
    end
    return MPO(os,s)
end

function KitaevChainMF(s::Vector{<:Index}, mu::Real, tt::Real, delta::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns
        os+=-im*mu/2, "Gamma1", j, "Gamma2", j
    end
    for j in 1:Ns-1
        os+=im*(tt+delta)/2, "Gamma2", j, "Gamma1", j+1
        os+=im*(-tt+delta)/2, "Gamma1", j, "Gamma2", j+1
    end
    return MPO(os,s)
end

function KitaevChainMF(s::Vector{<:Index}, mu::Real,tt::Real,delta::Real, V::Real)
    Ns=length(s)
    os=OpSum()
    os += V/4, "Id", 1
    for j in 1:Ns
        if j==1 || j==Ns
            os += im*(-mu/2+V/4), "Gamma1", j, "Gamma2", j
        else
            os += im*(-mu+V)/2, "Gamma1", j, "Gamma2", j
            os += V/4, "Id", j
        end
    end
    for j in 1:Ns-1
        os += V/4, "Id", j
        os += im*(tt+delta)/2, "Gamma2", j, "Gamma1", j+1
        os += im*(-tt+delta)/2, "Gamma1", j, "Gamma2", j+1
        os += V/4, "Gamma1" ,j, "Gamma1" ,j+1, "Gamma2" ,j, "Gamma2", j+1
    end
    return MPO(os,s)
end

let
    L,D=40,5
    μ, t, Δ = -4.0, 1.0, 1.0
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
    println("Running DMRG for $L sites Kitaev chain with μ = $μ, t = $t, Δ = $Δ")
    energy, psi=dmrg(Hkc,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------")
    #=
    H2=inner(Hkc,psi,Hkc,psi)
    varsq=H2-energy*energy
    density=expect(psi,"N")
    SvN=entangle_entropy(psi,Int(L/2))
    println("ground state energy is $energy  variance of energy is $varsq")
    println("entanglement entropy across middle bond is $SvN \n")
    for (j,nc) in enumerate(density)
        println("density on site $j is \t $nc")
    end
    =#
end  
