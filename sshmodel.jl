using ITensors, ITensorMPS
using LinearAlgebra

mutable struct EntangleObserver <: AbstractObserver
    bond::Int
    EntangleObserver(bond=1)=new(bond)
end
# construct evolution gates
function evolmatrix(arg::Float64)
    c,s =cosh(arg),sinh(arg)
    X = [c -s; -s c]
    A = [zeros(4);;[0.0 0.0]; X ; [0.0 0.0];;zeros(4)]
    A[1,1],A[4,4]=1.0,1.0
    return A
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

#construct Hamiltonian#
function SSH_OBC(s::Vector{Index{Int64}},t1::Float64,t2::Float64) 
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        t= isodd(j) ? t1 : t2
        os += t, "Cdag",j,   "C",j+1
        os += t, "Cdag",j+1, "C",j
    end
    Hamil=MPO(os, s)
    return Hamil
end
#inspect entanglement entropy after each sweep of DMRG
function ITensorMPS.measure!(O::EntangleObserver; psi,sweep_is_done,kwargs...)
    if sweep_is_done
        b=O.bond
        SvN=entangle_entropy(psi,b)
        println("  von Neumann SvN=$SvN")
    end
end
# evolution gates of ITE
function SSH_gate(s::Vector{Index{Int64}},t1::Float64,t2::Float64,dt::Float64)
    Ns=length(s)
    gates=ITensor[]
    A1=evolmatrix(t1*dt/2)
    A2=evolmatrix(t2*dt/2)
    for j in 1:Ns-1
        s1,s2=s[j],s[j+1]
        if isodd(j)
            U=ITensor(A1,(s1',s2',s1,s2))
        else
            U=ITensor(A2,(s1',s2',s1,s2))
        end
        push!(gates, U)
    end
    append!(gates,reverse(gates))
    return gates       
end
# imaginary time evolution algorithm
function ImaginaryTimeEvolve(Hamil::MPO, gates::Vector{ITensor} ,psi::MPS,evolsteps::Int;cutoff::Float64=1E-12, display::Bool=true)
    for p in 1:evolsteps
        psi=apply(gates,psi;cutoff)
        normalize!(psi)
        if (mod(p,10)==0)&&display 
            energy=inner(psi',Hamil, psi)
            H2=inner(Hamil,psi,Hamil, psi)
            varsq=H2-energy*energy
            println("After step $p  energy = $energy  variance = $varsq")
        end
    end
    energy=inner(psi',Hamil,psi)
    return energy, psi
end
# Exact diagonalization
function SSH_band(t1::Float64,t2::Float64,k::Float64)
    return sqrt(t1*t1+t2*t2+2*t1*t2*cos(k))
end
function SSH_spec(Lsize::Int,t1::Float64,t2::Float64)
    ks=range(-π,π,Int(Lsize/2))
    spectrum=SSH_band.(t1,t2,ks)
    sort(vcat(spectrum,-spectrum))
end

function SSH_spectrum(Lsize::Int,t1::Float64,t2::Float64;retstate::Bool=false)
    Ncell=Int(Lsize/2)
    arr=repeat([t1,t2],Ncell)
    deleteat!(arr,Lsize)
    Hamil=SymTridiagonal(zeros(Lsize),arr)
    if retstate==false
        return eigvals(Hamil)
    else
        return eigen(Hamil)
    end
end

function SSH_spectrum_pbc(Lsize::Int,t1::Float64,t2::Float64)
    Ncell=Int(Lsize/2)
    arr=repeat([t1,t2],Ncell)
    deleteat!(arr,Lsize)
    Hamil=diagm(1=>arr)
    Hamil[1,Lsize]=t2
    Hamil=Hermitian(Hamil)
    return eigvals(Hamil)
end

function SSH_ED(Lsize::Int,t1::Float64,t2::Float64,FermiLevel::Float64=0.0)
    spec=SSH_spectrum(Lsize,t1,t2)
    energy=sum(spec[spec[:].<=FermiLevel])
    return energy
end
#=
let
    L,D=40,5
    v,w=2.0,1.0
    vs=range(0.0,3.0,61)
    b_ent=Int(L/2)
    method="ite"
    #DMRG parameters
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    obs=EntangleObserver(b_ent)
    # ITE parameters
    cutoff=1.0E-12
    finaltemp=0.1
    tau=0.01
    steps=Int(1/finaltemp/tau)
    # initialize state
    state=[]
    for j in 1:L
        if isodd(j)
            push!(state,"1")
        else
            push!(state,"0")
        end
    end
    sites=siteinds("Fermion",L;conserve_qns=false)
    psi0=random_mps(sites;linkdims=D)
    Hssh=SSH_OBC(sites,v,w)
    sshgate=SSH_gate(sites,v,w,tau)

    SSH_ED(L,v,w)
    println("-------------------------------------------------------------")
    if method=="dmrg"
        println("Running DMRG for $L sites SSH model with v =$v, w=$w")
        energy,psi=dmrg(Hssh,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    elseif method=="ite"
        println("Running ITE for $L sites SSH model with v =$v, w=$w")
        energy,psi =ImaginaryTimeEvolve(Hssh,sshgate,psi0,steps,cutoff)
    end
    println("-------------------------------------------------------------")
 
    H2=inner(Hssh,psi,Hssh,psi)
    varsq=H2-energy*energy
    density=expect(psi,"N")
    SvN=entangle_entropy(psi,Int(L/2))
    println("ground state energy is $energy  variance of energy is $varsq")
    println("entanglement entropy across middle bond is $SvN \n")
    for (j,nc) in enumerate(density)
        println("density on site $j is \t $nc")
    end
    
    data=transpose(hcat(SSH_spectrum.(L,vs,w)...))
    plot(vs,data,w=1,leg=false,xlabel=L"v",framestyle=:box)
    yaxis!("energy",minorgrid=true)
    vline!([1.0],line=(1,:dash))
    annotate!([(0.5,-3,(L"w=1"))])
    savefig("figures/spectrum.pdf")

end

let 
    L,D=40,6
    w=1.0
    vs=range(0.0,3.0,11)
    b_ent=Int(L/2)
    #DMRG parameters
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    obs=EntangleObserver(b_ent)
    # ITE parameters
    cutoff=1.0E-12
    finaltemp=0.1
    tau=0.02
    steps=Int(1/finaltemp/tau)
    # initialize state
    state=[]
    for j in 1:L
        if isodd(j)
            push!(state,"1")
        else
            push!(state,"0")
        end
    end
    sites=siteinds("Fermion",L;conserve_qns=false)
    psi0=random_mps(sites;linkdims=D)

    energys=[]
    SvNs=[]
    for v in vs
        energy_ED=SSH_ED(L,v,w)
        Hssh=SSH_OBC(sites,v,w)
        sshgate=SSH_gate(sites,v,w,tau)
        energy_DMRG,psi_DMRG=dmrg(Hssh,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=0)
        energy_ITE,psi_ITE=ImaginaryTimeEvolve(Hssh,sshgate,psi0,steps;cutoff=cutoff,display=false)
        SvN_DMRG=entangle_entropy(psi_DMRG,b_ent)
        SvN_ITE=entangle_entropy(psi_ITE,b_ent)
        push!(energys,[energy_DMRG,energy_ITE,energy_ED])
        push!(SvNs,[SvN_DMRG,SvN_ITE])
    end
    energys=transpose(hcat(energys...))
    SvNs=transpose(hcat(SvNs...))
    plot(vs,energys,line=[:solid :solid :dash],lab=["DMRG" "ITE" "ED"])

end
=#

