include("../majorana.jl")

function majorana_hamiltonian(s::Vector{Index{Int}},t1::Real,t2::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += im*t1, "Gamma1", j, "Gamma2", j+1
        else
            os += im*t2, "Gamma2", j, "Gamma1", j+1
        end
    end
    return MPO(os, s)
end

function mfermion_hamiltonian(s::Vector{Index{Int}},t1::Real,t2::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        os += im*t1, "Gamma1", j, "Gamma2", j
        os += im*t2, "Gamma2", j, "Gamma1", j+1
    end
    os += im*t1, "Gamma1", Ns, "Gamma2", Ns
    return MPO(os,s)
end

function jordan_hamiltonian(s::Vector{Index{Int}},t1::Real,t2::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        os += -t1, "Z", j
        os += -t2, "X", j, "X", j+1
    end
    os += -t1, "Z", Ns
    return MPO(os,s)
end

function cfermion_hamiltonian(s::Vector{Index{Int}},t1::Real,t2::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns
        os += 2*t1, "Cdag", j, "C", j
        os += -t1, "Id", j
    end
    for j in 1:Ns-1
        os += -t2, "Cdag", j, "C", j+1
        os += -t2, "Cdag", j+1, "C", j
        os += t2, "C", j, "C", j+1
        os += t2, "Cdag", j+1, "Cdag", j
    end
    return MPO(os,s)
end

function spin_hamiltonian(s::Vector{Index{Int}},t1::Real,t2::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += t1, "X", j, "Y", j+1
        else
            os += t2, "Y", j, "X", j+1
        end
    end
    return MPO(os,s)
end

function isHermitian(Op::ITensor;tol=1.0e-12)
    Op2=apply(Op,Op)
    sqtr=scalar(Op2*delta(inds(Op)))
    norm2=inner(Op,Op)
    err=abs(sqtr-norm2)
    return err<tol
end

function isHermitian(A::MPO;tol=1.0e-12)
    Ns=length(A)
    A2=apply(A,A)
    V=ITensor(1.0)
    for j in 1:Ns
        ss=siteinds(A,j)
        V*=A2[j]*delta(ss)
    end
    sqtr=scalar(V)
    norm2=inner(A,A)
    err=abs(sqtr-norm2)
    return err<tol
end

let
    L,D=40,4
    v,w=2.0,1.0
    N=2*L 
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("MF",N)
    psi0=random_mps(sites;linkdims=D)
    H=majorana_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in Majorana basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------")
end

let 
    L,D=40,4
    v,w=2.0,1.0
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    H=mfermion_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in Fermion basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    L,D=40,4
    v,w=2.0,1.0
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("S=1/2",L)
    psi0=random_mps(sites;linkdims=D)
    H=jordan_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites Majorana Fermion chain in spin basis for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    L,D=40,4
    v,w=2.0,1.0
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("Fermion",L)
    psi0=random_mps(sites;linkdims=D)
    H=cfermion_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites simplified Kitaev chain for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
end

let 
    L,D=40,4
    v,w=2.0,1.0
    N=2*L 
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
    sites=siteinds("S=1/2",N)
    psi0=random_mps(sites;linkdims=D)
    H=spin_hamiltonian(sites,v,w)
    println("-----------------------------------------------------------------")
    println("Running DMRG for $L sites spin chain for v = $v, w = $w")
    energy, psi=dmrg(H,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------")
end