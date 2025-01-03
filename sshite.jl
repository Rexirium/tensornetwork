using ITensors, ITensorMPS

let
    L,D=40,5
    v,w=2.0,1.0
    b_ent=Int(L/2)
    cutoff=1.0E-12
    finaltemp=0.1
    tau=0.01
    steps=Int(1/finaltemp/tau)
    
    function evolmatrix(arg::Float64)
        c,s =cosh(arg),sinh(arg)
        X = [c -s; -s c]
        A = [zeros(4);;[0.0 0.0]; X ; [0.0 0.0];;zeros(4)]
        A[1,1],A[4,4]=1.0,1.0
        return A
    end
    function SSH_gate(s::Vector{Index{Int64}},t1::Float64,t2::Float64,dt::Float64)
        Ns=length(s)
        gates=ITensor[]
        os=OpSum()
        A1=evolmatrix(t1*dt/2)
        A2=evolmatrix(t2*dt/2)
        for j in 1:Ns-1
            s1,s2=s[j],s[j+1]
            if isodd(j)
                os += t1, "Cdag",j,   "C",j+1
                os += t1, "Cdag",j+1, "C",j
                U=ITensor(A1,(s1',s2',s1,s2))
            else
                os += t2, "Cdag",j,   "C",j+1
                os += t2, "Cdag",j+1, "C",j
                U=ITensor(A2,(s1',s2',s1,s2))
            end
            push!(gates, U)
        end
        Hamil = MPO(os,s)
        append!(gates,reverse(gates))
        return Hamil, gates       
    end

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

    function ImaginaryTimeEvolve(Hamil::MPO, gates::Vector{ITensor} ,psi::MPS,evolsteps::Int,cutoff::Float64)
        for p in 1:evolsteps
            psi=apply(gates,psi;cutoff)
            normalize!(psi)
            if mod(p,10)==0
                energy=inner(psi',Hamil, psi)
                H2=inner(Hamil,psi,Hamil, psi)
                varsq=H2-energy*energy
                SvN=entangle_entropy(psi,b_ent)
                println("After step $p  energy = $energy  variance = $varsq")
                println("    von Neumann entropy is $SvN")
            end
        end
        return energy,varsq, psi
    end
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
    Hssh,sshgate=SSH_gate(sites,v,w,tau)
    psi0=random_mps(sites,state;linkdims=D)
    
    println("----------------------------------------------------------")
    println("Running ITE for $L sites SSH model with v =$v, w=$w")
    energy,varsq,psi =ImaginaryTimeEvolve(Hssh,sshgate,psi0,steps,cutoff)
    println("----------------------------------------------------------")
    density=expect(psi,"N")
    SvN=entangle_entropy(psi,Int(L/2))
    println("ground state energy is $energy  variance is $varsq")
    println("entanglement entropy across middle bond is $SvN \n")
    for (j,nc) in enumerate(density)
        println("density on site $j is \t$nc")
    end
end