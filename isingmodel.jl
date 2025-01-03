using ITensors, ITensorMPS

let 
    L,D=100,8
    J1,J2=2.0,1.0

    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-12)
    krydim=4

    function Ising_OBC(s::Vector{Index{Int64}},Jn::Float64, Jnn::Float64)
        Ns=length(s)
        os=OpSum()
        for j in 1:Ns-1
            os+=Jn, "Sz", j, "Sz", j+1
            if j<Ns-1
                os+=Jnn, "Sz", j, "Sz", j+2
            end
        end
        Hamil=MPO(os,s)
        return Hamil
    end

    sites=siteinds("S=1/2",L;conserve_qns=false)
    H_ising=Ising_OBC(sites,J1,J2)
    psi0=random_mps(sites;linkdims=D)

    println("-------------------------------------------------------------")
    println("Running DMRG for $L sites Ising spin chain with n.n. interaction $J1, n.n.n. interaction $J2")
    energy, psi=dmrg(H_ising,psi0,sw;eigsolve_krylovdim=krydim,outputlevel=1)
    println("-------------------------------------------------------------")
end