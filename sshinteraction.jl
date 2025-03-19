using HDF5
include("sshmodel.jl")

#DMRG parameters
sw=Sweeps(15)
setmaxdim!(sw, 200)
setcutoff!(sw, 1E-14)
krydim=4
#=
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
    println("Running DMRG for $L sites SSH model for v = $v, w = $w and V= $V, W = $W")
    energy, psi=dmrg(Hssh,psi0,sw;observer=obs,eigsolve_krylovdim=krydim,outputlevel=1)
    println("-----------------------------------------------------------------") 
    data= hcat((obs.data)...)
    plot(data[1,:],xlabel="sweep",ylabel=L"S_\mathrm{vN}",
        ylim=(0.66,0.72),label="entanglement entropy",color=:red,w=2,legend=(0.6,0.8))
    bar!(twinx(),data[2,:],ylabel=L"D",ylim=(0,100), 
        label="maximum bond dimension",fillalpha=0.5,legend=(0.6,0.9))
end
=#

let 
    L, D = 40, 5
    v, w = 2.0, 1.0
    num = 10
    b = Int(L/2)
    Vs = LinRange(0.0, 2.0, num+1)
    Ws = LinRange(0.0, 2.0, num+1)

    sites = siteinds("Fermion", L; conserve_qns=false)
    psi0 = random_mps(sites; linkdims=D)

    energies = zeros((num+1, num+1))
    energytag = "energies"*"v$v"*"w$w"
    entropies = zeros((num+1, num+1))
    entropytag = "entropies"*"v$v"*"w$w"

    for (i, V) in enumerate(Vs)
        for (j, W) in enumerate(Ws)
            Hint = SSH_obc(sites, v, w, [V*v, W*w])
            energy, psi = dmrg(Hint, psi0, sw; outputlevel=0)
            entropy = entangle_entropy(psi, b)
            energies[i,j] = energy
            entropies[i,j] = entropy
        end
    end
    
    h5open("sshintdata.h5", "cw") do file
        write(file, energytag, energies)
        write(file, entropytag, entropies)
    end
end