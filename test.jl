using ITensors, ITensorMPS
using LinearAlgebra
using Plots, LaTeXStrings

include("ExactDiagonal.jl")
include("entanglement.jl")

function makeHamiltonian(s::Vector{Index{Int}}, hamil::AbstractMatrix)
    length(s) == size(hamil, 1) || error("Wrong length of MPS indices")
    Ns = length(s)
    os = OpSum()
    for i in 1:Ns, j in i:Ns
        if abs(hamil[i, j]) <= 1.0E-16
            continue
        end
        if i == j
            os += hamil[i, i], "N", i
        else
            os += hamil[i, j], "Cdag", i, "C", j
            os += hamil[i, j]', "Cdag", j, "C", i
        end
    end
    return MPO(os, s)
end

function makeSSH(Ls::Int, t1::Number, t2::Number)
    arr = repeat([t1, t2], Int(Ls/2))[1:Ls-1]
    H = diagm(1=>arr, -1=>conj(arr))
    return Hermitian(H)
end
function makeSSH(s::Vector{Index{Int}}, t1::Number, t2::Number)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        t= isodd(j) ? t1 : t2
        os += t, "Cdag",j, "C",j+1
        os += t', "Cdag",j+1, "C",j
    end
    return MPO(os, s)
end

function makeKitaevChain(Ls::Int, mu::Real, tt::Number, delta::Number)
    arrt = -tt .* ones(Ls-1)
    arrtc = -tt' .* ones(Ls-1)
    arrd = delta .* ones(Ls-1)
    A = diagm(0=>(-mu .* ones(Ls)), 1=>arrt, -1=>arrtc)
    B = diagm(1=>arrd, -1=>(-arrd))
    return A, B
end

function makeKitaevChain(s::Vector{Index{Int}}, mu::Real, tt::Number, delta::Number)
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

let 
    # model setup
    L, D = 20, 4
    v, w = 0.5, 1.0
    bs = range(0, L)
    center = L÷2
    # DMRG parameters and initialization
    sw =Sweeps(10)
    setmaxdim!(sw, 200)
    setcutoff!(sw, 1.0E-14)
    krydim = 4
    obs = BondsObserver()
    sites = siteinds("Fermion", L)
    psi0 = random_mps(sites; linkdims = D)
    # define hamiltonian and corresponding MPO
    mu =  repeat([1.0, -1.0], center)
    tnn = ones(L-1)
    tnnn = - 0.0 .* ones(L-2)
    H = Hermitian(diagm(0=>mu, 1=>tnn, 2=>tnnn))
    H_MPO = makeHamiltonian(sites, H)

    energy_DMRG, psi_DMRG = dmrg(H_MPO, psi0, sw, observer = obs; eigsolve_krylovdim = krydim)
    energy_ED, psi_ED = groundstate(H)
    maxbond = max((obs.bonds)...)
    spectrum = spectrum_BdG(H, retstate = false)

    density_DMRG = expect(psi_DMRG, "N")
    density_ED = density_vec(psi_ED[1])

    corr_DMRG = real.(correlation_matrix(psi_DMRG, "Cdag", "C")[center, :])
    corr_ED = real.(correlation_mat(psi_ED[1], Val(OpC))[center, :])


    entangle_DMRG = zeros(L+1)
    entangle_ED = zeros(L+1)
    for (i, b) in enumerate(bs)
        entangle_DMRG[i] = entangle_entropy(psi_DMRG, b)
        entangle_ED[i] = entangle_entropy(psi_ED[1], b)
    end
    
    println("DMRG energy is $energy_DMRG, ED energy is $energy_ED")
    println("maxbond during DMRG is $maxbond")
    hm = heatmap(Matrix(H), yflip = true, xtick=false, ytick=false, title="hamiltonian")
    sc = scatter(spectrum, xlabel = L"n", ylabel="energy", framestyle = :box, legend = false, title="spectrum")
    b1 = bar(density_DMRG, ylabel = "density", ylim=(0.0, 1.0),
        legend = false, framestyle = :box, color = :red, title = "DMRG density")
    b2 = bar(density_ED, xlabel = "site j", ylabel = "density", ylim=(0.0, 1.0),
        legend = false, framestyle = :box, color = :blue, title = "ED density")
    pl = plot(bs, [entangle_DMRG, entangle_ED], xlabel="site", xlim=(0, L), ylabel="SvN", framestyle=:box,
        label=["DMRG" "ED"], lw=2, title="entanglement entropy")
    pc = plot([corr_DMRG, corr_ED], xlabel="site j", xlim=(1, L), ylabel="C($center, j)", framestyle=:box,
        label=["DMRG" "ED"], lw=2, title="correlation")
    plot(hm, sc, b1, b2, pl, pc, layout = @layout([a{0.4h} b{0.4h}; [c; d] [e; f]]), size = (1000, 800))
end


