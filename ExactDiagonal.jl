using LinearAlgebra
using BenchmarkTools
using ITensors, ITensorMPS

mutable struct NumState
    sitenum::Int
    stateind::Vector{<:Int}
    statevec::Vector

    NumState(sitenum, stateind, statevec) =  new(sitenum, stateind, statevec)
end

function readbit(num::Int, pos::Int)
    mask = 1<<(pos-1)
    return (num & mask)>>(pos-1)
end
# count number of ones between bit i and j (excluding i and j), i<j 
# and return the sign according to the odd or even of the number 
function signbetween(num::Int, i::Int, j::Int)
    mask = 1<<(j-i-1) -1
    segnum = (num>>i) & mask
    return (-1)^count_ones(segnum)
end

function numconserve_basis(Ls::Int, Ns::Int)
    Ns<=Ls || error("more fermions than sites!")
    basis = collect(range(0, 1<<Ls -1))
    basis = basis[count_ones.(basis).== Ns]
    return basis
end

function spectrum_BdG(H::AbstractMatrix; retstate::Bool=true)
    if retstate
        return eigen(H)
    else
        return eigvals(H)
    end
end

function spectrum_BdG(A::AbstractMatrix, B::AbstractMatrix; retstate::Bool=true)
    size(A) == size(B) || error("incompactible size of A and B")
    Ls = size(A, 1)
    H = [A B; -conj(B) -conj(A)]
    mu = tr(A)/Ls
    H = Hermitian(H)
    if retstate
        λ, T = eigen(H)
        spec = λ[1:Ls] .+ mu
        U = T[1:Ls, 1:Ls]
        V = T[(Ls+1): 2Ls, 1:Ls]
        return spec, U, V
    else
        spec = eigvals(H)[1:Ls] .+ mu
        return spec
    end
end

function single_excitation(H::AbstractMatrix, n::Int)
    Ls = size(H, 1)
    spec, U = eigen(H)
    stateind = 1 .<< range(0, Ls-1)
    statevec = U[:, n]
    return spec[n], NumState(Ls, stateind, statevec)
end

function numconserve_gs(U::AbstractMatrix, Ls::Int, Ns::Int)
    basis = numconserve_basis(Ls, Ns)
    dim = length(basis)
    statevec = zeros(ComplexF64, dim)
    for (i, num) in enumerate(basis)
        mask = BitArray(digits(num, base=2, pad=Ls))
        statevec[i] = det(U[mask, 1:Ns])
    end
    return NumState(Ls, basis, statevec)
end

function groundstate(H::AbstractMatrix; etol::Float64= 1.0E-14)
    Ls = size(H, 1)
    spec, U = eigen(H)
    neginds1 = findall(x -> x<= -etol, spec)
    neginds0 = findall(x -> x<=0.0, spec)
    neginds2 = findall(x -> x<= etol, spec)
    energy = sum(spec[neginds0])
    N1, N2 = length(neginds1), length(neginds2)
    deg = N2 - N1 + 1
    states = ntuple(x -> numconserve_gs(U, Ls, N1 + x-1), deg)
    return energy, states
end
# calculate the density on site i , cdag(i) c(i) , i<j
function expectation(state::NumState, i::Int)
    stind = state.stateind
    maski = readbit.(stind, i)
    stvec = state.statevec
    return real(dot(stvec, maski.*stvec))
end
# calculate the operator cdag(i) c(j) , i<j
function expectation(state::NumState, i::Int, j::Int)
    i < j || error("i is supposed to be smaller than j")
    stind = copy(state.stateind)
    stvec = copy(state.statevec)
    maski = readbit.(stind, i) .== 0
    maskj = readbit.(stind, j) .== 1
    maskold = maski .& maskj
    masknew = .~(maski .| maskj)
    fsign = signbetween.(stind[maskold], i, j)
    stvec[masknew].= stvec[maskold].*fsign
    stvec[.~masknew].= 0.0
    return dot(state.statevec, stvec)
end

function expectation(state::NumState, A::AbstractMatrix)
    Ls = state.sitenum
    Ls == size(A, 1) || error("Wrong matrix dimension of operator")
    resd = 0.0
    resu = 0.0+0.0im
    for i in 1:Ls
        resd += A[i,i] * expectation(state, i)
    end
    for i in 1:Ls-1
        for j in i+1:Ls
            resu += A[i,j]* expectation(state, i, j)
        end
    end
    return real(resd) + 2 * real(resu)
end

function density_vec(state::NumState)
    Ls = state.sitenum
    density = zeros(Ls)
    for j in 1:Ls
        density[j] = expectation(state, j)
    end
    return density
end

function correlation_mat(state::NumState)
    Ls = state.sitenum
    corr = zeros(ComplexF64, Ls, Ls)
    for i in 1:Ls
        for j in i:Ls
            if j == i
                corr[i, j] = expectation(state, i)
            else
                corr[i, j] = expectation(state, i, j)
            end
        end
    end
    return Matrix(Hermitian(corr))
end

function Hamiltonian(s::Vector{Index{Int}}, t1::Number, t2::Number) 
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        t= isodd(j) ? t1 : t2
        os += t, "Cdag",j, "C",j+1
        os += t', "Cdag",j+1, "C",j
    end
    return MPO(os, s)
end

        

v, w =2.0+1.0im, 1.0 -0.5im
L = 10
sw = Sweeps(10)
setmaxdim!(sw, 200)

arr = repeat([v, w], Int(L/2))[1:L-1]
H = Hermitian(diagm(1=>arr))

ss = siteinds("Fermion", L)
H_mpo = Hamiltonian(ss, v, w)
psi0 = random_mps(ss; linkdims = 4)

EED, psiED = groundstate(H)
EDMRG, psiDMRG = dmrg(H_mpo, psi0, sw; outputlevel=0)

energy = expectation(psiED[1], H)

EED, energy
