using LinearAlgebra
using BenchmarkTools

mutable struct NumState
    sitenum::Int
    stateind::Vector{<:Int}
    statevec::Vector{<:Number}

    NumState(sitenum, stateind, statevec) =  new(sitenum, stateind, statevec)
end

struct OpN end
struct OpC end

function readbit(num::Int, pos::Int)
    return (num >> (pos -1)) & 1 == 1
end
# count number of ones between bit i and j (excluding i and j), i<j 
# and return the sign according to the odd or even of the number 
function signbetween(num::Int, i::Int, j::Int)
    mask = 1<<(j-i-1) -1
    segnum = (num>>i) & mask
    return (-1)^count_ones(segnum)
end

function splitbasis(num::Int, b::Int)
    b >=0 || return 0, num
    left = num >> b
    right = num & ((1<<b) - 1)
    return right, left
end

function numconserve_basis(Ls::Int, Ns::Int)
    Ns > Ls && error("more fermions than sites!")
    Ns == 0 && return Int[0]
    basis = Int[]
    sizehint!(basis, binomial(Ls, Ns))
    maxind = (1<<Ls) - 1
    ind = (1<<Ns) - 1
    while ind <= maxind
        push!(basis, ind)
        u = ind & (-ind)
        v = ind + u
        next = v + ((v ⊻ ind) ÷ u) >> 2
        next > maxind && break
        ind = next
    end
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
    H = [A -conj(B); B -transpose(A)]
    H = Hermitian(H)
    if retstate
        λ, T = eigen(H)
        U = T[1:Ls, Ls+1 : 2Ls]
        V = T[(Ls+1): 2Ls, Ls+1 : 2Ls]
        return λ, U, V
    else
        spec = eigvals(H)
        return spec
    end
end

function matrixize(state::NumState, b::Int)
    basis = state.stateind
    mask = (1<<b) -1
    rinds = basis .>> b
    linds = basis .& mask
    rset = unique(rinds)
    lset = unique(linds)
    M, N = length(rset), length(lset)
    rdict = Dict(val=>idx for (idx, val) in enumerate(rset))
    ldict = Dict(val=>idx for (idx, val) in enumerate(lset))
    T = eltype(state.statevec)
    mat = zeros(T, M, N)
    for k in 1:length(state.stateind)
        i = rdict[rinds[k]]
        j = ldict[linds[k]]
        mat[i,j] = state.statevec[k]
    end
    return mat
end
# |ψ⟩ = ∑_ij m_ij | j i ⟩
function reduce_densitymat(state::NumState, b::Int; traceout = :righht)
    mat = matrixize(state, b)
    if traceout == :right
        return mat * mat'
    else
        return transpose(mat) * conj(mat)
    end
end

function entangle_entropy(state::NumState, b::Int)
    b <= 0 && return 0.0
    mat = matrixize(state, b)
    Σ = svd(mat).S
    SvN = 0.0
    rank = length(Σ)
    for n in 1:rank
        p = Σ[n]*Σ[n]
        if p > 0.0
            SvN -= p*log(p)
        end
    end
    return SvN
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
    T = eltype(U)
    statevec = zeros(T, dim)
    for (i, num) in enumerate(basis)
        mask = readbit.(num, 1:Ls)
        subU = @view U[mask, 1:Ns]
        statevec[i] = det(subU)
    end
    return NumState(Ls, basis, statevec)
end

function groundstate(H::AbstractMatrix; degtol::Float64= 1.0E-3)
    Ls = size(H, 1)
    spec, U = eigen(H)
    N0 = searchsortedlast(spec, 0.0)
    N1 = searchsortedlast(spec, -degtol) 
    N2 = searchsortedlast(spec, degtol)
    deg = N2 - N1 + 1
    energy = sum(spec[1:N0])
    states = Vector{NumState}(undef, deg)
    for i in 1:deg
        states[i] = numconserve_gs(U, Ls, N1 + i-1)
    end
    return energy, states
end

function groundstate_energy(H::AbstractMatrix)
    spec = eigvals(H)
    N0 = searchsortedlast(spec, 0.0)
    return sum(spec[1:N0])
end

function groundstate_energy(A::AbstractMatrix, B::AbstractMatrix)
    size(A) == size(B) || error("incompactible size of A and B")
    Ns = size(A, 1)
    H = [A -conj(B); B -transpose(A)]
    H = Hermitian(H)
    spec = eigvals(H, 1:Ns)
    return sum(spec)
end
# calculate the density on site i , cdag(i) c(i) , i<j
function expectation(state::NumState, i::Int)
    stind = state.stateind
    maski = readbit.(stind, i)
    stvec = state.statevec
    return real(dot(stvec, maski.*stvec))
end
# calculate the operator cdag(i) c(j) , i<j
function expectation(state::NumState, i::Int, j::Int, ::Val{OpC})
    i == j && return expectation(state, i)
    i, j = min(i,j), max(i,j)
    stind = copy(state.stateind)
    stvec = copy(state.statevec)
    maski = .~readbit.(stind, i)
    maskj = readbit.(stind, j) 
    maskold = maski .& maskj
    masknew = .~(maski .| maskj)
    fsign = signbetween.(stind[maskold], i, j)
    stvec[masknew].= stvec[maskold].*fsign
    stvec[.~masknew].= 0.0
    return dot(state.statevec, stvec)
end

function expectation(state::NumState, i::Int, j::Int, ::Val{OpN})
    i ==j && return expect(state, i)
    stind = state.stateind
    maski = readbit.(stind, i)
    maskj = readbit.(stind, j)
    stvec = state.statevec
    return real(dot(maskj.*stvec, maski.*stvec))
end

function expectation(state::NumState, A::AbstractMatrix)
    Ls = state.sitenum
    Ls == size(A, 1) || error("Wrong matrix dimension of operator")
    resd = 0.0
    resu = complex(0.0, 0.0)
    for i in 1:Ls
        resd += A[i,i] * expectation(state, i)
    end
    for i in 1:Ls-1
        for j in i+1:Ls
            resu += A[i,j]* expectation(state, i, j, Val(OpC))
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

function correlation_mat(state::NumState, val)
    Ls = state.sitenum
    corr = zeros(eltype(state.statevec), Ls, Ls)
    for i in 1:Ls
        for j in i:Ls
            if j == i
                corr[i, j] = expectation(state, i)
            else
                corr[i, j] = expectation(state, i, j, val)
            end
        end
    end
    return Matrix(Hermitian(corr))
end


