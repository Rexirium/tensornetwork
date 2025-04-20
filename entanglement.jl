mutable struct EntangleObserver <: AbstractObserver
    bond::Int
    data::Vector{Any}
    EntangleObserver(bond=1, data=[])=new(bond,data)
end

mutable struct BondsObserver <: AbstractObserver
    bonds::Vector{Int}
    BondsObserver(bonds = Int[]) = new(bonds)
end

# calculate entanglement entropy
function entangle_entropy(psi::MPS, b::Int)
    b <= 0 && return 0.0
    psi_tmp = orthogonalize(psi,b)
    llink = linkinds(psi_tmp,b-1)
    lsite = siteinds(psi_tmp,b)
    U,S,V = svd(psi_tmp[b], (llink..., lsite...))
    SvN = 0.0
    for n = 1:dim(S,1)
        p = S[n, n]*S[n, n]
        SvN -= p*log(p)
    end
    return SvN
end

function entangle_entropy(state::NumState, b::Int)
    b <= 0 && return 0.0
    mat = matrixize(state, b)
    Σ = svd(mat).S
    SvN = 0.0
    rank = length(Σ)
    for n in 1:rank
        p = Σ[n]*Σ[n]
        SvN -= p*log(p)
    end
    return SvN
end
#inspect entanglement entropy after each sweep of DMRG
function ITensorMPS.measure!(O::EntangleObserver; psi, sweep_is_done, kwargs...)
    if sweep_is_done
        b = O.bond
        SvN = entangle_entropy(psi, b)
        D = maxlinkdim(psi)
        push!(O.data, [SvN, D])
        #println("  von Neumann entropy SvN = $SvN, max link dimension D = $D")
    end
end

function ITensorMPS.measure!(O::BondsObserver; psi, sweep_is_done, kwargs...)
    if sweep_is_done
        D = maxlinkdim(psi)
        push!(O.bonds, D)
    end
end