 mutable struct EntangleObserver <: AbstractObserver
    bond::Int
    data::Vector{Any}
    EntangleObserver(bond=1,data=[])=new(bond,data)
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
        D=maxlinkdim(psi)
        push!(O.data,[SvN,D])
        #println("  von Neumann entropy SvN = $SvN, max link dimension D = $D")
    end
end