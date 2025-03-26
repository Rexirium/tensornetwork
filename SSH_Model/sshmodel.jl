using ITensors, ITensorMPS
using LinearAlgebra
include("../entanglement.jl")

# construct evolution gates
function evolmatrix(arg::Float64)
    c,s =cosh(arg),sinh(arg)
    X = [c -s; -s c]
    A = [zeros(4);;[0.0 0.0]; X ; [0.0 0.0];;zeros(4)]
    A[1,1],A[4,4]=1.0,1.0
    return A
end

#construct Hamiltonian#
function SSH_obc(s::Vector{Index{Int}},t1::Number,t2::Number) 
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        t= isodd(j) ? t1 : t2
        os += t, "Cdag",j, "C",j+1
        os += t', "Cdag",j+1, "C",j
    end
    return MPO(os, s)
end

function SSH_obc(s::Vector{Index{Int}},t1::Number,t2::Number, u::Real)
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += u, "N", j
            os += t1, "Cdag", j, "C", j+1
            os += t1', "Cdag", j+1, "C", j
        else
            os += -u, "N", j
            os += t2, "Cdag", j, "C", j+1
            os += t2', "Cdag", j+1, "C", j
        end
    end
    os += -u, "N", Ns
    return MPO(os,s)
end

function SSH_obc(s::Vector{Index{Int}},t1::Number,t2::Number, V::Vector{<:Real})
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += t1, "Cdag", j, "C", j+1
            os += t1', "Cdag", j+1, "C", j
            os += V[1], "N", j, "N", j+1
        else
            os += t2, "Cdag", j, "C", j+1
            os += t2', "Cdag", j+1, "C", j
            os += V[2], "N", j, "N", j+1
        end
    end
    return MPO(os,s)
end

function SSH_obc(s::Vector{Index{Int}},t1::Number,t2::Number,u::Real,V::Vector{<:Real})
    Ns=length(s)
    os=OpSum()
    for j in 1:Ns-1
        if isodd(j)
            os += u, "N", j
            os += t1, "Cdag", j, "C", j+1
            os += t1', "Cdag", j+1, "C", j
            os += V[1], "N", j, "N", j+1
        else
            os += -u, "N", j
            os += t2, "Cdag", j, "C", j+1
            os += t2', "Cdag", j+1, "C", j
            os += V[2], "N", j, "N", j+1
        end
    end
    os += -u, "N", Ns
    return MPO(os,s)
end

# evolution gates of ITE
function SSH_gate(s::Vector{Index{Int}},t1::Number,t2::Number,dt::Float64)
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
function SSH_block(t1::Number,t2::Number,k::Float64;form=:array)
    hx=real(t1)+abs(t2)*cos(k+angle(t2))
    hy=-imag(t1)+abs(t2)*sin(k+angle(t2))
    if form ==:matrix
        return [0 hx-im*hy; hx+im*hy 0]
    else
        return [hx,hy]
    end
end

function SSH_block(t1::Number,t2::Number,u::Real,k::Float64;form=:array)
    hx=real(t1)+abs(t2)*cos(k+angle(t2))
    hy=-imag(t1)+abs(t2)*sin(k+angle(t2))
    if form ==:matrix
        return [u hx-im*hy ; hx+im*hy -u]
    else
        return [hx,hy,u]
    end
end

function SSH_band(t1::Number,t2::Number,k::Float64)
    at1,at2=abs(t1),abs(t2)
    phase=k+angle(t1)+angle(t2)
    return sqrt(at1*at1+at2*at2+2*at1*at2*cos(phase))
end
function SSH_band(t1::Number,t2::Number,u::Real,k::Float64)
    at1,at2=abs(t1),abs(t2)
    phase=k+angle(t1)+angle(t2)
    return sqrt(at1*at1+at2*at2+2*at1*at2*cos(phase)+u*u)
end

function SSH_spec(Lsize::Int,t1::Number,t2::Number)
    ks=range(-π,π,Int(Lsize/2))
    spectrum=SSH_band.(t1,t2,ks)
    sort(vcat(spectrum,-spectrum))
end
function SSH_spec(Lsize::Int,t1::Number,t2::Number,u::Real)
    ks=range(-π,π,Int(Lsize/2))
    spectrum=SSH_band.(t1,t2,u,ks)
    sort(vcat(spectrum,-spectrum))
end

function SSH_spectrum_obc(Lsize::Int,t1::Number,t2::Number;retstate::Bool=false)
    Ncell=Int(Lsize/2)
    arr=repeat([t1,t2],Ncell)[1:Lsize-1]
    Hamil=Hermitian(SymTridiagonal(zeros(Lsize),arr))
    if retstate==false
        return eigvals(Hamil)
    else
        return eigen(Hamil)
    end
end
function SSH_spectrum_obc(Lsize::Int,t1::Number,t2::Number,u::Real;retstate::Bool=false)
    Ncell = Int(Lsize/2)
    darr=repeat([u,-u],Ncell)
    uarr=repeat([t1,t2],Ncell)[1:Lsize-1]
    Hamil=Hermitian(SymTridiagonal(darr,uarr))
    if retstate==false
        return eigvals(Hamil)
    else
        return eigen(Hamil)
    end
end

function SSH_spectrum_pbc(Lsize::Int,t1::Number,t2::Number;retstate::Bool=false)
    Ncell=Int(Lsize/2)
    arr=repeat([t1,t2],Ncell)[1:Lsize-1]
    Hamil=diagm(1=>arr)
    Hamil[1,Lsize]=t2'
    Hamil=Hermitian(Hamil)
    if retstate==false
        return eigvals(Hamil)
    else
        return eigen(Hamil)
    end
end
function SSH_spectrum_pbc(Lsize::Int,t1::Number,t2::Number,u::Real;retstate::Bool=false)
    Ncell=Int(Lsize/2)
    darr=repeat([u,-u],Ncell)
    uarr=repeat([t1,t2],Ncell)[1:Lsize-1]
    Hamil=diagm(0=>darr,1=>uarr)
    Hamil[1,Lsize]=t2'
    Hamil=Hermitian(Hamil)
    if retstate==false
        return eigvals(Hamil)
    else
        return eigen(Hamil)
    end
end

function SSH_ED(Lsize::Int,t1::Number,t2::Number,mu::Real=0.0)
    spec=SSH_spectrum_obc(Lsize,t1,t2)
    energy=sum(spec[spec[:].<=mu])
    return energy
end
function SSH_ED(Lsize::Int,t1::Number,t2::Number,u::Real,mu::Real=0.0)
    spec=SSH_spectrum_obc(Lsize,t1,t2,u)
    energy=sum(spec[spec[:].<=mu])
    return energy
end


