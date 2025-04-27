using ITensors, ITensorMPS
using Random

function ITensors.space(
    ::SiteType"MF"
)
  return 2
end

ITensors.state(::StateName"Emp",::SiteType"MF")=[1.0 0.0]
ITensors.state(::StateName"Occ",::SiteType"MF")=[0.0 1.0]
ITensors.state(::StateName"0",st::SiteType"MF")=state(StateName("Emp"),st)
ITensors.state(::StateName"1",st::SiteType"MF")=state(StateName("Occ"),st)

function ITensors.op!(Op::ITensor,::OpName"Gamma1",::SiteType"MF",s::Index)
    Op[s'=>1, s=>2] = 1.0
    return Op[s'=>2, s=>1] = 1.0
end
function ITensors.op!(Op::ITensor,::OpName"Gamma2",::SiteType"MF",s::Index)
    Op[s'=>1, s=>2] = -1.0im
    return Op[s'=>2,s=>1] = 1.0im
end
function ITensors.op!(Op::ITensor,::OpName"Gamma",::SiteType"MF",s::Index)
    Op[s'=>1,s=>1]=1.0
    return Op[s'=>2,s=>2]=-1.0
end
function ITensors.op!(Op::ITensor, ::OpName"N", ::SiteType"MF", s::Index)
    return Op[s'=>2,s=>2] = 1.0
end
function ITensors.op!(Op::ITensor,::OpName"F",::SiteType"MF",s::Index)
    Op[s'=>1,s=>1] = 1.0
    return Op[s'=>2,s=>2] = -1.0
end

function ITensors.op(::OpName"Gamma1", ::SiteType"Fermion", s::Index)
  c = op("C", s)
  cdag = op("Cdag", s)
  return c + cdag
end
function ITensors.op(::OpName"Gamma2", ::SiteType"Fermion", s::Index)
  c = op("C", s)
  cdag = op("Cdag", s)
  return -im * (c - cdag)
end

function ITensors.op(::OpName"N", ::SiteType"S=1/2", s::Index)
  return 0.5 * op("Id", s) - op("Sz", s)
end

ITensors.has_fermion_string(::OpName"Gamma1",::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma2", ::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma",::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma1",::SiteType"Fermion")=true
ITensors.has_fermion_string(::OpName"Gamma2", ::SiteType"Fermion")=true

function ITensorMPS.random_mps(siteo::Vector{<:Index}, sitee::Vector{<:Index}; linkdims=1)
  psio = random_mps(siteo; linkdims)
  psie = random_mps(sitee; linkdims)
  return [psio, psie]
end

function ITensorMPS.random_mps(rng::Random.AbstractRNG, siteo::Vector{<:Index}, sitee::Vector{<:Index}; linkdims=1)
  psio = random_mps(rng, siteo; linkdims)
  psie = random_mps(rng, sitee; linkdims)
  return [psio, psie]
end

function ITensorMPS.dmrg(Hs::Vector{MPO}, psi0s::Vector{MPS}, sweeps::Sweeps; kwargs...)
  ns = length(psi0s)
  Es = 0.0
  psis = MPS[]
  for i = 1:ns
    E, psi = dmrg(Hs[i], psi0s[i], sweeps; kwargs...)
    Es += E
    push!(psis, psi)
  end
  return Es, psis
end