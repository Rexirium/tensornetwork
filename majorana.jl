using ITensors, ITensorMPS

function ITensors.space(
    ::SiteType"MF";
    conserve_qns=false,
    conserve_nf=conserve_qns,
    conserve_nfparity=conserve_qns,
    qnname_nf="Nf",
    qnname_nfparity="NfParity",
    qnname_sz="Sz",
    conserve_sz=false,
    conserve_parity=nothing,
)
  if !isnothing(conserve_parity)
    conserve_nfparity = conserve_parity
  end
  if conserve_sz == true
    conserve_sz = "Up"
  end
  if conserve_nf && conserve_sz == "Up"
    zer = QN((qnname_nf, 0, -1), (qnname_sz, 0)) => 1
    one = QN((qnname_nf, 1, -1), (qnname_sz, 1)) => 1
    return [zer, one]
  elseif conserve_nf && conserve_sz == "Dn"
    zer = QN((qnname_nf, 0, -1), (qnname_sz, 0)) => 1
    one = QN((qnname_nf, 1, -1), (qnname_sz, -1)) => 1
    return [zer, one]
  elseif conserve_nfparity && conserve_sz == "Up"
    zer = QN((qnname_nfparity, 0, -2), (qnname_sz, 0)) => 1
    one = QN((qnname_nfparity, 1, -2), (qnname_sz, 1)) => 1
    return [zer, one]
  elseif conserve_nfparity && conserve_sz == "Dn"
    zer = QN((qnname_nfparity, 0, -2), (qnname_sz, 0)) => 1
    one = QN((qnname_nfparity, 1, -2), (qnname_sz, -1)) => 1
    return [zer, one]
  elseif conserve_nf
    zer = QN(qnname_nf, 0, -1) => 1
    one = QN(qnname_nf, 1, -1) => 1
    return [zer, one]
  elseif conserve_nfparity
    zer = QN(qnname_nfparity, 0, -2) => 1
    one = QN(qnname_nfparity, 1, -2) => 1
    return [zer, one]
  end
  return 2
end

ITensors.state(::StateName"Emo",::SiteType"MF")=[1.0 0.0]
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
function ITensors.op!(Op::ITensor,::OpName"F",::SiteType"MF",s::Index)
    Op[s'=>1,s=>1] = 1.0
    return Op[s'=>2,s=>2] = -1.0
end

function ITensors.op!(Op::ITensor,::OpName"Gamma1",::SiteType"Fermion",s::Index)
    Op[s'=>1, s=>2] = 1.0
    return Op[s'=>2, s=>1] = 1.0
end
function ITensors.op!(Op::ITensor,::OpName"Gamma2",::SiteType"Fermion",s::Index)
    Op[s'=>1, s=>2] = -1.0im
    return Op[s'=>2,s=>1] = 1.0im
end

ITensors.has_fermion_string(::OpName"Gamma1",::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma2", ::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma",::SiteType"MF")=true
ITensors.has_fermion_string(::OpName"Gamma1",::SiteType"Fermion")=true
ITensors.has_fermion_string(::OpName"Gamma2", ::SiteType"Fermion")=true