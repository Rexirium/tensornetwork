using ITensors, ITensorMPS

ITensors,space(::SiteType"MF")=2
ITensors.op(::OpName"Gamma",::SiteType"MF")=[1 0 ; 0 -1]

let
    L,D=40,5
    v,w=2.0,1.0
    N=2*L 
    sw=Sweeps(15)
    setmaxdim!(sw,100)
    setcutoff!(sw,1E-14)
    krydim=4
end
