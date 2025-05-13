using Plots
using LaTeXStrings,Printf
include("ssh_hamil.jl")


L, num = 20, 100
cs = palette(:default)[1:3]


p1=let 
    v=1.0
    ws=LinRange(0, 2.5, num+1)
    spec = SSH_spectrum_obc.(L, v, ws)
    data = transpose(hcat(spec...))
    plot(ws,data,w=1.5,leg=false, xlabel=L"w", ylabel=L"E", framestyle=:box, c=cs[1])
    vline!([1.0], line=(1,:dash))
    annotate!([(0.5, 3, ("a) "*(L"v=1")))])
end
p2=let 
    w = 1.0
    vs = LinRange(0, 2.5, num+1)
    spec = SSH_spectrum_obc.(L,vs,w)
    data = transpose(hcat(spec...))
    plot(vs, data, w=1.5, leg=false, xlabel=L"v", ylabel=L"E", framestyle=:box, c=cs[2])
    vline!([1.0],line=(1,:dash))
    annotate!([(0.5, 3, ("b) "*(L"w=1")))])
end
p3=let 
    w = 1.0
    th = LinRange(0.0, 1.0, num+1)
    vx = 1.5 .- cos.(π*th)
    ux = sin.(π*th)
    spec = SSH_spectrum_obc.(L,vx,w,ux)
    data = transpose(hcat(spec...))
    plot(th, data, w=1.5, leg=false, xlabel=L"\theta", ylabel=L"E", framestyle=:box, c=:black)
    vline!([1/3], line=(1,:dash))
    annotate!([(0.2, 3, ("c) "*(L"w=1")))])
end
p4=let
    v, w = 0.5, 1.0
    ls = [11, 10, 2]
    lay = @layout([a;b;c])
    colors = repeat([:red, :blue], Int(L/2))
    spec, states = SSH_spectrum_obc(L,v,w;retstate=true)
    b1 = bar(states[:, ls[1]], xtick = false, framestyle = :box, ylim =(-0.7, 0.7), 
        ytick =[-0.6, 0.0, 0.6], legend = false, color = colors)
    annotate!([(30, 0.5, "ϵ = +0.0")])
    b2 = bar(states[:, ls[2]], xtick = false, framestyle = :box, ylim = (-0.4, 0.7), 
        ytick =[-0.4, 0.0, 0.6], legend =false, color =colors)
    annotate!([(30, 0.5, "ϵ = -0.0")])
    b3 = bar(states[:, ls[3]], xlabel = "site "*(L"j"),framestyle = :box, ylim =(-0.3, 0.4), 
        ytick =[-0.3, 0.0, 0.3], legend = false, color = colors)
    annotate!([(30, 0.3, "ϵ = -1.4")])
    plot(b1, b2, b3, layout = lay)
end

layout=@layout([a{0.45h} b{0.45h};c d])
P = plot(p1, p2, p3, p4, layout=layout, size=(800, 600))
savefig(P, "SSH_Model/sshfigs/ssh_spectrum.svg")

