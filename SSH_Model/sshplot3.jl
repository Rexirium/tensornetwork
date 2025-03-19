using HDF5, Interpolations
using Plots, LinearAlgebra
using LaTeXStrings, Printf

num = 10
Num = 50
Vs = LinRange(0.0,2.0, num+1)
Ws = LinRange(0.0,2.0, num+1)
Vi = LinRange(0.0, 2.0, Num+1) # interpolated V axis
Wi = LinRange(0.0,2.0, Num+1)  # interpolated W axis
# color gradient
cg1=cgrad(:thermal, rev = true, categorical=false)
cg2=cgrad(:viridis, rev = false, categorical=false)

# interaction in topological trivial phase
P1 = let 
    file = h5open("sshintdata.h5", "r")
    energies = read(file, "energiesv2.0w1.0")
    entropies = read(file, "entropiesv2.0w1.0")
    close(file)
    # cubic interpolation used
    energies_itp = cubic_spline_interpolation((Vs,Ws), transpose(energies))
    entropies_itp = cubic_spline_interpolation((Vs,Ws), transpose(entropies))
    energi = energies_itp(Vi,Wi)
    entropi = entropies_itp(Vi, Wi)

    p1 = heatmap(Vi, Wi, energi,  clim=(-42.5,-32.5), colorbar=false,
        c=cg1 ,framestyle=:box, title="energy "*(L"v=2,w=1"))
    p2 = heatmap(Vi, Wi, entropi, xlabel=L"V/v", clim=(0.0,1.2), colorbar=false,
    c=cg2, framestyle=:box, title="entropy "*(L"v=2,w=1"))
    plot(p1,p2, layout = @layout([a; b]),size = (400,800))
end

# interaction in topological non-trivial phase
P2 = let 
    file = h5open("sshintdata.h5", "r")
    energies = read(file, "energiesv1.0w2.0")
    entropies = read(file, "entropiesv1.0w2.0")
    close(file)
    # cubic interpolation used
    energies_itp = cubic_spline_interpolation((Vs,Ws), transpose(energies))
    entropies_itp = cubic_spline_interpolation((Vs,Ws), transpose(entropies))
    energi = energies_itp(Vi,Wi)
    entropi = entropies_itp(Vi, Wi)

    p1 = heatmap(Vi, Wi, energi, ylabel=L"W/w", clim=(-42.5,-32.5),
        c=cg1 ,framestyle=:box, title="energy "*(L"v=1,w=2"))
    p2 = heatmap(Vi, Wi, entropi, xlabel=L"V/v", ylabel=L"W/w", clim=(0,1.2),
        c=cg2, framestyle=:box, title="entropy "*(L"v=1,w=2"))
    plot(p1,p2, layout = @layout([a; b]),size = (600,800))
end

# interaction in critical point
P3 = let 
    file = h5open("sshintdata.h5", "r")
    energies = read(file, "energiesv1.0w1.0")
    entropies = read(file, "entropiesv1.0w1.0")
    close(file)
    # cubic interpolation used
    energies_itp = cubic_spline_interpolation((Vs,Ws), transpose(energies))
    entropies_itp = cubic_spline_interpolation((Vs,Ws), transpose(entropies))
    energi = energies_itp(Vi,Wi)
    entropi = entropies_itp(Vi, Wi)

    p1 = heatmap(Vi, Wi, energi, xlabel=L"V/v",
        c=cg1 ,framestyle=:box, title="energy "*(L"v=1,w=1"))
    p2 = heatmap(Vi, Wi, entropi, xlabel=L"V/v", ylabel=L"W/w", 
        c=cg2, framestyle=:box, title="entropy "*(L"v=1,w=1"))
    plot(p1,p2, layout = @layout([a b]),size = (800,400))
end

lay = @layout([a{0.45w} b{0.55w}])
plot(P1,P2, layout=lay,size = (1000,800))
