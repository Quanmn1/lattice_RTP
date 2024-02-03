# for installing required packages
# using Pkg   #for installing packages

# Pkg.add("Plots")
# Pkg.add("Random")

using Plots     #for plotting heatmap
using Random    #for Mersenne twister
using Printf


cd(@__DIR__)
include("lattice_rtp_module.jl")   #functions
include("lattice_functions.jl")

#random number generator
rng = MersenneTwister(1234)

v, alpha = 10, 1       #diffusion coefficient,  velocity,  flipping rate
Lx, Ly = 20, 20       #x and y extension of the system
rho_0, rho_m = 10, 20         #temperature and the average density

pr = RTP.Param(v, alpha, Lx, Ly, rho_0, rho_m)  #define a set of parameters

pos0 = zeros(Int64, pr.N, 2)    #initial position
s0 = zeros(Int8, pr.N, 2)          #initial spins

#initial condition: random
for i in 1:pr.N
    pos0[i,1], pos0[i,2] = rand(rng, 1:Lx), rand(rng, 1:Ly)
    dir = rand(rng, 1:4)
    s0[i,1], s0[i,2] = RTP.dirs[dir][1], RTP.dirs[dir][2]
end

alpha_min = 0.4
alpha_max = 2.0
alpha_step = 0.8

lp_min = v/alpha_max
lp_max = v/alpha_min
# phase diagram:
# x axis: density
# y axis: persistance length 
plt = plot(
    ylim = (0, 1.2*lp_max),
    xlim = (0, 1.5*rho_m),
    title = "Phase diagram",
)

vline!([rho_m])

for alpha in range(alpha_min,step=alpha_step,stop=alpha_max)
    pr.alpha = alpha
    lp = pr.v/pr.alpha
    t = 0.0

    st = RTP.State(t, pos0, s0,  Lx, Ly)    #define state variables

    #state, param, t_gap, n_frame, rng, fps
    make_movie!(st, pr, 10, 100, rng, 20)

    density_counts = density_histogram(st.rho)
    rho_max = length(density_counts) - 1
    densities = 0:rho_max
    plt2 = bar(densities, density_counts)
    savefig(plt2, @sprintf("Density_histogram_%.2f.png", lp))

    divider = Int(round(0.6*rho_m))
    rho_g = argmax(density_counts[1:divider])-1
    rho_l = argmax(density_counts[divider+1:rho_max+1])-1 + divider

    scatter!(plt, [rho_g], [lp], color=:blue)
    scatter!(plt, [rho_l], [lp], color=:red)
end

savefig(plt, "Phase_diagram.png")
