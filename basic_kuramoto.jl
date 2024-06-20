# instantiate and precompile environment
using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

# load dependencies
using DifferentialEquations
using Distributions: Uniform
using Random; Random.seed!(11)
using Statistics: mean

# include plotting functions
include(joinpath(@__DIR__, "plotting.jl"))

# define ODE function to update in-place
function kuramoto_f!(du, u, p, t)
    ω = p.ω # natural freqs
    K = p.K # coupling strengths

    N = length(u)
    z = mean(exp.(im .* u)) # order parameter (centroid of all oscillators)
    ψ = angle(z)

    @inbounds for i in 1:N
        du[i] = ω[i] + K*abs(z)*sin(ψ - u[i])
    end

    return
end

# number of oscillators
N = 100

# problem parameters (natural freqs, coupling strengths)
p = (ω = rand(Uniform(0.8, 1.1), N), 
     K = 100/99
)

# initial condition
u0 = asin.(rand(Uniform(-0.5, 1.0), N))

# timespan over which to solve
tspan = (0.0, 10.0)

# define & solve problem
prob = ODEProblem(kuramoto_f!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# plot sin(u) for all oscillators
plot_activity(sol, :steelblue2; dt=0.1)

# animate phase for all oscillators
colors = color_by_freq(p.ω)
animate_phase(sol, colors, 50, "basic_pairwise_N100.mp4")