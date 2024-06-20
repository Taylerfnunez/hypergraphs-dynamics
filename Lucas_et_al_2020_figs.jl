# instantiate and precompile environment
using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

# load dependencies
using Combinatorics: combinations
using DifferentialEquations
using Distributions
using Random; Random.seed!(11)

# include functions from other scripts
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "hoi_kuramoto.jl"))

# star-clique model from Lucas et al. (2020)
function star_clique()
    # 6-node star and 7-node clique
    M = vcat(
        [[1,6]],
        [[2,6]],
        [[3,6]],
        [[4,6]],
        [[5,6]],
        [[6,7]],
        collect(combinations(7:13, 2)),
        collect(combinations(7:13, 3)),
        collect(combinations(7:13, 4)),
        collect(combinations(7:13, 5)),
        collect(combinations(7:13, 6)),
        [collect(7:13)]
    )
    γ = Dict([(d, 1.0) for d in 1:6]) # coupling strengths at each order
    ω = ones(13) # natural freqs

    p = (ω = ω, γ = γ, M = M)

    u0 = asin.(1.5*rand(13).-0.5)
    tspan = (0.0, 15.0)
    
    prob = ODEProblem(hoi_kuramoto_f!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    # color plot by clique and star membership
    color=[fill(:red, 6); fill(:blue, 7)]
    return plot_activity(sol, color)
end

function fig_2a()
    # could also do this with basic kuramoto model
    N = 100
    M = collect(combinations(1:N, 2))
    γ = Dict([(1, 1.0)])
    ω = ones(N)

    p = (ω = ω, γ = γ, M = M)

    u0 = asin.(1.5*rand(N).-0.5)
    tspan = (0.0, 10.0)
    
    prob = ODEProblem(hoi_kuramoto_f!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    color = fill(:steelblue2, N)
    return plot_activity(sol, color; dt=0.1)
end

function fig_2b()
    N = 100
    M = collect(combinations(1:N, 3))
    γ = Dict([(2, 1.0)])
    ω = ones(N)

    p = (ω = ω, γ = γ, M = M)

    u0 = asin.(1.5*rand(N).-0.5)
    tspan = (0.0, 10.0)

    prob = ODEProblem(hoi_kuramoto_f!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    color = fill(:steelblue2, N)
    return plot_activity(sol, color; dt=0.1)
end

function fig_3a()
    N = 100
    M = vcat(collect(combinations(1:N, 2)), collect(combinations(1:N, 3)))
    γ = Dict([(1, -0.5), (2, 0.5)])
    ω = ones(N)

    p = (ω = ω, γ = γ, M = M)

    u0 = asin.(rand(Uniform(-0.8, 0.8), N))
    tspan = (0.0, 10.0)

    prob = ODEProblem(hoi_kuramoto_f!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    color = fill(:steelblue2, N)
    return plot_activity(sol, color; dt=0.1)
end

function fig_3b()
    N = 100
    M = vcat(collect(combinations(1:N, 2)), collect(combinations(1:N, 3)))
    γ = Dict([(1, 0.5), (2, -0.5)])
    ω = ones(N)

    p = (ω = ω, γ = γ, M = M)

    u0 = asin.(rand(Uniform(-0.8, 0.8), N))
    tspan = (0.0, 10.0)

    prob = ODEProblem(hoi_kuramoto_f!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    color = fill(:steelblue2, N)
    return plot_activity(sol, color; dt=0.1)
end