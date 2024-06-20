using CairoMakie
import ColorSchemes.viridis
using GLMakie
using Distributions: Uniform

function color_by_freq(ω)
    ω01 = (ω .- minimum(ω))./(maximum(ω) - minimum(ω)) # scale freqs between 0 and 1 for coloring
    return get(viridis, ω01)
end

function plot_activity(sol, color; dt::Union{Nothing, Float64}=nothing)
    # set timestep dt to plot finer resolution
    
    N = length(sol[:, 1])

    fig = Figure(fontsize=18)
    ax = Axis(fig[1,1])

    if dt == nothing
        for i in 1:N
            lines!(ax, sol.t, sin.(sol[i, :]), color=color[i], alpha=0.2)
        end
    else
        t = sol.t[1]:dt:sol.t[end]
        for i in 1:N
            sin_u = sin.(sol(collect(t), idxs=i).u)
            lines!(ax, t, sin_u, color=color[i], alpha=0.2)
        end
    end
    ax.xlabel=L"t"
    ax.xgridvisible=false
    ax.ylabel=L"\sin(\theta_i)"
    ax.ygridvisible=false

    return fig
end

function plot_phase_coherence(sol)
    r = [abs(mean(exp.(im .* u))) for u in sol.u]

    fig = Figure(fontsize=18)
    ax = Axis(fig[1,1], limits=(nothing, (0.0, 1.0)))
    scatter!(ax, sol.t, r)
    ax.xlabel=L"t"
    ax.xgridvisible=false
    ax.ylabel=L"\text{Global order parameter}"
    ax.ygridvisible=false

    return fig
end

function animate_phase(sol, color, framerate::Int, filename::String)
    time = Observable(sol.t[1])

    u = @lift(sol($time))
    r = rand(Uniform(0.8, 1.0), length(sol[1])) # offsets pts for visibility

    # figure to plot into
    fig = Figure()
    ax = PolarAxis(fig[1,1], 
                   spinevisible=false
    )
    rlims!(ax, 0.0, 1.1)
    hidedecorations!(ax)
    scatter!(ax, u, r, markersize=15, color=color)
    suptitle = Label(fig[0, 1], 
                     @lift("t = $(round($time, digits = 1))"), 
                     tellwidth=false
    )
    
    framerate = framerate
    timestamps = range(sol.t[1], sol.t[end], step=1/framerate)
    
    record(fig, filename, timestamps; framerate=framerate) do t
        time[] = t
    end
    return
end