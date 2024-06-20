using Statistics: mean

# computes degree at order d for each node - Eq. (4) in Lucas et al. (2020)
function degrees(M::Vector{Vector{Int64}})
    # input: M is a vector of index vectors representing *distinct* HOIs

    degs = sort(unique(length.(M))).-1 # degrees of HOIs present
    N = length(unique(reduce(vcat, M))) # number of nodes
    
    # build dict indexed by HOI order
    K = Dict([(d, zeros(Int, N)) for d in degs])
    for m in M # loop through HOIs
        d = length(m)-1
        for i in m
            K[d][i] += 1 # count number of order d HOIs involving node i
        end
    end

    return K
end

# calculates du/dt for higher-order Kuramoto model - Eq. (1) in Lucas et al. (2020)
function hoi_kuramoto_f!(du, u, p, t)
    ω = p.ω # natural freqs
    γ = p.γ # dict of coupling strengths at each order
    M = p.M # vector of vectors representing *distinct* HOIs

    # get average degree at each order across all nodes
    K = degrees(M)
    avg_K = Dict([(i, k) for (i, k) in zip(keys(K), mean.(collect(values(K))))])

    N = length(u)
    @inbounds for i in 1:N
        du[i] = ω[i] # natural frequency
        
        # loop through HOIs involving node i
        hois = filter(m -> in(i, m), M)
        for hoi in hois
            d = length(hoi) - 1
            du[i] += (γ[d] / avg_K[d]) * sin(sum(u[hoi]) - (d+1)*u[i])
        end
    end
    return
end