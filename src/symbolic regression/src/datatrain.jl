# Step 1. Develop Training Data

using CSV
using DataFrames
using Distributions
using LRP
using Plots
using Random
using Serialization
using StatsBase

# Fetch input-output (data) for instances with variation in customer location and demand
# NOTE: check-point seeding ensures consistent customer location and demand values across different instances
function io(customer_distribution, demand_distribution, seed)
    dir = joinpath(dirname(pwd()), "data/instances/train")
    
    # Step 1.1 Define all the input parameters for the distribution environment
    rng= MersenneTwister(seed+0)
    lˣ = sample(rng, 1:9)
    lʸ = sample(rng, 1:9)
    ρˣ = sample(rng, -1:0.1:1) * lˣ
    ρʸ = sample(rng, -1:0.1:1) * lʸ
    qᵛ = sample(rng, 1:10) * 10
    n  = sample(rng, 100:500)
    rng= MersenneTwister(seed+1)
    Xᶜ = zeros(n)
    Yᶜ = zeros(n)
    if isequal(customer_distribution, "unfrm")
        Xᶜ = rand(rng, Uniform(-0.5,0.5), n) .* lˣ
        Yᶜ = rand(rng, Uniform(-0.5,0.5), n) .* lʸ
    elseif isequal(customer_distribution, "tnglr")
        Xᶜ = rand(rng, TriangularDist(-0.5, 0.5, rand(rng, Uniform(-0.5,0.5))), n) .* lˣ
        Yᶜ = rand(rng, TriangularDist(-0.5, 0.5, rand(rng, Uniform(-0.5,0.5))), n) .* lʸ
    elseif isequal(customer_distribution, "norml")
        Xᶜ = rand(rng, Truncated(Normal(rand(rng, Uniform(-0.5,0.5)), rand(rng, Uniform(0.1,0.4))), -0.5, 0.5), n) .* lˣ
        Yᶜ = rand(rng, Truncated(Normal(rand(rng, Uniform(-0.5,0.5)), rand(rng, Uniform(0.1,0.4))), -0.5, 0.5), n) .* lʸ
    elseif isequal(customer_distribution, "expnl")
        Xᶜ = sample(rng, [-1,1]) .* (rand(rng, Truncated(Exponential(rand(rng, Uniform(0.,0.5))), 0., 0.5), n) .- 0.25) .* 2lˣ
        Yᶜ = sample(rng, [-1,1]) .* (rand(rng, Truncated(Exponential(rand(rng, Uniform(0.,0.5))), 0., 0.5), n) .- 0.25) .* 2lʸ
    elseif isequal(customer_distribution, "const")
        nˣ = floor(Int, sqrt(n))
        while n % nˣ ≠ 0 nˣ -= 1 end
        nʸ = n ÷ nˣ
        n  = isone(nˣ) ? n - 1 : n
        Xᶜ = [(-0.5 + 1/(nˣ - 1) * (i - 1)) for i ∈ 1:nˣ for _ ∈ 1:nʸ] .* lˣ
        Yᶜ = [(-0.5 + 1/(nʸ - 1) * (i - 1)) for _ ∈ 1:nˣ for i ∈ 1:nʸ] .* lʸ
    elseif isequal(customer_distribution, "bndry")
        Xᶜ = zeros(n)
        Yᶜ = zeros(n)
        k  = 0
        while k != n
            x = rand(rng, Uniform(-0.5,0.5)) .* lˣ
            y = rand(rng, Uniform(-0.5,0.5)) .* lʸ
            z = (abs(x) / (lˣ/2)) * (abs(y) / (lʸ/2))
            if z > rand(rng)
                k += 1
                Xᶜ[k] = x
                Yᶜ[k] = y
            end
        end
    end
    rng= MersenneTwister(seed+2)
    Qᶜ = zeros(n)
    if isequal(demand_distribution, "unfrm")
        Qᶜ = rand(rng, Uniform(0, 2), n)
    elseif isequal(demand_distribution, "tnglr")
        Qᶜ = rand(rng, TriangularDist(0, 2, rand(rng, Uniform(0,2))), n)
    elseif isequal(demand_distribution, "norml")
        Qᶜ = rand(rng, Truncated(Normal(rand(rng, Uniform(0,2)), rand(rng, Uniform(0.2,0.8))), 0, 2), n)
    elseif isequal(demand_distribution, "expnl")
        Qᶜ = sample(rng, [rand(rng, Truncated(Exponential(rand(rng, Uniform(0,2))), 0, 2), n), 2 .- rand(rng, Truncated(Exponential(rand(rng, Uniform(0,2))), 0, 2), n)])
    elseif isequal(demand_distribution, "const")
        Qᶜ = ones(n)
    elseif isequal(demand_distribution, "bndry")
        k  = 0
        while k != n
            q = rand(rng, Uniform(0,2))
            z = abs(q-1)
            if z > rand(rng)
                k += 1
                Qᶜ[k] = q
            end
        end
    end
    
    # Step 1.2 Visualize the distribution environment
    fig = plot(legend=:none)
    scatter!(fig, [ρˣ], [ρʸ], markershape=:rect, markersize=6, markerstrokewidth=0.8, color=palette(:default)[2])
    for i ∈ 1:n scatter!(fig, [Xᶜ[i]], [Yᶜ[i]], markershape=:circle, markersize=4, markerstrokewidth=0.5, color=palette(:default)[1], alpha=Qᶜ[i]/maximum(Qᶜ)) end
    vline!([-lˣ/2,lˣ/2], color="grey")
    hline!([-lʸ/2,lʸ/2], color="grey")
    display(fig)
    
    # Step 1.3 Develop input CSV files for LRP.jl
    dfᵈ = DataFrame(in=[1], x=[ρˣ], y=[ρʸ], qd=[sum(Qᶜ)], ts=[0], te=[0], co=[0], cf=[0], phi=[0])
    dfᶜ = DataFrame(in=2:(n+1), x=Xᶜ, y=Yᶜ, qc=Qᶜ, tc=zeros(n), tr=zeros(n), te=zeros(n), tl=zeros(n))
    dfᵛ = DataFrame(iv=[1], jv=[1], id=[1], qv=[qᵛ], lv=[2n*(lˣ+lʸ)], sv=[1], tf=[0], td=[0], tc=[0], tw=[0], r=[n], cd=[1], ct=[0], cf=[0])
    A   = [sqrt(((isone(j) ? ρˣ : Xᶜ[j-1]) - (isone(i) ? ρˣ : Xᶜ[i-1]))^2 + ((isone(j) ? ρʸ : Yᶜ[j-1]) - (isone(i) ? ρʸ : Yᶜ[i-1]))^2) + eps(1e7) * rand(rng) for i ∈ 1:(n+1), j ∈ 1:(n+1)]
    dfᵃ = DataFrame(A[2:end,:], string.(A[1,:]))
    mkdir("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed")
    
    CSV.write("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/depot_nodes.csv", dfᵈ)
    CSV.write("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/customer_nodes.csv", dfᶜ)
    CSV.write("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/vehicles.csv", dfᵛ)
    CSV.write("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/arcs.csv", dfᵃ)

    # Step 1.4 Solve the instance
    rng= MersenneTwister(seed+3)
    s₁ = initialize(rng, "$dir/loc-dem/$customer_distribution/$demand_distribution/$seed"; method=:global, dir=joinpath(dirname(@__DIR__), "instances"))
    display(visualize(s₁))
    x = max(100, lastindex(s₁.C))
    χ = ALNSparameters(
        j   =   50                      ,
        k   =   5                       ,
        n   =   x                       ,
        m   =   100x                    ,
        Ψᵣ  =   [
                    :randomcustomer!    ,
                    :randomroute!       ,
                    :randomvehicle!     ,
                    :relatedcustomer!   ,
                    :relatedroute!      ,
                    :relatedvehicle!    ,
                    :worstcustomer!     ,
                    :worstroute!        ,
                    :worstvehicle!
                ]                       ,
        Ψᵢ  =   [
                    :best!              ,
                    :precise!           ,
                    :perturb!           ,
                    :regret2!           ,
                    :regret3!
                ]                       ,
        Ψₗ  =   [
                    :intramove!         ,
                    :intraswap!         ,
                    :intraopt!          ,
                    :intermove!         ,
                    :interswap!         ,
                    :interopt!
                ]                       ,
        σ₁  =   15                      ,
        σ₂  =   10                      ,
        σ₃  =   3                       ,
        μ̲   =   0.1                     ,
        c̲   =   4                       ,
        μ̅   =   0.4                     ,
        c̅   =   60                      ,
        ω̅   =   0.05                    ,
        τ̅   =   0.5                     ,
        ω̲   =   0.01                    ,
        τ̲   =   0.01                    ,
        θ   =   0.9985                  ,
        ρ   =   0.1
    )
    s₂ = ALNS(rng, χ, s₁)
    display(visualize(s₂))
    serialize("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/solution.dat", s₂)
    
    # Step 1.5 Extract characteristics of the distribution environment
    s  = deserialize("$dir/loc-dem/$customer_distribution/$demand_distribution/$seed/solution.dat")
    K  = [findmin([isequal(i,j) ? Inf : sqrt((Xᶜ[j]-Xᶜ[i])^2 + (Yᶜ[j]-Yᶜ[i])^2) for j ∈ 1:n])[2] for i ∈ 1:n]
    L  = [sqrt((Xᶜ[K[i]]-Xᶜ[i])^2 + (Yᶜ[K[i]]-Yᶜ[i])^2) for i ∈ 1:n]
    r  = (maximum(Xᶜ) - minimum(Xᶜ)) * (maximum(Yᶜ) - minimum(Yᶜ)) / (lˣ * lʸ)
    ρ̄ᵘ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n])
    ρ̄ʷ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n], weights(Qᶜ))
    m  = ceil(Int, sum(Qᶜ) / qᵛ)
    μˢ = mean(L)
    σˢ = std(L)
    γˢ = (sum([(L[i] - μˢ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
    κˢ = (sum([(L[i] - μˢ)^4 for i ∈ 1:n]) / n) # ^ (1/4)
    μᵈ = mean(Qᶜ)
    σᵈ = std(Qᶜ)
    γᵈ = (sum([(Qᶜ[i] - μᵈ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
    κᵈ = (sum([(Qᶜ[i] - μᵈ)^4 for i ∈ 1:n]) / n) # ^ (1/4)
    z  = sum([r.l for d ∈ s.D for v ∈ d.V for r ∈ v.R])
    df = CSV.read("$dir/data.csv", DataFrame)
    push!(df, (nrow(df)+1, "$customer_distribution-$demand_distribution$seed", n, lˣ, lʸ, ρˣ, ρʸ, qᵛ, r, ρ̄ᵘ, ρ̄ʷ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ, z), promote=true)
    CSV.write("$dir/data.csv", df)

    return
end

seeds = [
    1036, 1444, 1896,
    2104, 2403, 2706, 2984,
    3020, 3605, 3811,
    4169, 4479, 4903,
    5020, 5193, 5438, 5722,
    6056, 6503, 6810,
    7033, 7531, 7950,
    8044, 8448, 8620, 8892,
    9092, 9378, 9666
]
for customer_distribution ∈ ["unfrm", "tnglr", "norml", "expnl", "const", "bndry"]
    for demand_distribution ∈ ["unfrm", "tnglr", "norml", "expnl", "const", "bndry"]
        for seed ∈ seeds io(customer_distribution, demand_distribution, seed) end
    end
end

