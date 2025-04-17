"""
    remove!([rng::AbstractRNG], q::Int, s::Solution, method::Symbol)

Returns solution removing `q` customer nodes from solution s using the given `method`.

Available methods include,
- Random Customer Node Removal  : `:randomcustomer!`
- Random Route Removal          : `:randomroute!`
- Random Vehicle Removal        : `:randomvehicle!`
- Random Depot Node Removal     : `:randomdepot!` 
- Related Customer Node Removal : `:relatedcustomer!`
- Related Route removal         : `:relatedroute!`
- Related Vehicle Removal       : `:relatedvehicle!`
- Related Depot Node Removal    : `:relateddepot!`
- Worst Customer Node Removal   : `:worstcustomer!`
- Worst Route Removal           : `:worstroute!`
- Worst Vehicle Removal         : `:worstvehicle!`
- Worst Depot NodeRemoval       : `:worstdepot!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, q::Int, s::Solution, method::Symbol)::Solution = isdefined(MetaHeuristic, method) ? getfield(MetaHeuristic, method)(rng, q, s) : getfield(Main, method)(rng, q, s)
remove!(q::Int, s::Solution, method::Symbol) = remove!(Random.GLOBAL_RNG, q, s, method)



# -------------------------------------------------- NODE REMOVAL --------------------------------------------------
"""
    randomcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes
selected randomly.
"""
function randomcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    W = isactive.(C)                # W[iⁿ]: selection weight of customer node C[iⁿ]
    # Step 2: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    while n < q
        iⁿ = sample(rng, eachindex(C), OffsetWeights(W))
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        W[iⁿ] = 0
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



"""
    relatedcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes
most related to a randomly selected pivot customer node.
"""
function relatedcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    X = fill(-Inf, eachindex(C))    # X[iⁿ]: relatedness of customer node C[iⁿ] with pivot customer node C[i]
    W = isactive.(C)                # W[iⁿ]: selection weight of customer node C[iⁿ]
    # Step 2: Randomly select a pivot customer node
    i = sample(rng, eachindex(C), OffsetWeights(W))
    # Step 3: For each customer node, evaluate relatedness to this pivot customer node
    m = sample(rng, s.φ ? [:q, :l, :t] : [:q, :l])
    for iⁿ ∈ eachindex(C) X[iⁿ] = isone(W[iⁿ]) ? relatedness(m, C[iⁿ], C[i], s) : -Inf end
    X[i] = Inf
    # Step 4: Remove q most related customer nodes
    n = 0
    while n < q
        iⁿ = argmax(X)
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        X[iⁿ] = -Inf
    end
    # Step 5: Remove redundant vehicles and routes
    postremove!(s)
    # Step 6: Return solution
    return s
end



"""
    worstcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes 
with highest removal cost (savings).
"""
function worstcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R if isactive(r)]
    L = [c for c ∈ C if isactive(c)]
    X = fill(-Inf, eachindex(L))   # X[i] : removal cost of customer node L[i]
    ϕ = ones(Int, eachindex(R))    # ϕʳ[j]: binary weight for route R[j]
    # Step 2: Iterate until q customer nodes have been removed
    n = 0
    while n < q
        # Step 2.1: For every closed customer node evaluate removal cost
        z = f(s)
        for (i,c) ∈ pairs(L)
            r = c.r
            j = findfirst(isequal(r), R)
            if isopen(c) || iszero(ϕ[j]) continue end
            # Step 2.1.1: Remove closed customer node c between tail node nᵗ and head node nʰ in route r
            nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            # Step 2.1.2: Evaluate the removal cost
            z′ = f(s) * (1 + rand(rng, Uniform(-0.2, 0.2)))
            Δ  = z′ - z
            X[i] = -Δ
            # Step 2.1.3: Re-insert customer node c between tail node nᵗ and head node nʰ in route r
            insertnode!(c, nᵗ, nʰ, r, s)
        end
        # Step 2.2: Remove the customer node with highest removal cost (savings)
        i  = argmax(X)
        c  = L[i]
        r  = c.r
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        # Step 2.3: Update cost and selection weight vectors
        ϕ .= 0
        X[i] = -Inf
        for (j,r) ∈ pairs(R) 
            φʳ = isequal(r, c.r)
            φᵛ = isequal(r.iᵛ, v.iᵛ) && isless(c.r.tⁱ, r.tⁱ) && isequal(s.φ, true)
            φᵈ = false
            φˢ = φʳ || φᵛ || φᵈ
            if isequal(φˢ, false) continue end
            ϕ[j] = 1
        end
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



# -------------------------------------------------- ROUTE REMOVAL --------------------------------------------------
"""
    randomroute!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after iteratively selecting a random route and 
removing customer nodes from it until exactly `q` customer nodes are
removed.
"""
function randomroute!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    W = isopt.(R) .* isactive.(R)   # W[iʳ]: selection weight for route R[iʳ]
    # Step 2: Iteratively select a random route and remove customer nodes from it until exactly q customer nodes are removed
    n = 0
    while n < q
        iʳ = sample(rng, eachindex(R), Weights(W))
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            if !isopt(r) || isdormant(r) break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end
        W[iʳ] = 0
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



"""
    relatedroute!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer 
nodes from the routes most related to a randomly selected 
pivot route.
"""
function relatedroute!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    X = fill(-Inf, eachindex(R))    # X[iʳ]: relatedness of route R[iʳ] with pivot route R[i]
    W = isopt.(R) .* isactive.(R)   # W[iʳ]: selection weight for route R[iʳ]
    # Step 2: Randomly select a pivot route
    i = sample(rng, eachindex(R), Weights(W))  
    # Step 3: For each route, evaluate relatedness to this pivot route
    m = sample(rng, s.φ ? [:l, :t] : [:l])
    for iʳ ∈ eachindex(R) X[iʳ] = isone(W[iʳ]) ? relatedness(m, R[iʳ], R[i], s) : -Inf end
    X[i] = Inf
    # Step 4: Remove exactly q customers from most related route to this pivot route
    n = 0
    while n < q
        iʳ = argmax(X)
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            if !isopt(r) || isdormant(r) break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end 
        X[iʳ] = -Inf
        W[iʳ] = 0
    end
    postremove!(s)
    # Step 5: Return solution
    return s
end



"""
    worstroute!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer 
nodes from low-utilization routes.
"""
function worstroute!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    X = fill(Inf, eachindex(R))     # X[iʳ]: utilization of route R[iʳ]
    W = isopt.(R) .* isactive.(R)   # W[iʳ]: selection weight for route R[iʳ]
    # Step 2: Evaluate utilization of each route
    for (iʳ,r) ∈ pairs(R)
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        X[iʳ] = isone(W[iʳ]) ? r.q/v.qᵛ : Inf
    end
    # Step 3: Iteratively select low-utilization route and remove customer nodes from it until exactly q customer nodes are removed
    n = 0
    while n < q
        iʳ = argmin(X)
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            if !isopt(r) || isdormant(r) break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end
        X[iʳ] = Inf
        W[iʳ] = 0
    end
    postremove!(s)
    # Step 4: Return solution
    return s
end    



# -------------------------------------------------- VEHICLE REMOVAL --------------------------------------------------
"""
    randomvehicle!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after iteratively selecting a random vehicle and 
removing customer nodes from its routes until at least `q` customer nodes 
are removed.
"""
function randomvehicle!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    V = [v for d ∈ D for v ∈ d.V]
    W = isopt.(V)                   # W[iᵛ]: selection weight for vehicle V[iᵛ]
    # Step 2: Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        iᵛ = sample(rng, eachindex(V), Weights(W))
        v  = V[iᵛ]
        d  = D[v.iᵈ]
        for r ∈ v.R
            if n ≥ q break end
            while true
                if !isopt(r) || isdormant(r) break end
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        W[iᵛ] = 0
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



"""
    relatedvehicle!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing at least `q` customer nodes
from the routes of the vehicles most related to a randomly 
selected pivot vehicle.
"""
function relatedvehicle!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    V = [v for d ∈ D for v ∈ d.V]
    X = fill(-Inf, eachindex(V))    # X[iᵛ]: relatedness of vehicle V[iᵛ] with pivot vehicle V[i]
    W = isopt.(V)                   # W[iᵛ]: selection weight for vehicle V[iᵛ]
    # Step 2: Select a random vehicle
    i = sample(rng, eachindex(V), Weights(W))
    # Step 3: For each vehicle, evaluate relatedness to this pivot vehicle
    m = sample(rng, s.φ ? [:l, :t] : [:l])
    for iᵛ ∈ eachindex(V) X[iᵛ] = isone(W[iᵛ]) ? relatedness(m, V[iᵛ], V[i], s) : -Inf end
    X[i] = Inf
    # Step 4: Remove at least q customers from the most related vehicles to this pivot vehicle
    n = 0
    while n < q
        iᵛ = argmax(X)
        v  = V[iᵛ]
        d  = D[v.iᵈ] 
        for r ∈ v.R
            if n ≥ q break end
            while true
                if !isopt(r) || isdormant(r) break end
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        X[iᵛ] = -Inf
        W[iᵛ] = 0
    end
    postremove!(s)
    # Step 5: Return solution
    return s
end



"""
    worstvehicle!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing at least `q` customer 
nodes from routes of low-utilization vehicles.
"""
function worstvehicle!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    V = [v for d ∈ D for v ∈ d.V]
    X = fill(Inf, eachindex(V))     # X[iᵛ]: utilization of vehicle V[iᵛ]
    W = isopt.(V)                   # W[iᵛ]: selection weight for vehicle V[iᵛ]
    # Step 2: Evaluate utilization for each vehicle
    for (iᵛ,v) ∈ pairs(V) X[iᵛ] = isone(W[iᵛ]) ? v.q/(length(v.R) * v.qᵛ) : Inf end
    # Step 3: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        iᵛ = argmin(X)
        v  = V[iᵛ]
        d  = D[v.iᵈ]
        for r ∈ v.R
            if n ≥ q break end
            while true
                if !isopt(r) || isdormant(r) break end
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        X[iᵛ] = Inf
        W[iᵛ] = 0
    end
    postremove!(s)
    # Step 4: Return solution
    return s
end



# -------------------------------------------------- DEPOT REMOVAL --------------------------------------------------
"""
    randomdepot!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after iteratively selecting a random depot node and 
removing customer nodes from its routes until at least `q` customer nodes 
are removed.
"""
function randomdepot!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    W = isopt.(D)                   # W[iᵈ]: selection weight for depot node D[iᵈ]
    # Step 2: Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        iᵈ = sample(rng, eachindex(D), Weights(W))
        d  = D[iᵈ]
        for v ∈ d.V
            if n ≥ q break end
            for r ∈ v.R
                while true
                    if !isopt(r) || isdormant(r) break end
                    nᵗ = d
                    c  = C[r.iˢ]
                    nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                    removenode!(c, nᵗ, nʰ, r, s)
                    n += 1
                    if isequal(nʰ, d) break end
                end
            end
        end
        W[iᵈ] = 0
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



"""
    relateddepot!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing at least `q` customer nodes 
from the routes of the depots most related to a randomly selected 
pivot depot node.
"""
function relateddepot!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    X = fill(-Inf, eachindex(D))    # X[iᵛ]: relatedness of depot node D[iⁿ] with pivot depot node D[i]
    W = isclose.(D)                 # W[iᵈ]: selection weight for depot node D[iᵈ]
    # Step 2: Select a random closed depot node
    i = sample(rng, eachindex(D), Weights(W))
    # Step 3: Evaluate relatedness of this depot node to every depot node
    m = sample(rng, s.φ ? [:q, :l, :t] : [:q, :l])
    for iᵈ ∈ eachindex(D) X[iᵈ] = iszero(W[iᵈ]) ? relatedness(m, D[iᵈ], D[i], s) : -Inf end
    X[i] = Inf
    # Step 4: Remove at least q customer nodes most related to this pivot depot node
    n = 0
    while n < q
        iᵈ = argmax(X)
        d  = D[iᵈ]
        for v ∈ d.V
            if n ≥ q break end
            for r ∈ v.R
                while true
                    if !isopt(r) || isdormant(r) break end
                    nᵗ = d
                    c  = C[r.iˢ]
                    nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                    removenode!(c, nᵗ, nʰ, r, s)
                    n += 1
                    if isequal(nʰ, d) break end
                end
            end
        end
        X[iᵈ] = -Inf
        W[iᵈ] = 0
    end
    postremove!(s)
    # Step 5: Return solution
    return s
end



"""
    worstdepot!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing at least `q` customer 
nodes from routes of low-utilization depot nodes.
"""
function worstdepot!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    X = fill(Inf, eachindex(D))     # X[iᵈ]: utilization of vehicle D[iᵈ]
    W = isopt.(D)                   # W[iᵈ]: selection weight for vehicle D[iᵈ]
    # Step 2: Evaluate utilization for each depot
    for (iᵈ,d) ∈ pairs(D) X[iᵈ] = isone(W[iᵈ]) ? d.q/d.qᵈ : Inf end
    # Step 3: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        iᵈ = argmin(X)
        d  = D[iᵈ]
        for v ∈ d.V
            if n ≥ q break end
            for r ∈ v.R
                while true
                    if !isopt(r) || isdormant(r) break end
                    nᵗ = d
                    c  = C[r.iˢ]
                    nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                    removenode!(c, nᵗ, nʰ, r, s)
                    n += 1
                    if isequal(nʰ, d) break end
                end
            end
        end
        X[iᵈ] = Inf
        W[iᵈ] = 0
    end
    postremove!(s)
    # Step 4: Return solution
    return s
end