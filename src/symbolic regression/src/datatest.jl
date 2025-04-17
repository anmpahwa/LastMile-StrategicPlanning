# Develop Testing Data

using CSV
using DataFrames
using StatsBase

# Fetch input-output (data) for test instances
function io(ins)
    println(ins)

    dir = joinpath(pwd(), "data/instances/test")
    dfᵖ = CSV.read("$dir/$ins.vrp", DataFrame, silencewarnings=true)
    dfˢ = CSV.read("$dir/$ins.sol", DataFrame, header=false, silencewarnings=true, delim=" ")

    try
        # Metadata Section
        r  = findfirst(contains("DIMENSION"), dfᵖ[:,1])
        n  = parse(Int, strip(dfᵖ[r,2])) - 1
        r  = findfirst(contains("CAPACITY"), dfᵖ[:,1])
        qᵛ = parse(Int, strip(dfᵖ[r,2]))

        # Depot Section
        i = findfirst(contains("DEPOT_SECTION"), dfᵖ[:,1])
        j = findfirst(contains("EOF"), dfᵖ[:,1])
        k = strip(dfᵖ[i+1,1])

        # Node Coordinate Section
        i  = findfirst(contains("NODE_COORD_SECTION"), dfᵖ[:,1])
        j  = findfirst(contains("DEMAND_SECTION"), dfᵖ[:,1])
        ρˣ = [parse(Float64, x) for (i,x,y) ∈ split.(strip.(dfᵖ[(i+1):(j-1),1])) if  isequal(i, k)][1]
        ρʸ = [parse(Float64, y) for (i,x,y) ∈ split.(strip.(dfᵖ[(i+1):(j-1),1])) if  isequal(i, k)][1]
        Xᶜ = [parse(Float64, x) for (i,x,y) ∈ split.(strip.(dfᵖ[(i+1):(j-1),1])) if !isequal(i, k)]
        Yᶜ = [parse(Float64, y) for (i,x,y) ∈ split.(strip.(dfᵖ[(i+1):(j-1),1])) if !isequal(i, k)]
        lˣ = maximum(Xᶜ)-minimum(Xᶜ)
        lʸ = maximum(Yᶜ)-minimum(Yᶜ)
        r  = 1
        ρ̄ᵘ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n])
        K  = [findmin([isequal(i,j) ? Inf : sqrt((Xᶜ[j]-Xᶜ[i])^2 + (Yᶜ[j]-Yᶜ[i])^2) for j ∈ 1:n])[2] for i ∈ 1:n]
        L  = [sqrt((Xᶜ[K[i]]-Xᶜ[i])^2 + (Yᶜ[K[i]]-Yᶜ[i])^2) for i ∈ 1:n]
        μˢ = mean(L)
        σˢ = std(L)
        γˢ = (sum([(L[i] - μˢ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
        κˢ = (sum([(L[i] - μˢ)^4 for i ∈ 1:n]) / n) # ^ (1/4)

        # Demand Section
        i  = findfirst(contains("DEMAND_SECTION"), dfᵖ[:,1])
        j  = findfirst(contains("DEPOT_SECTION"), dfᵖ[:,1])
        Qᶜ = [parse(Float64, q) for (i,q) ∈ split.(strip.(dfᵖ[(i+1):(j-1),1])) if !isone(parse(Int, i))]
        ρ̄ʷ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n], weights(Qᶜ))
        μᵈ = mean(Qᶜ)
        σᵈ = std(Qᶜ)
        γᵈ = (sum([(Qᶜ[i] - μᵈ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
        κᵈ = (sum([(Qᶜ[i] - μᵈ)^4 for i ∈ 1:n]) / n) # ^ (1/4)
        m  = ceil(sum(Qᶜ) / qᵛ)

        # Optimal Solution
        k = parse(Int, k)
        A = [sqrt(((isone(j) ? ρˣ : Xᶜ[j-1]) - (isone(i) ? ρˣ : Xᶜ[i-1]))^2 + ((isone(j) ? ρʸ : Yᶜ[j-1]) - (isone(i) ? ρʸ : Yᶜ[i-1]))^2) for i ∈ 1:(n+1), j ∈ 1:(n+1)]
        z = 0.
        for r ∈ 1:nrow(dfˢ)
            if !contains(dfˢ[r,1], "Route") break end
            i = k
            for c ∈ 3:ncol(dfˢ) 
                if ismissing(dfˢ[r,c]) break end
                j  = dfˢ[r,c] + 1
                z += A[i,j]
                i  = j
            end
            j  = k
            z += A[i,j]
        end

        # Save output
        df = CSV.read("$dir/data.csv", DataFrame)
        push!(df, (nrow(df)+1, ins, n, lˣ, lʸ, ρˣ, ρʸ, qᵛ, r, ρ̄ᵘ, ρ̄ʷ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ, z), promote=true)
        CSV.write("$dir/data.csv", df)
        
    catch
        # Metadata Section
        r  = findfirst(contains("DIMENSION"), skipmissing(dfᵖ[:,1]))
        n  = parse(Int, strip(dfᵖ[r,2])) - 1
        r  = findfirst(contains("CAPACITY"), skipmissing(dfᵖ[:,1]))
        qᵛ = parse(Int, strip(dfᵖ[r,2]))

        # Depot Section
        i = findfirst(contains("DEPOT_SECTION"), skipmissing(dfᵖ[:,1]))
        j = findfirst(contains("EOF"), skipmissing(dfᵖ[:,1]))
        k = parse(Int, strip(dfᵖ[i+1,2]))

        # Node Coordinate Section
        i  = findfirst(contains("NODE_COORD_SECTION"), skipmissing(dfᵖ[:,1]))
        j  = findfirst(contains("DEMAND_SECTION"), skipmissing(dfᵖ[:,1]))
        ρˣ = [parse(Float64, x) for (i,x) ∈ enumerate(dfᵖ[(i+1):(j-1),2]) if  isequal(i, k)][1]
        ρʸ = [y for (i,y) ∈ enumerate(dfᵖ[(i+1):(j-1),3]) if  isequal(i, k)][1]
        Xᶜ = [parse(Float64, x) for (i,x) ∈ enumerate(dfᵖ[(i+1):(j-1),2]) if !isequal(i, k)]
        Yᶜ = [y for (i,y) ∈ enumerate(dfᵖ[(i+1):(j-1),3]) if !isequal(i, k)]
        lˣ = maximum(Xᶜ)-minimum(Xᶜ)
        lʸ = maximum(Yᶜ)-minimum(Yᶜ)
        r  = 1
        ρ̄ᵘ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n])
        K  = [findmin([isequal(i,j) ? Inf : sqrt((Xᶜ[j]-Xᶜ[i])^2 + (Yᶜ[j]-Yᶜ[i])^2) for j ∈ 1:n])[2] for i ∈ 1:n]
        L  = [sqrt((Xᶜ[K[i]]-Xᶜ[i])^2 + (Yᶜ[K[i]]-Yᶜ[i])^2) for i ∈ 1:n]
        μˢ = mean(L)
        σˢ = std(L)
        γˢ = (sum([(L[i] - μˢ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
        κˢ = (sum([(L[i] - μˢ)^4 for i ∈ 1:n]) / n) # ^ (1/4)


        # Demand Section
        i  = findfirst(contains("DEMAND_SECTION"), skipmissing(dfᵖ[:,1]))
        j  = findfirst(contains("DEPOT_SECTION"), skipmissing(dfᵖ[:,1]))
        Qᶜ = [parse(Float64, q) for (i,q) ∈ enumerate(dfᵖ[(i+1):(j-1),2]) if !isone(i)]
        ρ̄ʷ = mean([sqrt((ρˣ-Xᶜ[i])^2 + (ρʸ-Yᶜ[i])^2) for i ∈ 1:n], weights(Qᶜ))
        μᵈ = mean(Qᶜ)
        σᵈ = std(Qᶜ)
        γᵈ = (sum([(Qᶜ[i] - μᵈ)^3 for i ∈ 1:n]) / n) # ^ (1/3)
        κᵈ = (sum([(Qᶜ[i] - μᵈ)^4 for i ∈ 1:n]) / n) # ^ (1/4)
        m  = ceil(sum(Qᶜ) / qᵛ)

        # Optimal Solution
        A = [sqrt(((isone(j) ? ρˣ : Xᶜ[j-1]) - (isone(i) ? ρˣ : Xᶜ[i-1]))^2 + ((isone(j) ? ρʸ : Yᶜ[j-1]) - (isone(i) ? ρʸ : Yᶜ[i-1]))^2) for i ∈ 1:(n+1), j ∈ 1:(n+1)]
        z = 0.
        for r ∈ 1:nrow(dfˢ)
            if !contains(dfˢ[r,1], "Route") break end
            i = k
            for c ∈ 3:ncol(dfˢ) 
                if ismissing(dfˢ[r,c]) break end
                j  = dfˢ[r,c] + 1
                z += A[i,j]
                i  = j
            end
            j  = k
            z += A[i,j]
        end

        # Save output
        df = CSV.read("$dir/data.csv", DataFrame)
        push!(df, (nrow(df)+1, ins, n, lˣ, lʸ, ρˣ, ρʸ, qᵛ, r, ρ̄ᵘ, ρ̄ʷ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ, z), promote=true)
        CSV.write("$dir/data.csv", df)
    end
end

sets = ["set A", "set B", "set C", "set D", "set E", "set F", "set G", "set H", "set I", "set J", "set K", "set L", "set M"]
for set ∈ sets 
    instances = unique!(["$set/$file" for (file,type) ∈ split.(readdir(joinpath(pwd(), "data/instances/test/$set")), ".")])
    for ins ∈ instances io(ins) end
end