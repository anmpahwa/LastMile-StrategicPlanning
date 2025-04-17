using Random
using DataFrames
using CSV

# Fetch input-output (data) for test instances
function io(ins, seed)
    println(ins)

    dir = joinpath(dirname(pwd()), "MetaHeuristic/instances")
    df  = CSV.read("$dir/$ins.vrp", DataFrame, silencewarnings=true)
    
    mkdir(joinpath(dirname(pwd()), "MetaHeuristic/instances/$ins"))

    try
        # Metadata Section
        r  = findfirst(contains("DIMENSION"), df[:,1])
        n  = parse(Int, strip(df[r,2])) - 1
        r  = findfirst(contains("CAPACITY"), df[:,1])
        qᵛ = parse(Int, strip(df[r,2]))

        # Depot Section
        i = findfirst(contains("DEPOT_SECTION"), df[:,1])
        j = findfirst(contains("EOF"), df[:,1])
        k = strip(df[i+1,1])

        # Node Coordinate Section
        i  = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
        j  = findfirst(contains("DEMAND_SECTION"), df[:,1])
        ρˣ = [parse(Float64, x) for (i,x,y) ∈ split.(strip.(df[(i+1):(j-1),1])) if  isequal(i, k)][1]
        ρʸ = [parse(Float64, y) for (i,x,y) ∈ split.(strip.(df[(i+1):(j-1),1])) if  isequal(i, k)][1]
        Xᶜ = [parse(Float64, x) for (i,x,y) ∈ split.(strip.(df[(i+1):(j-1),1])) if !isequal(i, k)]
        Yᶜ = [parse(Float64, y) for (i,x,y) ∈ split.(strip.(df[(i+1):(j-1),1])) if !isequal(i, k)]
        lˣ = maximum(Xᶜ)-minimum(Xᶜ)
        lʸ = maximum(Yᶜ)-minimum(Yᶜ)
        
        # Demand Section
        i  = findfirst(contains("DEMAND_SECTION"), df[:,1])
        j  = findfirst(contains("DEPOT_SECTION"), df[:,1])
        Qᶜ = [parse(Float64, q) for (i,q) ∈ split.(strip.(df[(i+1):(j-1),1])) if !isone(parse(Int, i))]

        # Arcs
        A = [sqrt(((isone(j) ? ρˣ : Xᶜ[j-1]) - (isone(i) ? ρˣ : Xᶜ[i-1]))^2 + ((isone(j) ? ρʸ : Yᶜ[j-1]) - (isone(i) ? ρʸ : Yᶜ[i-1]))^2) for i ∈ 1:(n+1), j ∈ 1:(n+1)]

        # Save Files
        dfᵈ = DataFrame(in=[1], x=[ρˣ], y=[ρʸ], qd=[sum(Qᶜ)], ts=[0], te=[0], co=[0], cf=[0], phi=[0])
        dfᶜ = DataFrame(in=2:(n+1), x=Xᶜ, y=Yᶜ, qc=Qᶜ, tc=zeros(n), tr=zeros(n), te=zeros(n), tl=zeros(n))
        dfᵛ = DataFrame(iv=[1], jv=[1], id=[1], qv=[qᵛ], lv=[2n*(lˣ+lʸ)], sv=[1], tf=[0], td=[0], tc=[0], tw=[0], r=[n], cd=[1], ct=[0], cf=[0])
        dfᵃ = DataFrame(A, :auto)
        CSV.write("$dir/$ins/depot_nodes.csv", dfᵈ)
        CSV.write("$dir/$ins/customer_nodes.csv", dfᶜ)
        CSV.write("$dir/$ins/vehicles.csv", dfᵛ)
        CSV.write("$dir/$ins/arcs.csv", dfᵃ)
        
    catch
        # Metadata Section
        r  = findfirst(contains("DIMENSION"), skipmissing(df[:,1]))
        n  = parse(Int, strip(df[r,2])) - 1
        r  = findfirst(contains("CAPACITY"), skipmissing(df[:,1]))
        qᵛ = parse(Int, strip(df[r,2]))

        # Depot Section
        i = findfirst(contains("DEPOT_SECTION"), skipmissing(df[:,1]))
        j = findfirst(contains("EOF"), skipmissing(df[:,1]))
        k = parse(Int, strip(df[i+1,2]))

        # Node Coordinate Section
        i  = findfirst(contains("NODE_COORD_SECTION"), skipmissing(df[:,1]))
        j  = findfirst(contains("DEMAND_SECTION"), skipmissing(df[:,1]))
        ρˣ = [parse(Float64, x) for (i,x) ∈ enumerate(df[(i+1):(j-1),2]) if  isequal(i, k)][1]
        ρʸ = [y for (i,y) ∈ enumerate(df[(i+1):(j-1),3]) if  isequal(i, k)][1]
        Xᶜ = [parse(Float64, x) for (i,x) ∈ enumerate(df[(i+1):(j-1),2]) if !isequal(i, k)]
        Yᶜ = [y for (i,y) ∈ enumerate(df[(i+1):(j-1),3]) if !isequal(i, k)]
        lˣ = maximum(Xᶜ)-minimum(Xᶜ)
        lʸ = maximum(Yᶜ)-minimum(Yᶜ)

        # Demand Section
        i  = findfirst(contains("DEMAND_SECTION"), skipmissing(df[:,1]))
        j  = findfirst(contains("DEPOT_SECTION"), skipmissing(df[:,1]))
        Qᶜ = [parse(Float64, q) for (i,q) ∈ enumerate(df[(i+1):(j-1),2]) if !isone(i)]

        # Arcs
        A = [sqrt(((isone(j) ? ρˣ : Xᶜ[j-1]) - (isone(i) ? ρˣ : Xᶜ[i-1]))^2 + ((isone(j) ? ρʸ : Yᶜ[j-1]) - (isone(i) ? ρʸ : Yᶜ[i-1]))^2) for i ∈ 1:(n+1), j ∈ 1:(n+1)]
        
        # Save Files
        dfᵈ = DataFrame(in=[1], x=[ρˣ], y=[ρʸ], qd=[sum(Qᶜ)], ts=[0], te=[0], co=[0], cf=[0], phi=[0])
        dfᶜ = DataFrame(in=2:(n+1), x=Xᶜ, y=Yᶜ, qc=Qᶜ, tc=zeros(n), tr=zeros(n), te=zeros(n), tl=zeros(n))
        dfᵛ = DataFrame(iv=[1], jv=[1], id=[1], qv=[qᵛ], lv=[2n*(lˣ+lʸ)], sv=[1], tf=[0], td=[0], tc=[0], tw=[0], r=[n], cd=[1], ct=[0], cf=[0])
        dfᵃ = DataFrame(A, :auto)
        CSV.write("$dir/$ins/depot_nodes.csv", dfᵈ)
        CSV.write("$dir/$ins/customer_nodes.csv", dfᶜ)
        CSV.write("$dir/$ins/vehicles.csv", dfᵛ)
        CSV.write("$dir/$ins/arcs.csv", dfᵃ)
        
    end

end

#instances = unique!(["$file" for (file,type) ∈ split.(readdir(joinpath(dirname(pwd()), "MetaHeuristic/instances")), ".")])
instances = ["Golden_8", "X-n459-k26", "Golden_16", "X-n502-k39"]
for ins ∈ instances io(ins, 1104) end