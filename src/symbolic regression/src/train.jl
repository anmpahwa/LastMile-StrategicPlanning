# Train and Test Model

using CSV
using DataFrames
using MLJ
using SymbolicRegression

dir = joinpath(dirname(pwd()), "data/instances")
df  = CSV.read("$dir/train/data.csv", DataFrame)

#  Train
## Format 1
X = (;  ρ  = df[:,:pu],
        N  = df[:,:n],
        ⎷δ = sqrt.(df[:,:n] ./ (df[:,:lx] .* df[:,:ly] .* df[:,:r])),
        m  = df[:,:n] .* df[:,:ud] ./ df[:,:qv],
        μˢ = df[:,:us] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        σˢ = df[:,:ss] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        γˢ = df[:,:gs] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        κˢ = df[:,:ks] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]),
        μᵈ = df[:,:ud] ./ df[:,:qv],
        σᵈ = df[:,:sd] ./ df[:,:qv],
        γᵈ = df[:,:gd] ./ df[:,:qv],
        κᵈ = df[:,:kd] ./ df[:,:qv]
    )
y = df[:,:z]
structure = TemplateStructure{(:k, )}(
    ((; k, ), (ρ, N, ⎷δ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ)) -> 2 * ρ * m + k(μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ) * N / ⎷δ
    )
model = SRRegressor(
            niterations=500,
            populations=20,
            binary_operators=[+, -, *, /],
            should_optimize_constants=true,
            should_simplify=true,
            expression_type=TemplateExpression,
            expression_options=(; structure)
        )
r = report(fit!(machine(model, X, y)))
r.equation_strings[r.best_idx]

## Format 2
X = (;  ρ  = df[:,:pu],
        N  = df[:,:n],
        ⎷δ = sqrt.(df[:,:n] ./ (df[:,:lx] .* df[:,:ly] .* df[:,:r])),
        m  = df[:,:m],
        μˢ = df[:,:us] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        σˢ = df[:,:ss] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        γˢ = df[:,:gs] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        κˢ = df[:,:ks] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]),
        μᵈ = df[:,:ud] ./ df[:,:qv],
        σᵈ = df[:,:sd] ./ df[:,:qv],
        γᵈ = df[:,:gd] ./ df[:,:qv],
        κᵈ = df[:,:kd] ./ df[:,:qv]
    )
y = df[:,:z]
structure = TemplateStructure{(:k, )}(
    ((; k, ), (ρ, N, ⎷δ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ)) -> 2 * ρ * m + k(μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ) * N / ⎷δ
    )
model = SRRegressor(
            niterations=500,
            populations=20,
            binary_operators=[+, -, *, /],
            should_optimize_constants=true,
            should_simplify=true,
            expression_type=TemplateExpression,
            expression_options=(; structure)
        )
r = report(fit!(machine(model, X, y)))
r.equation_strings[r.best_idx]

## Format 3
X = (;  ρ  = df[:,:pw],
        N  = df[:,:n],
        ⎷δ = sqrt.(df[:,:n] ./ (df[:,:lx] .* df[:,:ly] .* df[:,:r])),
        m  = df[:,:n] .* df[:,:ud] ./ df[:,:qv],
        μˢ = df[:,:us] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        σˢ = df[:,:ss] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        γˢ = df[:,:gs] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        κˢ = df[:,:ks] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]),
        μᵈ = df[:,:ud] ./ df[:,:qv],
        σᵈ = df[:,:sd] ./ df[:,:qv],
        γᵈ = df[:,:gd] ./ df[:,:qv],
        κᵈ = df[:,:kd] ./ df[:,:qv]
    )
y = df[:,:z]
structure = TemplateStructure{(:k, )}(
    ((; k, ), (ρ, N, ⎷δ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ)) -> 2 * ρ * m + k(μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ) * N / ⎷δ
    )
model = SRRegressor(
            niterations=500,
            populations=20,
            binary_operators=[+, -, *, /],
            should_optimize_constants=true,
            should_simplify=true,
            expression_type=TemplateExpression,
            expression_options=(; structure)
        )
r = report(fit!(machine(model, X, y)))
r.equation_strings[r.best_idx]

## Format 4
X = (;  ρ  = df[:,:pw],
        N  = df[:,:n],
        ⎷δ = sqrt.(df[:,:n] ./ (df[:,:lx] .* df[:,:ly] .* df[:,:r])),
        m  = df[:,:m],
        μˢ = df[:,:us] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        σˢ = df[:,:ss] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        γˢ = df[:,:gs] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]), 
        κˢ = df[:,:ks] ./ sqrt.(df[:,:lx] .* df[:,:ly] .* df[:,:r]),
        μᵈ = df[:,:ud] ./ df[:,:qv],
        σᵈ = df[:,:sd] ./ df[:,:qv],
        γᵈ = df[:,:gd] ./ df[:,:qv],
        κᵈ = df[:,:kd] ./ df[:,:qv]
    )
y = df[:,:z]
structure = TemplateStructure{(:k, )}(
    ((; k, ), (ρ, N, ⎷δ, m, μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ)) -> 2 * ρ * m + k(μˢ, σˢ, γˢ, κˢ, μᵈ, σᵈ, γᵈ, κᵈ) * N / ⎷δ
    )
model = SRRegressor(
            niterations=500,
            populations=20,
            binary_operators=[+, -, *, /],
            should_optimize_constants=true,
            should_simplify=true,
            expression_type=TemplateExpression,
            expression_options=(; structure)
        )
r = report(fit!(machine(model, X, y)))
r.equation_strings[r.best_idx]