push!(LOAD_PATH, @__DIR__)

using Plots: plot

G = 6.67e-11


function equation_matrix(n::Integer)::Matrix{Float64}  # Builds a FEM matrix for Poisson equation with n internal elements

    special_1::Vector{Float64} = [[-3.0, 0, 1.0] ; [0.0 for _ in 1:n-3]]
    special_2::Vector{Float64} = [[0.0, -2.0, 0.0, 1.0] ; [0.0 for _ in 1:n-4]]
    special_n_1::Vector{Float64} = [[0.0 for _ in 1:n-4] ; [1.0, 0.0, -2.0, 0.0]]
    special_n::Vector{Float64} = [[0.0 for _ in 1:n-3] ; [-1.0, 0.0, -1.0]]
    
    return [
        special_1'
        special_2'

        [
            [[0.0 for _ in 1:k-3]; [1.0, 0.0, -2.0, 0.0, 1.0]; [0.0 for _ in k+3:n]]' for k in 3:n-2
        ]...

        special_n_1'
        special_n'
    ]'

end


function basis(x1, x2, x3)::Function  # returns a lambda which computes tent function values in range {xi-1..xi..xi+1}

    # Boundary cases:
    if x1 === nothing
        return (x -> (x3 - x) / (x3 - x2))

    elseif x3 === nothing
        return (x -> (x - x1) / (x2 - x1))

    # normal case:
    else
        function to_return(x) x ≤ x2 ? 
            ( (x - x1) / (x2 - x1) ) : 
            ( (x3 - x) / (x3 - x2) )
        end
        return to_return
    end
end


function solve_nodes(domain::Array{Float64, 1}, a_value::Float64, b_value::Float64, rhs::Function)::Array{Float64, 1}  # Solves FEM matrix equation having domain (a domain without boundary), Dirichlet boundary conditions and right-hand-side function given

    int_domain = internal_domain(domain)

    n = length(int_domain)

    right_side = int_domain .|> rhs |> collect

    right_side[1] -= a_value
    right_side[2] -= a_value

    right_side[n-1] -= b_value
    right_side[n] -= 2b_value

    return equation_matrix(n) \ right_side

end


function simple_mass_distribution(start::Float64, stop::Float64, mass::Float64)::Function  # returns a function ρ(x) = mass if x ∈ [start, stop], 0 otherwise

    return x -> start ≤ x ≤ stop ? mass : 0.0

end


function build_domain(start::Float64, stop::Float64, n::Integer)::Array{Float64, 1}

    return LinRange(start, stop, n+2) |> collect

end


function internal_domain(domain::Array{Float64, 1})::Array{Float64, 1}  # cuts off start and end of the domain to leave only the divided elements

    return domain[2:length(domain)-1]

end


function build_solution(domain::Array{Float64, 1}, k::Integer, node_values::Array{Float64, 1})::Tuple{Array{Float64, 1}, Array{Float64, 1}}

    n = length(domain)

    points = [collect(LinRange(domain[i], domain[i+1], k)) for i in 1:n-1]

    vals = reduce(
        vcat, 
        [
            (
                map(basis(i-1 > 0 ? domain[i-1] : nothing, domain[i], domain[i+1]), points[i]) .* node_values[i] .+ 
                map(basis(domain[i], domain[i+1], i+2 < n ? domain[i+2] : nothing), points[i]) .* node_values[i+1]
            ) for i in 1:n-1
        ]
    )

    return reduce(vcat, points), vals

end


function main(n::Integer, k::Integer)

    μ = simple_mass_distribution(1.0, 2.0, 5.0)

    domain = build_domain(0.0, 3.0, n)

    node_values = [5.0; solve_nodes(domain, 5.0, 4.0, μ); 4.0]

    points, values = build_solution(domain, k, node_values)

    plot(points, values)

end


main(100, 3)
