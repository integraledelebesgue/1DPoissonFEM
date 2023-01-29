using Plots: plot, plot!

const G = 6.67e-11  # Gravitational proportionality constant can be ommitted (set to 1.0) to see a paraboloidal solution shape 

const integral_steps = 1000  # Number of steps in integral() function


# Structure representing n-node 1-D Lagrange symplectic basis
struct Basis

    a::Float64  # Start point
    b::Float64  # End point
 
    n::Integer  # Number of nodes, not added division points. In fact only n-2 points are added.

    h::Float64  # Size of a single element

    basic_funs::Dict{Integer, Function}  # Lazy construction of basic functions
    nodes::Array{Float64, 1}  # All the n nodes, including start and end given

    # Constructor function
    Basis(start_point::Float64, end_point::Float64, n_divisions::Int) = new(
        start_point, 
        end_point, 
        n_divisions, 
        (end_point - start_point) / (n_divisions-1),
        Dict{Integer, Function}(),
        range(start_point, end_point, n_divisions)
    )

end


# Convention: Basis[i] is an i-th basic function
function Base.getindex(basis::Basis, i::Int64)::Function

    if get(basis.basic_funs, i, nothing) === nothing
        basis.basic_funs[i] = basic_i(basis, i)        
    end

    return basis.basic_funs[i]

end


# Technical safety improvement
function safeget(array::Array{Float64, 1}, i::Int64)::Number

    if i ≤ 0
        return array[1] - 1.0
    elseif i ≥ length(array) + 1
        return array[length(array)] + 1.0
    else 
        return array[i]
    end

end


# Create a lambda representing i-th basic function of Lagrange symplectic basis
function basic_i(basis::Basis, i::Int64)::Function

    left = safeget(basis.nodes, i-1)
    right = safeget(basis.nodes, i+1)

    return x -> 
        if left ≤ x ≤ basis.nodes[i]
            return 1/basis.h * (x - left)
        elseif basis.nodes[i] ≤ x ≤ right
            return 1/basis.h * (-x + right)
        else
            return 0
        end
    
end


# Construct an interpolation by taking a dot product of basic functions and node values 
function linear_combination(basis::Basis, coefficients::Array{Float64, 1}, shift::Function, result_length::Int64)::Tuple{Array{Float64, 1}, Array{Float64, 1}}

    domain = range(basis.a, basis.b, result_length)
    image = zeros(result_length)

    for i in eachindex(coefficients)
        image += map(basis[i], domain) .* coefficients[i]
    end

    image += map(shift, domain)

    return domain, image

end


# Arrange lhs bilinear form matrix of Galerkin system
function lhs_matrix(basis::Basis)::Matrix{Float64}

    n = basis.n

    lhs = zeros(n - 2, n - 2)
    diagonal = 2.0 / basis.h
    offset = -1.0 / basis.h

    for i in 2:n-3
        lhs[i, i] = diagonal
        lhs[i, i-1] = offset
        lhs[i, i+1] = offset
    end

    lhs[1, 1] = diagonal
    lhs[1, 2] = offset

    lhs[n-2, n-2] = diagonal
    lhs[n-2, n-3] = offset

    return lhs

end


# Simple rectangular method integration (Riemannian sum)
function integral(start_point::Float64, end_point::Float64, integrand::Function)::Float64

    integral_width = (end_point - start_point) / integral_steps

    return sum(map(integrand, range(start_point, end_point, integral_steps)) .* integral_width)

end


# Arrange rhs functional values of Galerkin system
function rhs_vector(basis::Basis, mass_distribution::Function, shift_slope::Float64)::Array{Float64, 1}

    return [
        (
            -4π * G *integral(basis.nodes[i-1], basis.nodes[i+1], x -> basis[i](x) * mass_distribution(x)) - 
            integral(basis.nodes[i-1], basis.nodes[i+1], x -> (
                    if basis.nodes[i-1] ≤ x ≤ basis.nodes[i]
                        return shift_slope / basis.h
                    elseif basis.nodes[i] < x ≤ basis.nodes[i+1] 
                        return -shift_slope / basis.h
                    else
                        return 0.0
                    end
                )
            )
        ) for i in 2:basis.n-1
    ]

end

# Set up the basis, arrange lhs bilinear form and rhs functional, solve Galerkin equations and construct a dense, plot-ready solution
function solve_FEM(a::Float64, b::Float64, val_a::Float64, val_b::Float64, n_nodes::Int64, mass_distribution::Function, result_length)::Tuple{Array{Float64, 1}, Array{Float64, 1}}

    basis::Basis = Basis(a, b, n_nodes)

    shift_slope::Float64 = (val_b - val_a) / (b - a)

    shift = x -> shift_slope * (x - a) + val_a

    node_values = [0.0; lhs_matrix(basis) \ rhs_vector(basis, mass_distribution, shift_slope); 0.0]

    return linear_combination(basis, node_values, shift, result_length)

end


# Tests used to validate particular functions

function test_basis()

    basis = Basis(0.0, 5.0, 5)
    display(basis.nodes)
    println(basis.h)

end

function test_basic_functions()

    domain_range = (0.0, 5.0)
    n_basis = 5

    basis = Basis(domain_range..., n_basis)
    domain = range(domain_range..., 100)

    plot(legend = false) |> display

    for i in 1:n_basis
        plot!(domain, map(basis[i], domain), seriestype = :scatter) |> display
    end

end

function test_linear_combination()

    basis = Basis(0.0, 5.0, 100)
    coefficients = map(sin, basis.nodes)

    shift = x -> 0.0

    domain, image = linear_combination(basis, coefficients, shift, 1000)

    plot(domain, image) |> display

end

function test_lhs_matrix()

    basis = Basis(0.0, 5.0, 10)
    lhs = lhs_matrix(basis)

    display(lhs)

end

function test_integral()

    fun = x -> x^2

    println(integral(0.0, 1.0, fun))

end


# Main function providing solution to the task
function main(n_nodes::Int64 = 100)

    n_points = 1000

    mass_distribution = x -> (1.0 ≤ x ≤ 2.0 ? 1.0 : 0.0)

    a = 0.0
    b = 3.0
    
    val_a = 5.0
    val_b = 4.0

    domain, image = solve_FEM(a, b, val_a, val_b, n_nodes, mass_distribution, n_points)

    plot(domain, image, legend = false)

end
