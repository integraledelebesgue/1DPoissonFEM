push!(LOAD_PATH, @__DIR__)

struct Basis

    a::Float64
    b::Float64

    n::Integer  # Number of inner points; a = x₀, b = xₙ₊₁

    h::Float64

end


function getindex(basis::Basis, i::Integer)

    return 

end


function ηᵢ(x, i, a, h)

    return ( x -> 
        if a + (i-1)h ≤ x ≤ a + ih
            return (x - a - (i-1)h)
        
        elseif a + ih ≤ x ≤ a + (i+1)h
            return (-x + a + (i+1)h)
        
        else
            return 0
        end
    )

end

