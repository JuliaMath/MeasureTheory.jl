const EmptyNamedTuple = NamedTuple{(),Tuple{}}

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

function fix(f, x)
    y = f(x)
    while x ≠ y
        (x,y) = (y, f(y))
    end

    return y
end

function constructor(::Type{T}) where {T} 
    C = T
    while C isa UnionAll
        C = C.body
    end

    return C.name.wrapper
end

macro trysupport(ex)
    ex = esc(ex)
    quote
        result = nothing
        try
            result = $ex
        catch e
            if e isa DomainError
                result = -Inf
            else
                rethrow()
            end
        end
        result
    end
end

export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))
