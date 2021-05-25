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

export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))
