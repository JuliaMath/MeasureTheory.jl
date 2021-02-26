const EmptyNamedTuple = NamedTuple{(),Tuple{}}

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

function fix(f, x)
    y = f(x)

    if y == x 
        return y
    else
        return fix(f, y)
    end
end

Dists.logpdf(d::AbstractMeasure, x) = logdensity(d,x)
