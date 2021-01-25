export ParameterizedMeasure
abstract type ParameterizedMeasure{N,T} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N,T}, prop::Symbol) where {N,T}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N,T}) where {N,T}
    return N
end

function Base.show(io::IO, μ::ParameterizedMeasure{(),Tuple{}}) 
    print(io, nameof(typeof(μ)), "()")
end

function Base.show(io::IO, μ::ParameterizedMeasure{N,T}) where {N,T}
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)))
    print(io, getfield(μ,:par))
end
