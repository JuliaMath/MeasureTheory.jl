struct Exp{X} <: AbstractFloat
    log::X
end

Base.iszero(::Exp) = false

Base.log(x::Exp) = x.log

function Base.show(io::IO, x::Exp)
    io = IOContext(io, :compact => true)
    print(io, "Exp(", log(x), ")")
end

Base.:*(a::Exp, b::Exp) = Exp(log(a) + log(b))

Base.:/(a::Number, b::Exp) = a * Exp(-log(b))

Base.inv(a::Exp) = Exp(-log(a))

Base.:^(a::Exp, b::Number) = Exp(b * log(a))

Base.:^(a::Exp, b::Integer) = Exp(b * log(a))

Base.:^(a::Exp, b) = Exp(b * log(a))

Base.promote_rule(T::Type{<:Real}, ::Type{Exp}) = T
Base.promote_rule(::Type{Exp}, T::Type{<:Real}) = T

Base.convert(::Type{<:Real}, x::Exp) = exp(log(x))



# Exp(3) / Exp(2)

# Exp(3) * Exp(2)

# 3 + Exp(2)

# 3 / Exp(4)

# Exp(2) ^ 3

# 1 / Exp(2)

# 2 + Exp(3) ^ 0
