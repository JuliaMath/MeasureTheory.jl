struct Const{T} <: Function
    t::T
end

function Base.show(io::IO, ::MIME"text/plain", c::Const)
    io = IOContext(io, :compact => true)
    print(io, "_ -> ")
    print(io, c.t)
end

(c::Const)(x) = c.t
