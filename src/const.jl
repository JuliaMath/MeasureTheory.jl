struct Const{T} <: Function
    t::T
end

(c::Const)(x) = c.t
