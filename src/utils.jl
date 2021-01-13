const EmptyNamedTuple = NamedTuple{(),Tuple{}}

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

using SymbolicUtils: Symbolic

const MaybeSym{T} = Union{T, Symbolic{T}}
