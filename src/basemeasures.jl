export basemeasure

"""
    basemeasure(μ)


"""
function basemeasure end

include("basemeasures/trivial.jl")
include("basemeasures/lebesgue.jl")



export isprimitive

"""
    isprimitive(μ)

Most measures are defined in terms of other measures, for example using a
density or a pushforward. Those that are not are considered (in this library,
it's not a general measure theory thing) to be _primitive_. The canonical
example of a primitive measure is `Lebesgue{X}()` for some `X`.

The default method is
    isprimitive(μ) = false

So when adding a new primitive measure, it's necessary to add a method for its type
that returns `true`.
"""
isprimitive(μ) = false
