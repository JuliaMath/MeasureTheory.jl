export basemeasure

"""
    basemeasure(μ)

Many measures are defined in terms of a logdensity relative to some base
measure. This makes it important to be able to find that base measure.

For measures not defined in this way, we'll typically have `basemeasure(μ) == μ`.
"""
function basemeasure end

include("basemeasures/trivial.jl")
include("basemeasures/lebesgue.jl")
include("basemeasures/counting.jl")



export isprimitive

"""
    isprimitive(μ)

Most measures are defined in terms of other measures, for example using a
density or a pushforward. Those that are not are considered (in this library,
it's not a general measure theory thing) to be _primitive_. The canonical
example of a primitive measure is `Lebesgue(X)` for some `X`.

The default method is
    isprimitive(μ) = false

So when adding a new primitive measure, it's necessary to add a method for its type
that returns `true`.
"""
isprimitive(μ) = false
