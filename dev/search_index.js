var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MeasureTheory","category":"page"},{"location":"#MeasureTheory","page":"Home","title":"MeasureTheory","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MeasureTheory]","category":"page"},{"location":"#MeasureTheory.Density","page":"Home","title":"MeasureTheory.Density","text":"struct Density{M,B}\n    μ::M\n    base::B\nend\n\nFor measures μ and ν with μ≪ν, the density of μ with respect to ν (also called the Radon-Nikodym derivative dμ/dν) is a function f defined on the support of ν with the property that for any measurable a ⊂ supp(ν), μ(a) = ∫ₐ f dν.\n\nBecause this function is often difficult to express in closed form, there are many different ways of computing it. We therefore provide a formal representation to allow comptuational flexibilty.\n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.DensityMeasure","page":"Home","title":"MeasureTheory.DensityMeasure","text":"struct DensityMeasure{F,B} <: AbstractMeasure\n    density :: F\n    base    :: B\nend\n\nA DensityMeasure is a measure defined by a density with respect to some other \"base\" measure \n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.Kernel","page":"Home","title":"MeasureTheory.Kernel","text":"kernel(f, M)\nkernel((f1, f2, ...), M)\n\nA kernel κ = kernel(f, m) returns a wrapper around a function f giving the parameters for a measure of type M, such that κ(x) = M(f(x)...) respective κ(x) = M(f1(x), f2(x), ...)\n\nIf the argument is a named tuple (;a=f1, b=f1), κ(x) is defined as M(;a=f(x),b=g(x)).\n\nReference\n\nhttps://en.wikipedia.org/wiki/Markov_kernel\n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.PowerMeasure","page":"Home","title":"MeasureTheory.PowerMeasure","text":"struct PowerMeasure{M,N}\n    μ::M\n    size::NTuple{N,Int}\nend\n\nA power measure is a product of a measure with itself. The number of elements in the product determines the dimensionality of the resulting support.\n\nNote that power measures are only well-defined for integer powers.\n\nThe nth power of a measure μ can be written μ^x.\n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.SuperpositionMeasure","page":"Home","title":"MeasureTheory.SuperpositionMeasure","text":"struct SuperpositionMeasure{X,NT} <: AbstractMeasure\n    components :: NT\nend\n\nSuperposition of measures is analogous to mixture distributions, but (because measures need not be normalized) requires no scaling.\n\nThe superposition of two measures μ and ν can be more concisely written as μ + ν.\n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.WeightedMeasure","page":"Home","title":"MeasureTheory.WeightedMeasure","text":"struct WeightedMeasure{R,M} <: AbstractMeasure\n    logweight :: R\n    base :: M\nend\n\n\n\n\n\n","category":"type"},{"location":"#MeasureTheory.:≃-Tuple{Any,Any}","page":"Home","title":"MeasureTheory.:≃","text":"≃(μ,ν)\n\nEquivalence of Measure\n\nMeasures μ and ν on the same space X are equivalent, written μ ≃ ν, if μ ≪ ν and ν ≪ μ.\n\nNote that this is often written ~ in the literature, but this is overloaded in probabilistic programming, so we use alternate notation.\n\n\n\n\n\n","category":"method"},{"location":"#MeasureTheory.:≪","page":"Home","title":"MeasureTheory.:≪","text":"≪(μ,ν)\n\nAbsolute continuity\n\nThe following are equivalent:\n\nμ ≪ ν\nμ is absolutely continuous wrt ν\nThere exists a function f such that μ = ∫f dν\n\n\n\n\n\n","category":"function"},{"location":"#MeasureTheory.basemeasure","page":"Home","title":"MeasureTheory.basemeasure","text":"basemeasure(μ)\n\n\n\n\n\n","category":"function"},{"location":"#MeasureTheory.isprimitive-Tuple{Any}","page":"Home","title":"MeasureTheory.isprimitive","text":"isprimitive(μ)\n\nMost measures are defined in terms of other measures, for example using a density or a pushforward. Those that are not are considered (in this library, it's not a general measure theory thing) to be primitive. The canonical example of a primitive measure is Lebesgue(X) for some X.\n\nThe default method is     isprimitive(μ) = false\n\nSo when adding a new primitive measure, it's necessary to add a method for its type that returns true.\n\n\n\n\n\n","category":"method"},{"location":"#MeasureTheory.logdensity","page":"Home","title":"MeasureTheory.logdensity","text":"logdensity(μ::Measure{X}, x::X)\n\nCompute the logdensity of the measure μ at the point x. This is the standard way to define logdensity for a new measure. the base measure is implicit here, and is understood to be basemeasure(μ).\n\nMethods for computing density relative to other measures will be\n\n\n\n\n\n","category":"function"},{"location":"#MeasureTheory.@measure-Tuple{Any}","page":"Home","title":"MeasureTheory.@measure","text":"@measure <declaration>\n\nThe <declaration> gives a measure and its default parameters, and specifies its relation to its base measure. For example,\n\n@measure Normal(μ,σ) ≃ Lebesgue{X}\n\ndeclares the Normal is a measure with default parameters μ and σ, and it is equivalent to its base measure, which is Lebesgue{X}\n\nYou can see the generated code like this:\n\njulia> MacroTools.prettify(@macroexpand @measure Normal(μ,σ) ≃ Lebesgue{X})\nquote\n    struct Normal{P, X} <: AbstractMeasure\n        par::P\n    end\n    function Normal(nt::NamedTuple)\n        P = typeof(nt)\n        return Normal{P, eltype(Normal{P})}\n    end\n    Normal(; kwargs...) = Normal((; kwargs...))\n    (basemeasure(μ::Normal{P, X}) where {P, X}) = Lebesgue{X}\n    Normal(μ, σ) = Normal(; Any[:μ, :σ])\n    ((:≪)(::Normal{P, X}, ::Lebesgue{X}) where {P, X}) = true\n    ((:≪)(::Lebesgue{X}, ::Normal{P, X}) where {P, X}) = true\nend\n\nNote that the eltype function needs to be defined separately by the user.\n\n\n\n\n\n","category":"macro"}]
}
