
using MLStyle
using Random: AbstractRNG
using KeywordCalls
using ConstructionBase
export @parameterized, @half

# A fold over ASTs. Example usage in `replace`
function foldast(leaf, branch; kwargs...)
    function go(ast)
        MLStyle.@match ast begin
            Expr(head, args...) => branch(head, map(go, args); kwargs...)
            x                   => leaf(x; kwargs...)
        end
    end

    return go
end

# function replace(p, f, expr)
#     leaf(s) = p(s) ? f(s) : s
#     branch(head, newargs) = Expr(head, newargs...)
#     foldast(leaf, branch)(expr)
# end

# from https://thautwarm.github.io/MLStyle.jl/latest/tutorials/capture.html
function capture(template, ex, action)
    let template = Expr(:quote, template)
        quote
            @match $ex begin
                $template => $action
                _         => nothing
            end
        end
    end
end

macro capture(template, ex, action)
    capture(template, ex, action) |> esc
end

function _parameterized(__module__, expr)
    @capture ($op($μ($(p...)), $base)) expr begin
        @assert op ∈ [:<<, :≪, :≃]
        μbase = Symbol(:__, μ, :_base)

        μ = esc(μ)
        base = esc(base)

        # @gensym basename
        q = quote
            struct $μ{N,T} <: MeasureBase.ParameterizedMeasure{N}
                par::NamedTuple{N,T}
            end

            const $μbase = $base
            MeasureBase.basemeasure(::$μ) = $μbase
        end

        if !isempty(p)
            # e.g. Normal(μ,σ) = Normal((μ=μ, σ=σ))
            pnames = QuoteNode.(p)
            push!(q.args, :($μ($(p...)) = $μ(NamedTuple{($(pnames...),)}(($(p...),)))))
        end

        return q
    end

    @capture $μ($(p...)) expr begin
        μ = esc(μ)

        q = quote
            struct $μ{N,T} <: MeasureBase.ParameterizedMeasure{N}
                par::NamedTuple{N,T}
            end
        end

        if !isempty(p)
            # e.g. Normal(μ,σ) = Normal((μ=μ, σ=σ))
            pnames = QuoteNode.(p)
            push!(q.args, :($μ($(p...)) = $μ(NamedTuple{($(pnames...),)}(($(p...),)))))
        end

        return q
    end

    @capture $μ($(p...)) expr begin
        μ = esc(μ)

        q = quote
            struct $μ{N,T} <: MeasureBase.ParameterizedMeasure{N}
                par::NamedTuple{N,T}
            end
        end

        if !isempty(p)
            # e.g. Normal(μ,σ) = Normal((μ=μ, σ=σ))
            pnames = QuoteNode.(p)
            push!(q.args, :($μ($(p...)) = $μ(NamedTuple{($(pnames...),)}(($(p...),)))))
        end

        return q
    end
end

"""
    @parameterized <declaration>
    
The <declaration> gives a measure and its default parameters, and specifies
its relation to its base measure. For example,
    @parameterized Normal(μ,σ)
declares the `Normal` is a measure with default parameters `μ and σ`. The result is equivalent to
```
struct Normal{N,T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}
end
KeywordCalls.@kwstruct Normal(μ,σ)
Normal(μ,σ) = Normal((μ=μ, σ=σ))
```
See [KeywordCalls.jl](https://github.com/cscherrer/KeywordCalls.jl) for details on `@kwstruct`.
"""
macro parameterized(expr)
    _parameterized(__module__, expr)
end

"""
    @half dist([paramnames])
Starting from a symmetric univariate measure `dist ≪ Lebesgue(ℝ)`, create a new
measure `Halfdist ≪ Lebesgue(ℝ₊)`. For example,
    @half Normal()
creates `HalfNormal()`, and 
    @half StudentT(ν)
creates `HalfStudentT(ν)`.
"""
macro half(ex)
    _half(__module__, ex)
end

function _half(__module__, ex)
    @match ex begin
        :($dist) => begin
            halfdist = esc(Symbol(:Half, dist))
            dist = esc(dist)
            quote
                $halfdist(args...) = half($dist(args...))
                $halfdist(; kwargs...) = half($dist(; kwargs...))
            end
        end
    end
end