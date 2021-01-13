
using MLStyle
using StatsFuns

export @measure

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

function replace(p, f, expr)
    leaf(s) = p(s) ? f(s) : s 
    branch(head, newargs) = Expr(head, newargs...)
    foldast(leaf, branch)(expr)
end

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

function _measure(expr)
    @capture $rel($μ($(p...)), $base) expr begin
        base = replace(s -> s∈p, s -> :(μ.$s), base)

        q = quote
            struct $μ{N,T} <: ParameterizedMeasure{N,T}
                par :: NamedTuple{N,T}
            end

            function MeasureTheory.basemeasure(μ::$μ{P}) where {P}
                return $base
            end
            
            # e.g. Normal(μ,σ) = Normal(;μ=μ, σ=σ)
            # Requires Julia 1.5
        
            $μ($(p...)) = $μ(;$(p...))
            
            # e.g. Normal(;μ=μ, σ=σ) = Normal((μ=μ, σ=σ))
            $μ(;kwargs...) = $μ((;kwargs...))
 
            # ((::$μ{P} ≪ ::typeof($base) ) where {P})  = true
        end    
    
        # if rel == (:≃)
        #     push!(q.args, :((::typeof($base) ≪ ::$μ{P}) where {P} = true))
        # end
    
        return q
    end
end

"""
    @measure <declaration>
    
The <declaration> gives a measure and its default parameters, and specifies
its relation to its base measure. For example,

    @measure Normal(μ,σ) ≃ Lebesgue{X}

declares the `Normal` is a measure with default parameters `μ and σ`, and it is
equivalent to its base measure, which is `Lebesgue{X}`

You can see the generated code like this:

    julia> MacroTools.prettify(@macroexpand @measure Normal(μ,σ) ≃ Lebesgue{X})
    quote
        struct Normal{P, X} <: AbstractMeasure
            par::P
        end
        function Normal(nt::NamedTuple)
            P = typeof(nt)
            return Normal{P, eltype(Normal{P})}
        end
        Normal(; kwargs...) = Normal((; kwargs...))
        (basemeasure(μ::Normal{P, X}) where {P, X}) = Lebesgue{X}
        Normal(μ, σ) = Normal(; Any[:μ, :σ])
        ((:≪)(::Normal{P, X}, ::Lebesgue{X}) where {P, X}) = true
        ((:≪)(::Lebesgue{X}, ::Normal{P, X}) where {P, X}) = true
    end

Note that the `eltype` function needs to be defined separately by the user.
"""
macro measure(expr)
    esc(_measure(expr))
end


# (@macroexpand @measure Normal(μ,σ) ≃ (1/sqrt2π) * Lebesgue(X)) |> MacroTools.prettify


using MLStyle

macro μσ_methods(ex)
    esc(_μσ_methods(ex))
end

function _μσ_methods(ex)
    @match ex begin
        :($dist($(args...))) => begin
            argnames = QuoteNode.(args)


            quote
                function logdensity(d::$dist{($(argnames...), :μ, :σ)}, x::X) where {X} 
                    z = (x - d.μ) / d.σ   
                    return logdensity($ex, z) - log(d.σ)
                end

                function logdensity(d::$dist{($(argnames...), :σ)}, x::X) where {X} 
                    z = x / d.σ   
                    return logdensity($ex, z) - log(d.σ) 
                end

                function logdensity(d::$dist{($(argnames...), :μ)}, x::X) where {X} 
                    z = x - d.μ
                    return logdensity($ex, z)
                end
            end 
        end
    end
end
