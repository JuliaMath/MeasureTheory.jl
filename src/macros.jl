
using MLStyle
using StatsFuns

export @measure

# A fold over ASTs. Example usage in `replaceX`
function foldast(leaf, branch; kwargs...)
    function go(ast)
        MLStyle.@match ast begin
            Expr(head, args...) => branch(head, map(go, args); kwargs...)
            x                   => leaf(x; kwargs...)
        end
    end

    return go
end

# Walk through expr, replacing every occurrence of :X with newX
function replaceX(newX, expr)
    leaf(s) = (s == :X) ? newX : s 
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
        base = replaceX(:(eltype($μ{P})), base)
        q = quote
            struct $μ{P,X} <: MeasureTheory.AbstractMeasure{X}
                par :: P    
            end
        
            function $μ(nt::NamedTuple)
                P = typeof(nt)
                return $μ{P, eltype($μ{P})}(nt)
            end
        
            $μ(;kwargs...) = $μ((;kwargs...))
        
            baseMeasure(μ::$μ{P,X}) where {P,X} = $base
        
            $μ($(p...)) = $μ(;$(p...))
        
            ((::$μ{P} ≪ ::typeof($base) ) where {P})  = true
        end    
    
        if rel == (:≃)
            push!(q.args, :((::typeof($base) ≪ ::$μ{P}) where {P} = true))
        end
    
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
        struct Normal{P, X} <: AbstractMeasure{X}
            par::P
        end
        function Normal(nt::NamedTuple)
            P = typeof(nt)
            return Normal{P, eltype(Normal{P})}
        end
        Normal(; kwargs...) = Normal((; kwargs...))
        (baseMeasure(μ::Normal{P, X}) where {P, X}) = Lebesgue{X}
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
