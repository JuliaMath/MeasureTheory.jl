
using MLStyle
import MacroTools

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
    q = @capture $rel($μ($(p...)), $base) expr begin
        q = quote
            struct $μ{P,X} <: AbstractMeasure{X}
                par :: P    
            end
        
            function $μ(nt::NamedTuple)
                P = typeof(nt)
                return $μ{P, eltype($μ{P})}
            end
        
            $μ(;kwargs...) = $μ((;kwargs...))
        
            baseMeasure(μ::$μ{P,X}) where {P,X} = $base
        
            $μ($(p...)) = $μ(;$(p...))
        
            :≪(::$μ{P,X}, ::$base) where {P,X} = true
        end    
    
        if rel == (:≃)
            push!(q.args, :(:≪(::$base, ::$μ{P,X}) where {P,X} = true))
        end
    
        q
    end

    return q 
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


MacroTools.prettify(@macroexpand @measure Normal(μ,σ) ≃ Lebesgue{X})
