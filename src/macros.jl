
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
        
            $μ($(p...)) = $μ(;$p)
        
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
"""
macro measure(expr)
    esc(_measure(expr))
end


MacroTools.prettify(@macroexpand @measure Normal(μ,σ) ≃ Lebesgue{X})
