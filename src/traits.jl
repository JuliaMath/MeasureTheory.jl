export IsMeasure

@trait IsMeasure{M,X} where {X = eltype(M)}

###################################################

export HasDensity
export logdensity

@trait IsMeasure{M, X} >: HasDensity{M, X} where {X = eltype(M)} begin
    logdensity :: [M, X] => Real
end

###################################################

export IsMeasurable
export measure

@trait IsMeasure{M,X} >: IsMeasurable{M,S,X} where {X = eltype(S)} begin
    measure :: [M, S] => Real
end

###################################################

export HasRand

@trait IsMeasure{M,X} >: HasRand{M,X} where {X = eltype(M)} begin
    rand :: [M] => eltype(M)
end

###################################################

# @trait Measure{M,X} >: FiniteMeasure{M,X} where {X = eltype{M}} begin

# end
