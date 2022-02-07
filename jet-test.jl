using MeasureTheory

dist = For(3) do j Normal(σ=j) end

b = basemeasure_depth(dist)