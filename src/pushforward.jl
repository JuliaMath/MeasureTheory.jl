struct PushForward{F,M}
    f::F
    μ::M
end

import Base

Base.broadcast(f::Function, μ::Measure) = PushForward(f,μ)
