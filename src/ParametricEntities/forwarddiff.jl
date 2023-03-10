@info "including forwarddiff.jl from WavePropBase/ParametricEntities"

function jacobian(psurf::ParametricEntity, s::SVector)
    return ForwardDiff.jacobian(psurf.parametrization,s)
end
jacobian(psurf::ParametricEntity, s) = jacobian(psurf, SVector(s))

function jacobian(f, s)
    return ForwardDiff.jacobian(f,s)
end

# higher order derivatives used in some Nystrom methods for lines
function derivative(l::ParametricElement, s)
    @assert domain(l) == ReferenceLine()
    return ForwardDiff.derivative(l, s[1])
end
function derivative2(l::ParametricElement, s)
    return ForwardDiff.derivative(s -> derivative(l, s[1]), s[1])
end
function derivative3(l::ParametricElement, s)
    return ForwardDiff.derivative(s -> derivative2(l, s[1]), s[1])
end
