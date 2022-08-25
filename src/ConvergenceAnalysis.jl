module ConvergenceAnalysis

using ImageTransformations
using Interpolations
using Scapin

function refine_and_apply(F_N, x)
    N = grid_size(F_N)
    apply(F_N, imresize(x, size(F_N, 2), N..., method = BSpline(Constant())))
end

export refine_and_apply

end # of module
