module VmProjective
include("source.jl")
include("compute.jl")
include("plotResults.jl")
include("homography.jl")

export computeall, getDNMdata, interpolations_derived, spline_convex, noneuclidean, getrelpermdata, myplotrelperm
export Data, DerivedData, PhysicalQuantities, RelpermQuantities, Results



end
