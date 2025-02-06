using BSplineKit, LaTeXStrings, Plots, DelimitedFiles
import CSV: read as csvread
import CurveFit: linear_fit
import DataFrames: DataFrame
import LsqFit: curve_fit as lsqfit_curve_fit
import LazySets: convex_hull
import StatsBase: coef, mean
using Interpolations
import LinearAlgebra: eigen, cross, diagm, hessenberg, svd
import FractalDimensions: slopefit, linear_region, linreg
using Dierckx

const I = AbstractVector{<:Real}

"Data return"
struct  Data{I}
    Sw::I
    Sn::I
    v::I
    vw::I
    vn::I
    Ca::Float64
end

"Derived quantities"
struct DerivedData{I}
    v_d::I
    vw_hat::I
    vn_hat::I
    vm::I
end


"Constants, DNM"
Base.@kwdef struct PhysicalQuantities{I}
    μ_w::Float64 # Pa*s
    μ_n::Float64
    μ_e::I
    Ca::Float64
    # v ~ mm/s
    # T is in dyne/mm = 0.01 N/m
end


"Constants, relperm"
struct RelpermQuantities
    σ::Float64 # mN/m
    μ_w::Float64 # Pa⋅s
    μ_n::Float64 # Pa⋅s
    K::Float64 # Darcy
    ϕ::Float64 # 1
    ΔP::Float64 #MPa
end

"Results"
struct Results{I}
    vm_linear_fit::I
    vm_Sw_f_line::I
    vm_Sw_f_line_c0::I
    vw_rec_hom::I
    vn_rec_hom::I
    Sw_tilde::I
    Sn_tilde::I
    k::I
    # Parameters
    a::Float64
    b::Float64
    c::Float64
end

"Read input data from csv-file."
function csv_read(file_string)
    return csvread(file_string, DataFrame; delim=" ", ignorerepeated=true, select=[1, 2, 3, 4, 13], types=Dict(1 => String, 2 => String, 3 => String, 4 => String, 13 => String), skipto=2, header=false)
end

"Find roots of ax^2 + bx + c. Found here: http://mth229.github.io/zeros.html"
function quadratic(a, b, c)
    discr = b^2 - 4 * a * c
    # println("Disc. =", discr)
    sq = (discr > 0) ? sqrt(discr) : sqrt(discr + 0im)
    [(-b - sq) / (2a), (-b + sq) / (2a)]
end

"Computes the tangent-line at index ind."
line(ind) = (Sw .- Sw[ind]) .* v_d[ind] .+ v[ind]

"The gradient line at a point p, computed from a spline-function f at p."
gradient_line(f, p) = (x -> (x - p) * map(diff(f, Derivative(1)), p) + f(p))

"Returns the value of the v-spline at a point x."
spl_v_f(x) = map(spl_v, x) # Helper function

"Converts from Sw/Sn to Sw (cross ratio permutation)"
to_Sw(x) = 1 - (1 - x) .^ -1

"Converts from Sw to Sw/Sn. Note that this is a cross-ratio permutation with an added sign!"
to_ratio(x) = ((1 - x).^-1 - 1) # Adds minus sign to get Sw./Sn!

"Computes cross ratio from affine coordinates x_i."
cross_ratio(x_1, x_2, x_3, x_4) = ((x_3 - x_1) * (x_4 - x_2)) / ((x_3 - x_2) * (x_4 - x_1)) 

"Linear fractional transformation in x, with c being the coefficient vector."
lin_frac(x, c) = (x*c[1] + c[2])/(x*c[3] + c[4])

" A 1D homography, dehomogenized."
proj_dehom(x, c_mat) = (c_mat * x)[1] ./ (c_mat * x)[2]

"Fixed points, Sw"
foo(x1, x2) = if (isreal(x1))
    x1_Sw, x2_Sw = to_Sw.([x1, x2])
else
    x1_Sw, x2_Sw = [0,0]
end

# Other convenient functions
rat_cr(x) = -x/(1 - x) # Different c.r.
conformal_ratio(x_1, x_2, x_3, x_4) = (x_3 - x_1)^2 * (x_4 - x_2)^2 / ((x_3 - x_2)^2 * (x_4 - x_1)^2) # Conformal ratio with affine coordinates.
