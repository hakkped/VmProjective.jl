"""
Get data for the dynamic network model.
"""
function getDNMdata(MTP::Vector{String})
    # These combinations do not exist
    if  MTP[2] == "4.0" && MTP[3] == "0.20"
        return
    end
    relative_2d =  "../data/2D_PB_CP_Eqn/vComov-Sw_M-"*MTP[1]*"_T-"*MTP[2]*"_P-"*MTP[3]*".d"
    file2D = joinpath(@__DIR__,relative_2d)
    file = file2D # Set desired dimension
    # String manipulation, extract M,T,P
    file_string = replace(split(file, "..")[end], ".d" => "" )
    findvalue(name, length ) = first((split( split(file_string, "..")[end],  name*"-"))[2], length)
    # M,T,P =  parse.(Float64, findvalue.(["M", "T", "P"], [4,3,4]))
    M,T,P = parse.(Float64, [MTP[1], MTP[2], MTP[3]])

    fv(file, column) = parse.(Float64, (csv_read(file))[:, column])[1:end] # Read from file
    Sw, v, vw, vn, Ca = fv.(Ref(file), collect(1:5)) # Create vectors
    Sn = 1 .- fv(file, 1) # Define Sn
    return Data(Sw, Sn, v, vw, vn, Ca)
end

"""
Get relperm data.
Fewer datasets.
"""
function getrelpermdata(file::String, values::RelpermQuantities)
    (; σ,μ_w,μ_n,K,ϕ,ΔP)  = values
    darcy = 9.869*10^-13
    M = readdlm(file; comments = true)

    # NOTE: check data to be sure what columns are which data! 
    # The order might have changed
    krw = float.(M[:, 2])
    krn = float.(M[:, 3])
    Sn = float.(M[:, 1])
    Sw = 1 .- Sn

    K = K*darcy
    ΔP = ΔP*10^6 # To Pa
    visc_conv = 10^-3 # From cP (mPa ⋅ s) to Pa⋅s
    vw_r = @. -(K*krw*ΔP)/(Sw*ϕ*μ_w*visc_conv)
    vn_r = @. -(K*krn*ΔP)/(Sn*ϕ*μ_n*visc_conv)
    ind = 1
    # Manual check if entries are NaN or 0.
    if isinf(vn_r[1]) || isinf(vw_r[1]) || isnan(vw_r[1]) || isnan(vn_r[1])
        ind = 2
    end
    if iszero(vw_r[end]) || iszero(vn_r[end])
        vw_r = vw_r[1:end-1]
        vn_r = vn_r[1:end-1]
        Sw = Sw[1:end-1]
        Sn = Sn[1:end-1]
        krw = krw[1:end-1]
        krn = krn[1:end-1]
    end
    # Sort in same order, with correct indices
    ord = sortperm(Sw[ind:end])
    Sw = sort(Sw[ind:end])
    Sn = (Sn[ind:end])[ord]
    vw_r = (vw_r[ind:end])[ord]
    vn_r = (vn_r[ind:end])[ord]
    krw = (krw[ind:end])[ord]
    krn = (krn[ind:end])[ord]

    μ_e = @. (Sw * μ_w + Sn * μ_n)
    v = @. Sw * vw_r + Sn * vn_r
    Ca = sum((μ_e .* v) ./ (σ * 10^-3)) ./ length(Sw)
    return Data(Sw, Sn, v, vw_r, vn_r, Ca), krw, krn
end

"""
  Compute interpolated and derived data.
"""
function interpolations_derived(data::Data)
    (; Sw, Sn, v, vw, vn, Ca) = data
    ord = sortperm(Sw)
    Sw = sort(Sw)
    Sn = Sn[ord]
    vw = vw[ord]
    vn = vn[ord]
    v = v[ord]

    # Some examples of interpolations

    # BSplineKit. Gives approximately same result as with Interpolations.
    # spl_v = BSplineKit.interpolate(Sw, v, BSplineOrder(4), BSplineKit.Natural()) 
    # spl_v = BSplineKit.interpolate(Sw, v, BSplineOrder(4))
    # spl_v = BSplineKit.extrapolate(spl_v, BSplineKit.Linear()) 
    # v_d = BSplineKit.diff(spl_v, Derivative(1))
    # v_d = v_d.(Sw)

    # Interpolations
    Interpolations.deduplicate_knots!(Sw)
    spl_v = Interpolations.linear_interpolation((Sw,), v) # Only uses linear interpolation
    # spl_v = Interpolations.interpolate((Sw,), v, Gridded(Interpolations.Linear())) # Only uses linear interpolation
    v_d = only.(Interpolations.gradient.(Ref(spl_v), Sw))

    # Dierckx
    # spl_v = Spline1D(Sw, v[ord])
    # v_d = Dierckx.derivative(spl_v, Sw)

    vn_hat = @. v - Sw * v_d
    vw_hat = @. v + Sn * v_d
    vm = @. v_d - vw + vn
    return DerivedData(v_d, vw_hat, vn_hat, vm)
end

"""
Convex envelopes, more spline functions etc.
Not used. 
"""
function spline_convex(data::Data)
    (; Sw, Sn, v, vw, vn, Ca) = data
    spl_v, spl_vw, spl_vn, spl_vn_min_vw = BSplineKit.interpolate.(Ref(Sw), [v, vw, vn, vn .- vw], Ref(BSplineOrder(4)), Ref(BSplineKit.Natural()))
    spl_v = BSplineKit.extrapolate(spl_v, BSplineKit.Linear())

    # # # # # # # # # #  Spline derivatives # # # # # # # # # #
    v_d_spl, v_d_w_spl, v_d_n_spl = BSplineKit.diff.([spl_v, spl_vw, spl_vn], Ref(Derivative(1)))
    v_dd = BSplineKit.diff(spl_v, Derivative(2))
    v_d, vw_spl, vn_spl, v_d_w, v_d_n = (v_d_spl.(Sw), spl_vw.(Sw), spl_vn.(Sw), v_d_w_spl.(Sw), v_d_n_spl.(Sw))

    # # # # # # # # # # Convex hull of v # # # # # # # # # #
    A(x) = [x, spl_v.(x)]
    hull = sort(convex_hull(A.(Sw))) # Returns sorted convex hull
    hull_2 = (convex_hull(A.(Sw))) # Returns convex hull
    low = Vector{Bool}(undef, length(Sw))

    v_hull_points = [x[2] for x in hull] # Returns hull points
    Sw_hull_points = [x[1] for x in hull] # Returns hull points

    spl_v_hull = BSplineKit.interpolate(Sw_hull_points, v_hull_points, BSplineOrder(2))
    v_d_spl_hull = diff(spl_v_hull, Derivative(1)) #Interpolated derivative
    v_conv = spl_v_hull.(Sw)
    v_d_conv = v_d_spl_hull.(Sw)
    non_convex_region = (v_d_spl_hull.(Sw) .== v_d)

    # Interpolate and take derivatives
    spl_v_conv = BSplineKit.interpolate(Sw, v_conv, BSplineOrder(8), BSplineKit.Natural())
    v_d_spl_conv = diff(spl_v_conv, Derivative(1))
    vw_hat_conv = v_conv .+ Sn .* v_d_spl_conv.(Sw) # v̂_w from convex hull
    vn_hat_conv = spl_v_hull.(Sw) .- Sw .* v_d_spl_conv.(Sw) # v̂_n from convex hull
    return nothing 
end


function computeall(MTP::Vector{String},μ_w::Float64, Asw::Float64, Asn::Float64, data::Data; relperm = false)

    (; Sw, Sn, v, vw, vn, Ca) = data
    data = Data(Sw, Sn, v,vw, vn, Ca)
    M,T,P = parse.(Float64, [MTP[1], MTP[2], MTP[3]])
    μ_n = M*μ_w
    μ_e = @. (Sw + Sn*M)*μ_w
    Ca = sum((μ_e.*v.*10^-3)./(T.*0.01))./length(Sw)
    physicalparameters = PhysicalQuantities(μ_w,μ_n,μ_e,Ca)
    (; v_d, vw_hat, vn_hat, vm) = interpolations_derived(data)

    # Linear fit of v_m(v_d)
    llr, b_v = linear_region(v_d, vm, tol=0.1) # llr is largest linear region. 0.05
    a_v, _ = linreg(v_d[llr], vm[llr]) # Find a from llr
    vm_linear_fit = @. a_v + b_v * v_d
    # Homographies
    τ = -vn./(vw)
    Sn_tilde = Sn ./ (τ)
    Sw_points = [[Sw[i], Sn[i]] for i in eachindex(Sw)]  # Convert graphs to points
    Sw_curve_points = [[Sw[i], Sn_tilde[i]] for i in eachindex(Sw)] # Work with parameters. NOTE that we have dropped a minus sign here, to be added back later.

    # Sw, Sn_tilde parameterize projective line. Translate, then x/(x+y)
    # This is an alternative to computing via the cross ratio in 𝐏^2.
    # Manual check for smallest pressure gradients.

    if P == "0.20" 
        Asw = 0.0 ; Asn = 0.08
    else
        Asw = 0.0 ; Asn = 0.0
    end

    Sw_projected = @. (Sw .- Asw)./(Sw .+ Sn_tilde .- Asn .- Asw ) # x/(x+y)
    Sw_line_points = [[(Sw_projected)[i], (1 .- Sw_projected)[i]] for i in eachindex(Sw)] # Projection in affine coordinates. The cross-ratio is preserved under the above projection.
    # Convert to matrices. Use copy to get a matrix-type and not transpose.
    Sw_line_points_mat = copy(([first.(Sw_line_points) last.(Sw_line_points)])')
    Sw_curve_points_mat = copy(([first.(Sw_curve_points) last.(Sw_curve_points)])')
    Sw_points_mat = copy(([first.(Sw_points) last.(Sw_points)])')

    # Compute homography.
    hom_line = homography1d(Sw_points_mat, Sw_line_points_mat)

    # Extract parameters of homography.
    a, c, b, d = vec(hom_line) # ∞ mapped to a/c, -d/c mapped to ∞
    hom_line_c0 = copy(hom_line)
    hom_line_c0[2,1] = 0
    ∞, to_∞ = to_Sw.([ -d/c, a/c]) # ∞, to_∞ = -d/c, a/c

# # # # # # # # # #  Fixed points of homography # # # # # # # # # #
# Note that these do not correspond exactly to where the saturations cross due to the homography estimation
    x1, x2 = quadratic(c, (d - a), -b)
    x1_Sw, x2_Sw = foo(x1, x2) # Fixed points, Sw
    # println(∞, ", ", to_∞)

    # Define saturations from homography
    Sw_tilde = to_Sw.(proj_dehom([Sw, Sn], hom_line))# Corresponds to a hom. in the parameter Sw/Sn, giving \tilde{Sw/Sn}
    Sw_tilde_c0 = to_Sw.(proj_dehom([Sw, Sn], hom_line_c0))
    Sn_tilde = 1 .- Sw_tilde # NOTE: redefines Sn_tilde by the approximated equivalent.
    Sn_tilde_c0 = 1 .- Sw_tilde_c0


# # # # # # # # # #  Compute vm from the above saturations. # # # # # # # # # #
    k = 1 .*to_ratio.(Sw_tilde) 
    k_c0 = to_ratio.(Sw_tilde_c0) 
    vm_Sw_f_line = ((k).*vw_hat.*Sn .- vn_hat.*Sw)./((k).*Sn.^2 .+ Sw.^2)
    vm_Sw_f_line_c0 = ((k_c0).*vw_hat.*Sn .- vn_hat.*Sw)./((k_c0).*Sn.^2 .+ Sw.^2)

    # Recreated seepage velocities from vm as defined above
    vw_rec_hom, vn_rec_hom = (vw_hat .- Sn .* vm_Sw_f_line, vn_hat .+ Sw .* vm_Sw_f_line)
    vw_rec_hom_c0, vn_rec_hom_c0 = (vw_hat .- Sn .* vm_Sw_f_line_c0, vn_hat .+ Sw .* vm_Sw_f_line_c0)

    # # # # # # # # # # Linear fit of βv. Should give the same as vm_linear_fit # # # # # # # # # #
    βv = @. (k*Sn - Sw)/(k*Sn^2 + Sw^2)*v
    a_c, b_c = linreg(v_d[llr], βv[llr]) # Use same region as for vm_linear_fit
    vm_hom_linear_fit = @. a_c + (1 + b_c) * v_d

    return Results(vm_linear_fit, vm_Sw_f_line, vm_Sw_f_line_c0,vw_rec_hom, vn_rec_hom, Sw_tilde, Sn_tilde, k, a/d, b/d, c/d), data, interpolations_derived(data)
   
end

"""
Overloaded version, for relperm calculations.
"""
function computeall(data::Data, physicalquantities::RelpermQuantities, Asw::Float64, Asn::Float64, krw::Vector{Float64}, krn::Vector{Float64})
    (; Sw, Sn, v, vw, vn, Ca) = data
    (; σ, μ_w, μ_n, K, ϕ, ΔP) = physicalquantities
    M,T,P = μ_n/μ_w, σ, ΔP
    μ_e = @. (Sw + Sn*M)*μ_w
    (; v_d, vw_hat, vn_hat, vm) = interpolations_derived(data)

    # Compute vm from relperm curves
    # krw_Sw_d = Dierckx.derivative(Spline1D(Sw, krw./Sw; bc="extrapolate"), Sw)
    # krn_Sw_d = Dierckx.derivative(Spline1D(Sw, krn./Sn; bc="extrapolate"), Sw)
    # v0 = - K*ΔP/(μ_w*ϕ)
    # vm = @. μ_w * v0 * ((Sw / μ_w) * krw_Sw_d + (Sn / μ_n) * krn_Sw_d)

    # Linear fit of v_m(v_d)
    llr, b_v = linear_region(v_d, vm, tol=0.1) # llr is largest linear region. 0.05
    a_v, _ = linreg(v_d[llr], vm[llr]) # Find a from llr
    vm_linear_fit = @. a_v + b_v * v_d
    # Homographies
    τ = -vn./(vw)
    # τ = -krn .* μ_w .* Sw ./ (krw .* μ_n .* Sn)
    Sn_tilde = Sn ./ (τ)
    Sw_points = [[Sw[i], Sn[i]] for i in eachindex(Sw)]  # Convert graphs to points
    Sw_curve_points = [[Sw[i], Sn_tilde[i]] for i in eachindex(Sw)] # Work with parameters. NOTE that we have dropped a minus sign here, to be added back later.
    Sw_projected = @. (Sw .- Asw)./(Sw .+ Sn_tilde .- Asn .- Asw ) # x/(x+y)
    Sw_line_points = [[(Sw_projected)[i], (1 .- Sw_projected)[i]] for i in eachindex(Sw)] # Projection in affine coordinates. The cross-ratio is preserved under the above projection.
    # Convert to matrices. Use copy to get a matrix-type and not transpose.
    Sw_line_points_mat = copy(([first.(Sw_line_points) last.(Sw_line_points)])')
    Sw_curve_points_mat = copy(([first.(Sw_curve_points) last.(Sw_curve_points)])')
    Sw_points_mat = copy(([first.(Sw_points) last.(Sw_points)])')

    # Compute homography.
    hom_line = homography1d(Sw_points_mat, Sw_line_points_mat)
    # Extract parameters of homography.
    a, c, b, d = vec(hom_line) # ∞ mapped to a/c, -d/c mapped to ∞
    hom_line_c0 = copy(hom_line)
    hom_line_c0[2,1] = 0
    ∞, to_∞ = to_Sw.([ -d/c, a/c]) # ∞, to_∞ = -d/c, a/c

    # # # # # # # # # #  Fixed points of homography # # # # # # # # # #
    # Note that these do not correspond exactly to where the saturations cross due to the homography estimation
    x1, x2 = quadratic(c, (d - a), -b)
    x1_Sw, x2_Sw = foo(x1, x2) # Fixed points, Sw
    println(∞, ", ", to_∞)

    # Define saturations from homography
    Sw_tilde = to_Sw.(proj_dehom([Sw, Sn], hom_line))# Corresponds to a hom. in the parameter Sw/Sn, giving \tilde{Sw/Sn}
    Sw_tilde_c0 = to_Sw.(proj_dehom([Sw, Sn], hom_line_c0))
    Sn_tilde = 1 .- Sw_tilde # NOTE: redefines Sn_tilde by the approximated equivalent.
    Sn_tilde_c0 = 1 .- Sw_tilde_c0


    # # # # # # # # # #  Compute vm from the above saturations. # # # # # # # # # #
    k = to_ratio.(Sw_tilde) 
    k_c0 = to_ratio.(Sw_tilde_c0) 
    vm_Sw_f_line = ((k).*vw_hat.*Sn .- vn_hat.*Sw)./((k).*Sn.^2 .+ Sw.^2)
    vm_Sw_f_line_c0 = ((k_c0).*vw_hat.*Sn .- vn_hat.*Sw)./((k_c0).*Sn.^2 .+ Sw.^2)
    # Recreated seepage velocities from vm as defined above
    vw_rec_hom, vn_rec_hom = (vw_hat .- Sn .* vm_Sw_f_line, vn_hat .+ Sw .* vm_Sw_f_line)
    # Test: recreate relperms instead
    vw_rec_hom, vn_rec_hom = (vw_hat .- Sn .* vm_Sw_f_line, vn_hat .+ Sw .* vm_Sw_f_line)
    vw_rec_hom_c0, vn_rec_hom_c0 = (vw_hat .- Sn .* vm_Sw_f_line_c0, vn_hat .+ Sw .* vm_Sw_f_line_c0)

    # # # # # # # # # # Linear fit of βv. Should give the same as vm_linear_fit # # # # # # # # # #
    βv = @. (k*Sn - Sw)/(k*Sn^2 + Sw^2)*v
    a_c, b_c = linreg(v_d[llr], βv[llr]) # Use same region as for vm_linear_fit
    vm_hom_linear_fit = @. a_c + (1 + b_c) * v_d

    return Results(vm_linear_fit, vm_Sw_f_line, vm_Sw_f_line_c0,vw_rec_hom, vn_rec_hom, Sw_tilde, Sn_tilde, k, a/d, b/d, c/d), data, interpolations_derived(data),  vw_rec_hom_c0, vn_rec_hom_c0
   
end

"""
Compute the non-Euclidean geometric quantities
For theory verification.
"""
function noneuclidean()
    l_∞ = [1, 1, -1]
    B_init = [0 1 1; 1 0 1; 1 1 2] # Dual conic, initial
    B_t = [1 0 0 ; 0 -1 0 ; 0 0 0] # Dual conic, new
    eigsB_init, Q =eigen(B_init)
    Qt = Q*diagm([1,1,sqrt(3)]) # Define new orthogonal matrix by normalizing eigenvalues
    A1 = [1 1 -1; 1 1 -1; -1 -1 1] # Primal conic
    ll = inv(Qt)' * A1 * inv(Qt) # Compute transformed primal from transformed dual
    # Decompose ll
    y = ll[:,1] # Random nonzero column. Scalar times y is the sought after line
    μ = sqrt(((y*y')[1,1]/ll[1,1])^-1) # Compute scaling parameter
    lt = μ .* y # Transformed l_∞
    l_∞_lt_inter = cross(l_∞, lt) # Intersection. (0,1,1) is a fixed point.
    return nothing
end
