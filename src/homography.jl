# The following code is obtained from the package ImageProjectiveGeometry.jl,
# found here:  https://github.com/peterkovesi/ImageProjectiveGeometry.jl?tab=readme-ov-file
# The code is include explicitly due to compatibility issues with Interpolations.jl.
# License:
#= 
    Copyright (c) 2016: Peter Kovesi.

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

"Normalization"
function hnormalise!(x::Union{Array{T,2},Vector{T}}) where T <: Real
    (rows, npts) = (size(x,1), size(x,2))  # This handles 2D arrays and vectors

    for n = 1:npts
        if abs(x[rows, n]) > eps(T)   # point not at infinity
            for r = 1:rows-1
                x[r, n] /= x[rows, n]
            end
            x[rows, n] = 1
        end
    end

    return x
end

"Normalize 1d homogeneous coordinates"
function normalise1dpts(ptsa::Array{T1,2}) where T1 <: Real
    # ? Can this be unified with normalise2dpts ?

    pts = copy(ptsa)

    if size(pts,1) != 2
        error("pts must be 2xN")
    end

    if any(pts[2,:] .== 0)
        @warn("Attempt to normalise a point at infinity")
        return pts, I(2)
    end

    # Ensure homogeneous coords have scale of 1
    hnormalise!(pts)

    # Get centroid and mean distance from centroid
    c = mean(view(pts, 1,:))
    meandist = mean(abs.(pts[1,:] .- c))
    scale = 1/meandist

    T = [scale    -scale*c
         0         1      ]

    pts = T*pts

    return pts, T
end

" Compute 1d homography"
function homography1d(x1::Array{T1,2}, x2::Array{T2,2}) where {T1 <: Real, T2 <: Real}

    if size(x1) != size(x2)
        error("x1 and x2 must have same dimensions")
    end

    (dim, Npts) = size(x1)

    # Attempt to normalise each set of points so that the origin
    # is at centroid and mean distance from origin is 1.
    (x1n, Tx1) = normalise1dpts(x1)
    (x2n, Tx2) = normalise1dpts(x2)

    # Note that it may have not been possible to normalise
    # the points if one was at infinity so the following does not
    # assume that scale parameter w = 1.
    A = zeros(2*Npts,4)

    for n = 1:Npts
        X = x1n[:,n]'
        x = x2n[1,n]
        w = x2n[2,n]
        A[n,:] = [-w*X x*X]
    end

    (U,S,V) = svd(A)

    # Extract homography
    H = reshape(V[:,4],2,2)'

    # Denormalise
    return H = Tx2\H*Tx1
end
