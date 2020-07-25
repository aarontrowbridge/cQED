module CircleFit

using LinearAlgebra, NLsolve

export FitCircle

const B̂ = [0  0 0 -2;
           0  1 0  0;
           0  0 1  0;
           -2 0 0  0]

struct FitCircle{T<:AbstractFloat}
    xs::Vector{T}
    ys::Vector{T}
    xc::T
    yc::T
    r₀::T
    N::Int

    function FitCircle{T}(xs::Vector{T}, ys::Vector{T}) where {T<:AbstractFloat}
        N = length(xs)
        M̂ = M(xs, ys, N)
        f!(F, η) = setindex!(F, det(M̂ - η[1] * B̂), 1)
        η′ = nlsolve(f!, [0.], method = :newton, ftol = 1e-12).zero[1]
        A′ = nullspace(M̂ - η′ * B̂)[:,1]
        A, B, C, D = A′
        xc = -B / 2A
        yc = -C / 2A
        r₀ = 1 / (2 * abs(A))
        new(xs, ys, xc, yc, r₀, N)
    end
end

function M(x::Vector{T}, y::Vector{T}, N::Int) where {T<:AbstractFloat}
    z = [xᵢ^2 + yᵢ^2 for (xᵢ, yᵢ) in zip(x, y)]

    Mxx = dot(x, x); Myy = dot(y, y); Mzz = dot(z, z)
    Mxz = dot(x, z); Myz = dot(y, z); Mxy = dot(x, y)
    Mx = sum(x); My = sum(y); Mz = sum(z)

    [Mzz Mxz Myz Mz;
     Mxz Mxx Mxy Mx;
     Myz Mxy Myy My;
     Mz  Mx  My  N]
end

function SNR(circ::FitCircle)
    rs = [sqrt((x - circ.xc)^2 + (y - circ.yc)^2) for (x, y) in zip(circ.xs, circ.ys)]
    σᵣ = sqrt(1 / (circ.N - 1) * sum(rs .- circ.r₀))
    circ.r₀ / σᵣ
end



end
