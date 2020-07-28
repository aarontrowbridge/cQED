module ResonanceFit

using CircleFit
using LsqFit

export S21Fit

struct S21Fit{T<:AbstractFloat}
    fᵣ::T
    Qi::T
    Ql::T
    Qc::Complex{T}
    θ₀::T

    function S21Fit{T}(fs::Vector{T},
                       zs::Vector{Complex{T}}) where {T<:AbstractFloat}
        circle = FitCircle{T}(real.(zs), imag.(zs))
        ϕ₀ = -asin(circle.yc / circle.r₀)
        zc = complex(circle.xc, circle.yc)
        zs .-= zc
        ϕ₀ = angle(zc)
        θ₀ᵢ = mod((ϕ₀ + π), π)
        fᵣᵢ = fs[argmin(abs.(zs))]
        θ₀, Ql, fᵣ = phase_angle_fit(fs, angle.(zs), θ₀ᵢ, fᵣᵢ)
        Qc = Ql / (2 * circle.r₀) * exp(im * ϕ₀)
        Qi = 1 / (1 / Ql - real(1 / Qc))
        new(fᵣ, Qi, Ql, Qc, θ₀)
    end
end

function phase_angle_fit(fs::Vector{T},
                         θs::Vector{T},
                         θ₀ᵢ::T,
                         fᵣᵢ::T) where {T<:AbstractFloat}
    θ(f, p) = p[1] .+ 2atan.(2 * p[2] * (1 .- f ./ p[3]))
    p₀ = [θ₀ᵢ, 1e4, fᵣᵢ]
    fit = curve_fit(θ, fs, θs, p₀)
    return fit.param
end

end
