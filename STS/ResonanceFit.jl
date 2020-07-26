module ResonanceFit

using CircleFit
using JuMP, Ipopt, GLPK

export S21Fit

struct S21Fit{T<:AbstractFloat}
    fᵣ::T
    Qi::T
    Ql::T
    Qc::Complex{T}

    function S21Fit{T}(fs::Vector{T}, zs::Vector{Complex{T}}) where {T<:AbstractFloat}
        circle = FitCircle{T}(real.(zs), imag.(zs))
        ϕ₀ = -asin(circle.yc / circle.r₀)
        zc = complex(circle.xc, circle.yc)
        zs .-= zc
        ϕ₀ = angle(zc)
        θᵢ = mod((ϕ₀ + π), π)
        fᵣ, Ql = phase_angle_fit(fs, zs, θᵢ)
        Qc = Ql / (2 * circle.r₀) * exp(im * ϕ₀)
        Qi = 1 / (1 / Ql - real(1 / Qc))
        new(fᵣ, Qi, Ql, Qc)
    end
end

function phase_angle_fit(fs::Vector{T}, zs::Vector{Complex{T}}, θᵢ::T) where {T<:AbstractFloat}
    xs, ys = real.(zs), imag.(zs)
    model = Model(Ipopt.Optimizer)
    @variable(model, fs[1] <= fᵣ <= fs[end])
    @variable(model, 850 <= Ql <= 1100)
    @variable(model, θ₀, start=θᵢ)
    @NLobjective(model, Min,
                 sum(abs2(atan(y / x) - θ₀ - 2atan(2Ql * (1 - f / fᵣ)))
                     for (f, x, y) in zip(fs, xs, ys)))
    optimize!(model)
    return value(fᵣ), value(Ql)
end

end
