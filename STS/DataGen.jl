module DataGen

using Distributions

export S21Params, S₂₁, Qloaded

mutable struct S21Params{T<:AbstractFloat}
    a::T
    α::T
    ϕ::T
    τ::T
    fᵣ::T
    Qc::T
    Qi::T

    S21Params{T}(; a=0.1, α=0.4π, ϕ=0.3π, τ=50, fᵣ=5, Qc=1e3, Qi=1e4) where {T<:AbstractFloat} =
        new(a, α, ϕ, τ, fᵣ, Qc, Qi)
end

function S₂₁(f::T, p::S21Params; SNR::T=100.) where {T<:AbstractFloat}
    Ql = 1 / (1 / p.Qi + real(1 / p.Qc))
    z = p.a * exp(im * (p.α - 2π * f * p.τ)) *
            (1 - (Ql / p.Qc) * exp(im * p.ϕ) / (1 + 2im * Ql * (f / p.fᵣ - 1)))
    z
end

Qloaded(Qi, Qc) = 1 / (1 / Qi + real(1 / Qc))


end
