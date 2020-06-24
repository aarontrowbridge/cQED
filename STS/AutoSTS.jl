module AutoSTS

using DelimitedFiles, Plots, LaTeXStrings

export Qubit, STS
export get_slice, plot_slice, extract_resonance_freqs

mutable struct Qubit
    f_c::Float64
    g::Float64
    Π::Float64
    I_ss::Float64
    f_ge_max::Float64
    d::Float64
end

mutable struct STS
    path::String
    data::Matrix{Float64}
    freq_lims::Tuple{Float64,Float64}
    volt_lims::Tuple{Float64,Float64}
    Δf::Float64
    ΔV::Float64
    f_range::StepRangeLen
    V_range::StepRangeLen

    function STS(path::String)
        data = readdlm(path, '\t', Float64)
        f = match(r"fr(.+?)_(.+?)_", path)
        V = match(r"_V(.+?)_(.+?)_", path)
        fmin, fmax = parse(Float64, f[1]), parse(Float64, f[2])
        Vmin, Vmax = parse(Float64, V[1]), parse(Float64, V[2])
        Δf = (fmax - fmin) / (size(data, 1) - 1)
        ΔV = (Vmax - Vmin) / (size(data, 2) - 1)
        f_range = fmin:Δf:fmax
        V_range = Vmin:ΔV:Vmax
        new(path, data, (fmin, fmax), (Vmin, Vmax), Δf, ΔV, f_range, V_range)
    end
end

Plots.heatmap(sts::STS) = heatmap(sts.V_range, sts.f_range, sts.data)

get_slice(sts::STS, V::Float64) =
    sts.data[:, floor(Int, (V - sts.volt_lims[1]) / sts.ΔV)]

plot_slice(sts::STS, V::Float64) =
    plot(sts.f_range, get_slice(sts, V),
         xlabel=L"f_{p}[GHz]]",
         ylabel=L"|S_{21}|")

function extract_resonance_freqs(sts::STS)
    f_r = []
    for col in eachcol(sts.data)
        i = argmin(col)
        f = sts.freq_lims[1] + (i - 1) * sts.Δf
        push!(f_r, f)
    end
    f_r
end













end
