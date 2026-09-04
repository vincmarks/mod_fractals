
#  Phase 2: fBm-Felder synthetisieren und die Mess-Pipeline gegenprüfen
#
#  Idee: ein 2D-Feld mit vorgegebenem Spektralexponenten β erzeugen, β über
#  das radial gemittelte Leistungsspektrum zurückmessen, und schauen ob das
#  wieder rauskommt was reingesteckt wurde.
#
#      S(k) ~ k^(-β)

# loading packages
using Pkg
Pkg.activate(@__DIR__)
# Pkg.add(["FFTW", "Plots", "Random", "Statistics", "LinearAlgebra", "WAV", "NCDatasets", "LaTeXStrings"])
Pkg.instantiate()
using FFTW
using Random
using Statistics
using LinearAlgebra
using Plots
using Plots.PlotMeasures      # mm, cm, px, pt für die Ränder
using Printf
using WAV
using NCDatasets
using LaTeXStrings

#######
# fbm erzeugen und grid

"""
    kgrid(N) -> Matrix{Float64}

Betrag des Wellenzahlvektors für ein N x N Feld
`fftfreq(N, N)` liefert 0, 1, …, N÷2-1, -N÷2, …, -1.
"""
function kgrid(N::Integer)
    kx = collect(fftfreq(N, N))          # Spaltenrichtung (Dimension 1)
    return sqrt.(kx .^ 2 .+ (kx') .^ 2)  # kx' = Zeilenrichtung (Dimension 2)
end

kgrid1(n::Integer) = abs.(collect(fftfreq(n, n)))

"""
    synth_fbm(N, β; rng) -> Matrix{Float64}

2D-Feld mit Leistungsspektrum ~ k^(-β).

Weißes Gauß-Rauschen rein, FFT, mit k^(-β/2) filtern, zurück-transformieren.
Weil das Ausgangsrauschen reell ist, ist seine FFT automatisch hermitesch --
das Ergebnis bleibt also bis auf Rundungsfehler reell, ganz ohne dass man
sich selbst um die Symmetrie kümmern müsste.
"""
function synth_fbm(N::Integer, β::Real; rng::AbstractRNG = Random.default_rng())
    white = randn(rng, N, N)
    F = fft(white)

    k = kgrid(N)
    filt = zeros(Float64, N, N)
    nz = k .> 0                          # k = 0 (DC) bleibt null
    filt[nz] .= k[nz] .^ (-β / 2)        # /2, weil |F|² den Exponenten verdoppelt!

    return real.(ifft(F .* filt))
end

"""1D-Variante, z.B. für Klangbeispiele."""
function synth_fbm1d(n::Integer, β::Real; rng::AbstractRNG = Random.default_rng())
    white = randn(rng, n)
    F = fft(white)

    k = kgrid1(n)
    filt = zeros(Float64, n)
    nz = k .> 0
    filt[nz] .= k[nz] .^ (-β / 2)

    return real.(ifft(F .* filt))
end

####
# hann fenster - aktuell noch nicht zu benutzten

"""Hann-Fenster, identisch zu `numpy.hanning` (Endpunkte exakt null)."""
hann(n::Integer) = 0.5 .* (1 .- cos.(2π .* (0:n-1) ./ (n - 1)))

"""
    radial_spectrum(field; window=false, combine=:mean) -> (k, S)

Radial gemitteltes Leistungsspektrum eines quadratischen Feldes.

`window`  -- 2D-Hann-Fenster gegen spektrale Leckage. Standardmäßig aus.
`combine` -- `:mean` ist die radiale Mittelung
"""
function radial_spectrum(field::AbstractMatrix; window::Bool = false,
                         combine::Symbol = :mean)
    f = field .- mean(field)

    if window
        w = hann(size(f, 1)) * hann(size(f, 2))'   # äußeres Produkt -> Matrix
        f = f .* w ./ sqrt(mean(w .^ 2))           # Leistungsnormierung
    end

    P = abs2.(fft(f)) ./ length(f)

    N = size(field, 1)
    kint = round.(Int, kgrid(N))
    kmax = N ÷ 2

    total = zeros(Float64, kmax)
    cnt   = zeros(Int, kmax)
    @inbounds for i in eachindex(kint)
        κ = kint[i]
        if 1 ≤ κ ≤ kmax
            total[κ] += P[i]
            cnt[κ]   += 1
        end
    end

    S = combine === :mean ? total ./ max.(cnt, 1) : total
    return collect(1:kmax), S
end

#######
# fitting

"""
    fit_beta(k, S; kmin, kmax, nbins=25) -> (β, achsenabschnitt, R²)

Steigung im Log-Log-Plot über logarithmisch gleichverteilte Bins.

Das Log-Binning ist nötig, weil sonst rund 90 % der Punkte im hochfrequenten
Bereich liegen und die Regression komplett dominieren würden.
"""
function fit_beta(k::AbstractVector, S::AbstractVector;
                  kmin::Real = 4, kmax::Real = length(S) ÷ 2, nbins::Integer = 25)

    m = (k .>= kmin) .& (k .<= kmax) .& (S .> 0)
    lk = log10.(k[m])
    ls = log10.(S[m])

    # logarithmisch gleich breite Bins; nextfloat, damit der Maximalwert
    # noch in den letzten Bin fällt
    edges = range(minimum(lk), nextfloat(maximum(lk)); length = nbins + 1)

    xs = Float64[]
    ys = Float64[]
    for b in 1:nbins
        sel = (lk .>= edges[b]) .& (lk .< edges[b+1])
        if any(sel)
            push!(xs, mean(lk[sel]))
            push!(ys, mean(ls[sel]))
        end
    end

    X = [ones(length(xs)) xs]
    icept, slope = X \ ys

    resid = ys .- X * [icept, slope]
    r2 = 1 - sum(abs2, resid) / sum(abs2, ys .- mean(ys)) # bestimmhteismaß

    return -slope, icept, r2
end

"""
    measure(field; kmin=4, kmax=nothing, kwargs...) -> (β, achsenabschnitt, R², k, S)

Komplette Messkette an einem Feld: Spektrum + Fit in einem Aufruf. Ohne
`kmax` input wird bei N/4 abgeschnitten, um den verrauschten Bereich
nahe der Nyquist-Frequenz gar nicht erst mitzufitten.
"""
function measure(field::AbstractMatrix; kmin::Real = 4, kmax = nothing, kwargs...)
    k, S = radial_spectrum(field; kwargs...)
    kmx = kmax === nothing ? maximum(k) ÷ 2 : kmax
    β, icept, r2 = fit_beta(k, S; kmin = kmin, kmax = kmx)
    return β, icept, r2, k, S
end
