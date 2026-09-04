# Code übersetzt aus Python
include("func.jl")
include("setup.jl")

ds = NCDataset(joinpath(@__DIR__, "pl201001.nc"))
q = ds["q"][:, :, :, :]   # (1440, 721, 3, 124)
t = ds["t"][:, :, :, :]   # same shape

lon = ds["longitude"][:]
lat = ds["latitude"][:]
level = ds["level"][:]
time = ds["time"][:]

#close(ds)

rhi_200 = rhi_calc(q[:, :, 1, :], t[:, :, 1, :], 20000)
rhi_250 = rhi_calc(q[:, :, 2, :], t[:, :, 2, :], 25000)
rhi_300 = rhi_calc(q[:, :, 3, :], t[:, :, 3, :], 30000)


########
# Alle Zeitpunkte aus rhi_300 durch die fBm-Pipeline aus setup.jl

# radial_spectrum/kgrid gehen von einem quadratischen Feld aus, also einen
# N×N-Ausschnitt rausschneiden. Für halbwegs isotrope Pixelabstände (1°
# Länge entspricht am Äquator in km ungefähr 1° Breite, an den Polen aber
# viel weniger) den Ausschnitt in Breitenrichtung um 0° zentrieren.
N_real = 512
j0 = clamp(argmin(abs.(lat .- 0)) - N_real ÷ 2, 1, length(lat) - N_real + 1)
i0 = (size(rhi_300, 1) - N_real) ÷ 2 + 1
n_times = size(rhi_300, 3)
k = collect(1:(N_real ÷ 2))

# Ein Spektrum pro Zeitpunkt berechnen und anschließend über den Monat mitteln.
spectra = Matrix{Float64}(undef, N_real ÷ 2, n_times)
betas = Vector{Float64}(undef, n_times)
r2_values = Vector{Float64}(undef, n_times)

for t_idx in 1:n_times
    slice = rhi_300[:, :, t_idx]
    field = slice[i0:i0+N_real-1, j0:j0+N_real-1]

    # echter, nicht-periodischer Ausschnitt -> Fenster einschalten
    β, _, r2, k, S = measure(field; window = true)
    spectra[:, t_idx] = S
    betas[t_idx] = β
    r2_values[t_idx] = r2
end

S_month = vec(mean(spectra; dims = 2))
β_month, icept_month, r2_month = fit_beta(k, S_month)
@printf("rhi_300, %d Zeitpunkte: β = %.3f  (R² = %.4f)\n",
        n_times, β_month, r2_month)
@printf("Zeitpunkt-β: Mittelwert = %.3f, Standardabweichung = %.3f\n",
        mean(betas), std(betas))

p = plot(k, S_month; seriestype = :scatter, xscale = :log10, yscale = :log10,
         markersize = 2, markerstrokewidth = 0, label = "S(k)");
kf = [4.0, N_real / 4]
plot!(p, kf, (10 ^ icept_month) .* kf .^ (-β_month); linestyle = :dash,
      label = @sprintf("Fit: β = %.3f  (R² = %.4f)", β_month, r2_month));
plot!(p; xlabel = "Wellenzahl k", ylabel = "Leistung S(k)",
      title = "rhi_300, Monatsmittel aus $n_times Zeitpunkten")
isdir("out") || mkdir("out")
savefig(p, joinpath("out", "05_rhi300_monatsspektrum_jan10.png"))
