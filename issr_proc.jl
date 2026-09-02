# Code übersetzt aus Python
include("func.jl")

ds = NCDataset(joinpath(@__DIR__, "pl201107.nc"))
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
# Ein Zeitpunkt aus rhi_200 durch die fBm-Pipeline aus setup.jl

include("setup.jl")

t_idx = 50
slice = rhi_200[:, :, t_idx]              # (lon, lat) = (1440, 721)
slice = rhi_300[:, :, t_idx] 
# radial_spectrum/kgrid gehen von einem quadratischen Feld aus, also einen
# N×N-Ausschnitt rausschneiden. Für halbwegs isotrope Pixelabstände (1°
# Länge entspricht am Äquator in km ungefähr 1° Breite, an den Polen aber
# viel weniger) den Ausschnitt in Breitenrichtung um 0° zentrieren.
N_real = 512
i0 = (size(slice, 1) - N_real) ÷ 2 + 1
j0 = clamp(argmin(abs.(lat .- 0)) - N_real ÷ 2, 1, length(lat) - N_real + 1)
field = slice[i0:i0+N_real-1, j0:j0+N_real-1]

# echter, nicht-periodischer Ausschnitt -> diesmal Fenster einschalten
β, icept, r2, k, S = measure(field; window = true)
@printf("rhi_200, Zeitschritt %d: β = %.3f  (R² = %.4f)\n", t_idx, β, r2)

p = plot(k, S; seriestype = :scatter, xscale = :log10, yscale = :log10,
         markersize = 2, markerstrokewidth = 0, label = "S(k)");
kf = [4.0, N_real / 4]
plot!(p, kf, (10 ^ icept) .* kf .^ (-β); linestyle = :dash,
      label = @sprintf("Fit: β = %.3f  (R² = %.4f)", β, r2));
plot!(p; xlabel = "Wellenzahl k", ylabel = "Leistung S(k)",
      title = "rhi_300, Zeitschritt $t_idx")
isdir("out") || mkdir("out")
savefig(p, joinpath("out", "05_rhi300_spektrum.png"))
