include("setup.jl")
  # GDAL-Backend, das Rasters zum Lesen von GeoTIFFs braucht
Rasters.checkmem!(false)   # Sandbox meldet fälschlich zu wenig freien RAM

OUT = "out"
isdir(OUT) || mkdir(OUT)

r = Raster(joinpath(@__DIR__, "alps_cop30.tif"); checkmem = false)

p_overview = plot(r; title = "Alpen – Höhendaten", clims = extrema(skipmissing(r)))
savefig(p_overview, joinpath(OUT, "alps_00_overview.png"))

########
# topographie

raw = r.data
gap = ismissing.(raw)
@printf("Ausschnitt: %d x %d px (ganze Kachel), %.2f%% Lücken -> auf Mittelwert gesetzt\n",
        size(raw, 1), size(raw, 2), 100 * count(gap) / length(raw))

const field = Float64.(coalesce.(raw, mean(skipmissing(raw))))
N = size(field, 1)

p_field = heatmap(field; c = :terrain, title = "")

savefig(p_field, joinpath(OUT, "alps.png"))

########
# Powersprektum, look at two regimes 
# large scales
KMIN = 11
KMAX = 108 


β, icept, r2, k, S = measure(field; kmin = KMIN, kmax = KMAX, window = true)
@printf("Alpen-Höhenfeld: β = %.3f  (R² = %.4f)\n", β, r2)

p = plot(k, S; seriestype = :scatter, xscale = :log10, yscale = :log10,
         markersize = 2, markerstrokewidth = 0, label = L"S(k)");
         kf = [Float64(KMIN), Float64(KMAX)]
plot!(p, kf, (10 ^ icept) .* kf .^ (-β); linestyle = :dashdot, linewidth = 3,
      label = L"\beta" * @sprintf(" = %.3f  (R² = %.4f)", β, r2));

plot!(p; xlabel = L"k", ylabel = L"S(k)",
      title = "Powerspektrum Alpen")
savefig(p, joinpath(OUT, "alps_powerspekt_low.png"))

KMIN = 108
KMAX = 432

β, icept, r2, k, S = measure(field; kmin = KMIN, kmax = KMAX, window = true)
@printf("Alpen-Höhenfeld: β = %.3f  (R² = %.4f)\n", β, r2)

p = plot(k, S; seriestype = :scatter, xscale = :log10, yscale = :log10,
         markersize = 2, markerstrokewidth = 0, label = L"S(k)");
         kf = [Float64(KMIN), Float64(KMAX)]
plot!(p, kf, (10 ^ icept) .* kf .^ (-β); linestyle = :dashdot, linewidth = 3,
      label = L"\beta" * @sprintf(" = %.3f  (R² = %.4f)", β, r2));

plot!(p; xlabel = L"k", ylabel = L"S(k)",
      title = "Powerspektrum Alpen")
savefig(p, joinpath(OUT, "alps_powerspektrum_high.png"))
