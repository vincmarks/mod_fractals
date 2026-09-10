"
Look at 11 years of ERA5 data
 and see if there is any trend
 in the fractal dimension of the data over time.
"

## Imports

include("func.jl")
include("setup.jl")

## Read in the data

ds = NCDataset(joinpath(@__DIR__, "pl201001.nc"))
q = ds["q"][461:670, 71:280, :, :]   # (211, 211, 3, 124) 72.5 N - 20 N; -65 E - -12.5 E
t = ds["t"][461:670, 71:280, :, :]   # same shape

lon = ds["longitude"][:]
lat = ds["latitude"][:]
level = ds["level"][:]
time = ds["time"][:]

## Calculate the RHI

#rhi_200 = rhi_calc(q[:, :, 1, :], t[:, :, 1, :], 20000)
#rhi_250 = rhi_calc(q[:, :, 2, :], t[:, :, 2, :], 25000)
rhi_300 = rhi_calc(q[:, :, 3, :], t[:, :, 3, :], 30000);

## Convert into Fourier space

# Ein 2D-Fourierspektrum für einen ausgewählten Zeitpunkt berechnen.

n_times = size(rhi_300, 3)

spectra = Matrix{Float64}(undef, 105, n_times)
betas = Vector{Float64}(undef, n_times)
r2_values = Vector{Float64}(undef, n_times)



for t_idx in 1:n_times
    field = rhi_300[:, :, t_idx]
    β, _, r2, k, S = measure(field; window = true)
    spectra[:, t_idx] = S
    betas[t_idx] = β
    r2_values[t_idx] = r2
end

S_month = vec(mean(spectra; dims = 2))
β_month, icept_month, r2_month = fit_beta(k, S_month)

## Plot the results

p = plot(k, S_month; seriestype = :scatter, xscale = :log10, yscale = :log10,
         markersize = 2, markerstrokewidth = 0, label = "S(k)");
kf = [4.0, 210 / 4]
plot!(p, kf, (10 ^ icept_month) .* kf .^ (-β_month); linestyle = :dash,
      label = @sprintf("Fit: β = %.3f  (R² = %.4f)", β_month, r2_month));
plot!(p; xlabel = "Wellenzahl k", ylabel = "Leistung S(k)",
      title = "rhi_300, Monatsmittel aus $n_times Zeitpunkten")


## Create the animation

lat_sel = lat[81:241]

anim = @animate for i in 1:124
    field = reverse(rhi_300[:, :, i]', dims = 1)

    heatmap(lon, reverse(lat_sel), field;
            xlabel = "Longitude",
            ylabel = "Latitude",
            title = "RHI, Zeitschritt $i",
            colorbar_title = "RHI",
            clims = (-0.2, 1.7))
end

gif(anim, fps = 10)