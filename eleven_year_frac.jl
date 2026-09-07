"
Look at 11 years of ERA5 data
 and see if there is any trend
 in the fractal dimension of the data over time.
"

include("func.jl")
include("setup.jl")

ds = NCDataset(joinpath(@__DIR__, "pl201001.nc"))
q = ds["q"][:, 81:241, :, :]   # (1440, 721, 3, 124)
t = ds["t"][:, 81:241, :, :]   # same shape

lon = ds["longitude"][:]
lat = ds["latitude"][:]
level = ds["level"][:]
time = ds["time"][:]

#rhi_200 = rhi_calc(q[:, :, 1, :], t[:, :, 1, :], 20000)
#rhi_250 = rhi_calc(q[:, :, 2, :], t[:, :, 2, :], 25000)
rhi_300 = rhi_calc(q[:, :, 3, :], t[:, :, 3, :], 30000);

anim = @animate for i in 1:124
    heatmap(rhi_300[:, :, i]';
               xlabel = "Longitude",
               ylabel = "Latitude",
               title = "RHI, Zeitschritt $i",
               colorbar_title = "RHI",
               clims = (-0.2, 1.7))
end

gif(anim, fps = 10)

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