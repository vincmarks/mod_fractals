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
q = ds["q"][461:671, 71:281, :, :]   # (211, 211, 3, 124) 72.5 N - 20 N; -65 E - -12.5 E
t = ds["t"][461:671, 71:281, :, :]   # same shape

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