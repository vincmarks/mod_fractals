# Code übersetzt aus Python

using Pkg
Pkg.add("NCDatasets")

using NCDatasets

include("func.jl")

ds = NCDataset("/Users/helenaschuh/Projects/Master/Mod_Prak/mod_fractals/pl201001.nc")



q = ds["q"][:, :, :, :]   # (1440, 721, 3, 124)
t = ds["t"][:, :, :, :]   # same shape

lon = ds["longitude"][:]
lat = ds["latitude"][:]
level = ds["level"][:]
time = ds["time"][:]

close(ds)

rhi_200 = rhi_calc(q[:, :, 1, :], t[:, :, 1, :], 20000)
