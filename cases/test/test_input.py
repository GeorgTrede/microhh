"""
Simplified Jaenschwalde setup.
Profiles are slightly idealised, from 04:00 UTC ERA5 data.
Wind is rotated to be perfectly westerly.
"""

import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

pl.close("all")
pl.ion()

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open("test.ini") as f:
    for line in f:
        if line.split("=")[0] == "ktot":
            kmax = int(line.split("=")[1])
        if line.split("=")[0] == "zsize":
            zsize = float(line.split("=")[1])
        if line.split("=")[0] == "ysize":
            ysize = float(line.split("=")[1])

# Enable resolved plume rise:
sw_plume_rise = True

# Vertical grid LES
dz = zsize / kmax
z = np.arange(0.5 * dz, zsize, dz)

z_u = np.array([0, 270, 3000, 5000])
v_u = np.array([2.3, 8.5, 0.6, 5.7])*0.75
u = np.interp(z, z_u, v_u)

v = np.zeros(kmax)
co2 = np.zeros(kmax)

# # Surface fluxes, again idealised from ERA5.
# t0 = 4 * 3600
# t1 = 16 * 3600
# td1 = 12 * 3600
# td2 = 14 * 3600

# time = np.linspace(t0, t1, 32)
# wthl = 0.17   * np.sin(np.pi * (time-t0) / td1)
# wqt  = 8.3e-5 * np.sin(np.pi * (time-t0) / td2)

# Write input NetCDF file
nc_file = nc.Dataset(
    "test_input.nc", mode="w", datamodel="NETCDF4", clobber=True
)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init")
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))
nc_co2 = nc_group_init.createVariable("co2", float_type, ("z"))
nc_co2_inflow = nc_group_init.createVariable("co2_inflow", float_type, ("z"))

nc_z[:] = z[:]
nc_u[:] = u[:]
nc_v[:] = v[:]
nc_co2[:] = co2[:]
nc_co2_inflow[:] = co2[:]

# nc_group_tdep = nc_file.createGroup("timedep")
# nc_group_tdep.createDimension("time_surface", time.size)
# nc_time_surface = nc_group_tdep.createVariable(
#     "time_surface", float_type, ("time_surface")
# )
# nc_thl_sbot = nc_group_tdep.createVariable("thl_sbot", float_type, ("time_surface"))
# nc_qt_sbot = nc_group_tdep.createVariable("qt_sbot" , float_type, ("time_surface"))

# nc_time_surface[:] = time
# nc_thl_sbot[:] = wthl
# nc_qt_sbot[:] = wqt

nc_file.close()