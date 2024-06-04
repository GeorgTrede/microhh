import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc
from scipy.integrate import odeint

pl.close("all")
pl.ion()

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open("./test.ini") as f:
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
dz = zsize / kmax  # type: ignore
z = np.arange(0.5 * dz, zsize, dz)

v_thl = np.array([285.7, 291.9, 293, 297.4, 307])
z_thl = np.array([0, 400, 2000, 2500, 5000])
thl = np.interp(z, z_thl, v_thl)

z_qt = np.array([0, 400, 2000, 2500, 5000])
v_qt = np.array([6.2, 4.93, 3.61, 1, 0.3]) / 1000
qt = np.interp(z, z_qt, v_qt)

# cubic interpolation of u
# z_u = np.array([0, 500, 1500, 2000, 3000, 4000, 5000])
# v_u = np.array([2.3, 6.5, 4.4, 4.0, 4.0, 4.4, 5.7])
# f = interp1d(z_u, v_u, kind="cubic")
# u = f(z)

# linear interpolation of u
# z_u = np.array([0, 1000, 2500, 5000])
# v_u = np.array([2.3, 6.5, 4.0, 5.7])
# u = np.interp(z, z_u, v_u)

# logarithmic profile of u
# u_star = 0.3
# z0 = 0.1
# u = u_star / 0.4 * np.log((z + 1) / z0)


# Monin-Obukhov profile of u
def dUdz(U, z, u_star):
    L = -15  # Obukhov length
    # \frac{dU}{dz} = \frac{ustar}{0.4 z} (1 - 15 \frac{z}{L})^{-1/4}
    return u_star / 0.4 / z * (1 - 15 * z / L) ** (-1 / 4)


u_star = 0.5
z0 = 0.1
u = odeint(dUdz, 0, z, args=(u_star,)).flatten()

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
nc_file = nc.Dataset("test_input.nc", mode="w", datamodel="NETCDF4", clobber=True)  # type: ignore

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init")
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_u_nudge = nc_group_init.createVariable("u_nudge", float_type, ("z"))
nc_nudge_fac = nc_group_init.createVariable("nudgefac", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
# nc_qt = nc_group_init.createVariable("qt", float_type, ("z"))
nc_co2 = nc_group_init.createVariable("co2", float_type, ("z"))
nc_co2_inflow = nc_group_init.createVariable("co2_inflow", float_type, ("z"))

nc_z[:] = z[:]
nc_u[:] = u[:]
nc_u_nudge[:] = u[:]
nc_nudge_fac[:] = np.ones(kmax) * 1000
nc_v[:] = v[:]
nc_thl[:] = thl[:]
# nc_qt[:] = qt[:]
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
