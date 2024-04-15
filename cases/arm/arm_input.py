import numpy as np
import netCDF4 as nc

float_type = "f8"

# Get number of vertical levels and size from .ini file
with open("arm.ini") as f:
    for line in f:
        if line.split("=")[0] == "ktot":
            kmax = int(line.split("=")[1])
        if line.split("=")[0] == "zsize":
            zsize = float(line.split("=")[1])

dz = zsize / kmax

# set the height
z = np.linspace(0.5 * dz, zsize - 0.5 * dz, kmax)
thl = np.zeros(np.size(z))
qt = np.zeros(np.size(z))
u = np.zeros(np.size(z))
ug = np.zeros(np.size(z))

for k in range(kmax):
    # temperature
    if z[k] <= 50.0:
        thl[k] = 299.0 + (z[k]) * (301.5 - 299.0) / (50.0)
        qt[k] = 15.20 + (z[k]) * (15.17 - 15.20) / (50.0)
    elif z[k] <= 350.0:
        thl[k] = 301.5 + (z[k] - 50.0) * (302.5 - 301.5) / (350.0 - 50.0)
        qt[k] = 15.17 + (z[k] - 50.0) * (14.98 - 15.17) / (350.0 - 50.0)
    elif z[k] <= 650.0:
        thl[k] = 302.5 + (z[k] - 350.0) * (303.53 - 302.5) / (650.0 - 350.0)
        qt[k] = 14.98 + (z[k] - 350.0) * (14.80 - 14.98) / (650.0 - 350.0)
    elif z[k] <= 700.0:
        thl[k] = 303.53 + (z[k] - 650.0) * (303.7 - 303.53) / (700.0 - 650.0)
        qt[k] = 14.80 + (z[k] - 650.0) * (14.70 - 14.80) / (700.0 - 650.0)
    elif z[k] <= 1300.0:
        thl[k] = 303.7 + (z[k] - 700.0) * (307.13 - 303.7) / (1300.0 - 700.0)
        qt[k] = 14.70 + (z[k] - 700.0) * (13.50 - 14.80) / (1300.0 - 700.0)
    elif z[k] <= 2500.0:
        thl[k] = 307.13 + (z[k] - 1300.0) * (314.0 - 307.13) / (2500.0 - 1300.0)
        qt[k] = 13.50 + (z[k] - 1300.0) * (3.00 - 13.50) / (2500.0 - 1300.0)
    elif z[k] <= 5500.0:
        thl[k] = 314.0 + (z[k] - 2500.0) * (343.2 - 314.0) / (5500.0 - 2500.0)
        qt[k] = 3.00

    # u-wind component
    u[:] = 10.0

    # ug-wind component
    ug[k] = 10.0

# normalize profiles to SI
qt /= 1000.0  # g to kg

# set the time series
time_surface = np.array([0.0, 4.0, 6.5, 7.5, 10.0, 12.5, 14.5])

H = np.array([-30.0, 90.0, 140.0, 140.0, 100.0, -10.0, -10])
LE = np.array([5.0, 250.0, 450.0, 500.0, 420.0, 180.0, 0])

advthl = np.array([0.0, 0.0, 0.0, -0.08, -0.16, -0.16])
radthl = np.array([-0.125, 0.0, 0.0, 0.0, 0.0, -0.1])
advqt = np.array([0.08, 0.02, -0.04, -0.10, -0.16, -0.30])

time_ls = np.array([0.0, 3.0, 6.0, 9.0, 12.0, 14.5])
thlls = np.zeros((time_ls.size, kmax))
qtls = np.zeros((time_ls.size, kmax))

# large scale forcings
for n in range(time_ls.size):
    tendthl = advthl[n] + radthl[n]
    tendqt = advqt[n]
    for k in range(kmax):
        if z[k] <= 1000.0:
            thlls[n, k] = tendthl
            qtls[n, k] = tendqt
        else:
            thlls[n, k] = tendthl - (z[k] - 1000.0) * (tendthl) / (5500.0 - 1000.0)
            qtls[n, k] = tendqt - (z[k] - 1000.0) * (tendqt) / (5500.0 - 1000.0)

time_ls *= 3600.0  # h to s
thlls /= 3600.0  # h to s
qtls /= 3600.0  # h to s
qtls /= 1000.0  # g to kg

# Calculate the surface fluxes in the correct units.
Rd = 287.0
cp = 1005.0
Lv = 2.5e6
p0 = 97000.0
rho = p0 / (Rd * thl[0] * (1.0 + 0.61 * qt[0]))
time_surface *= 3600.0  # h to s
sbotthl = H / (rho * cp)
sbotqt = LE / (rho * Lv)

# Save all the input data to NetCDF
nc_file = nc.Dataset("arm_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup("init")

nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt = nc_group_init.createVariable("qt", float_type, ("z"))
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_ug = nc_group_init.createVariable("u_geo", float_type, ("z"))
nc_thl[:] = thl[:]
nc_qt[:] = qt[:]
nc_u[:] = u[:]
nc_ug[:] = ug[:]

# Create a group called "timedep" for the timedep.
nc_group_timedep = nc_file.createGroup("timedep")
nc_group_timedep.createDimension("time_surface", time_surface.size)
nc_group_timedep.createDimension("time_ls", time_ls.size)

nc_time_surface = nc_group_timedep.createVariable(
    "time_surface", float_type, ("time_surface")
)
nc_thl_sbot = nc_group_timedep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot = nc_group_timedep.createVariable("qt_sbot", float_type, ("time_surface"))
nc_time_surface[:] = time_surface[:]
nc_thl_sbot[:] = sbotthl[:]
nc_qt_sbot[:] = sbotqt[:]

nc_time_ls = nc_group_timedep.createVariable("time_ls", float_type, ("time_ls"))
nc_thl_ls = nc_group_timedep.createVariable("thl_ls", float_type, ("time_ls", "z"))
nc_qt_ls = nc_group_timedep.createVariable("qt_ls", float_type, ("time_ls", "z"))
nc_time_ls[:] = time_ls[:]
nc_thl_ls[:, :] = thlls[:, :]
nc_qt_ls[:, :] = qtls[:, :]

nc_file.close()
