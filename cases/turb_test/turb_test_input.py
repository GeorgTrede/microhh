import netCDF4 as nc4
import numpy as np

float_type = "f8"

# Vertical dimension and coordinates
ktot = 48
dz = 5000./ktot
z = np.arange(dz/2, ktot*dz, dz)

v_thl = np.array([285.7, 291.9, 293, 297.4, 307])
z_thl = np.array([0, 400, 2000, 2500, 5000])
thl = np.interp(z, z_thl, v_thl)

z_qt = np.array([0, 400, 2000, 2500, 5000])
v_qt = np.array([6.2, 4.93, 3.61, 1, 0.3]) / 1000
qt = np.interp(z, z_qt, v_qt)

z_u = np.array([0, 270, 3000, 5000])
v_u = np.array([2.3, 8.5, 0.6, 5.7])
u = np.interp(z, z_u, v_u)

v = np.zeros(ktot)
co2 = np.zeros(ktot)

t0 = 4 * 3600
t1 = 16 * 3600
td1 = 12 * 3600
td2 = 14 * 3600

time = np.linspace(t0, t1, 32)
wthl = 0.17 * np.sin(np.pi * (time - t0) / td1)
wqt = 8.3e-5 * np.sin(np.pi * (time - t0) / td2)

# Save all the input data to NetCDF
nc_file = nc4.Dataset('turb_test_input.nc', mode='w', datamodel='NETCDF4', clobber=False) # type: ignore

def add_variable(nc_group, name, dims, data):
    """ Help function for adding a new variable """
    var = nc_group.createVariable(name, 'f8', dims)
    var[:] = data

# Create dimension and variable of vertical grid in main group:
nc_file.createDimension('z', ktot)
add_variable(nc_file, 'z', ('z'), z)

# Create a group called `init` for the initial profiles:
nc_group_init = nc_file.createGroup('init')
add_variable(nc_group_init, 'thl', ('z'), thl)
add_variable(nc_group_init, 'qt', ('z'), qt)
add_variable(nc_group_init, 'u', ('z'), u)
add_variable(nc_group_init, 'v', ('z'), v)
add_variable(nc_group_init, 'co2', ('z'), co2)
add_variable(nc_group_init, 'co2_inflow', ('z'), co2)

time_ls = np.linspace(0, 2*86400, 64)

nc_group_tdep = nc_file.createGroup("timedep")
nc_group_tdep.createDimension("time_surface", time.size)
nc_group_tdep.createDimension('time_ls', time_ls.size)

nc_time_surface = nc_group_tdep.createVariable("time_surface", float_type, ("time_surface"))
nc_time_surface[:] = time
nc_time_ls = nc_group_tdep.createVariable("time_ls", float_type, ("time_ls"))
nc_time_ls[:] = time_ls

nc_u_ls = nc_group_tdep.createVariable("u_ls", float_type, ("time_ls"))
nc_u_ls[:] = 4 * (1 + np.sin(2 * np.pi * time_ls / 3600))

nc_thl_sbot = nc_group_tdep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot = nc_group_tdep.createVariable("qt_sbot", float_type, ("time_surface"))

nc_time_surface[:] = time
nc_thl_sbot[:] = wthl
nc_qt_sbot[:] = wqt

nc_file.close()