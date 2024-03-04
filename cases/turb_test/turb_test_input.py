import netCDF4 as nc4
import numpy as np

# Vertical dimension and coordinates
ktot = 48
dz = 5000./ktot
z = np.arange(dz/2, ktot*dz, dz)

# Initial profiles
thl = 280. + 0.006*z

z_u = np.array([0, 270, 3000, 5000])
v_u = np.array([2.3, 8.5, 0.6, 5.7])
u = np.interp(z, z_u, v_u)

# Save all the input data to NetCDF
nc_file = nc4.Dataset('turb_test_input.nc', mode='w', datamodel='NETCDF4', clobber=False)

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
add_variable(nc_group_init, 'u', ('z'), u)

nc_file.close()