import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize

# Open the relevant data files
u_xz_file = nc.Dataset('u.xz.nc', 'r')
co2_xz_file = nc.Dataset('co2.xz.nc', 'r')
tke_file = nc.Dataset('test.default.0000000.nc', 'r')

# Extract the relevant data
u_xz_data = u_xz_file.variables['u']
co2_xz_data = co2_xz_file.variables['co2']
tke_data = tke_file['default']['tke']
z_xz = u_xz_file.variables['z'][:]
x_xz = u_xz_file.variables['xh'][:]

# Get the time dimensions
START_IDX = 1000
FACTOR = 10
time_xz = u_xz_file.variables['time'][START_IDX:]

# Create the figure and subplots
fig = plt.figure(figsize=(12, 9))
# make mosaic of 2x2 plots and one plot for TKE
axs = fig.subplot_mosaic([['u_xz', 'zoom_u_xz'], ['co2_xz', 'zoom_co2_xz'], ['tke_xz', 'tke_xz']], width_ratios=[1, 1])
ax1_xz, ax2_xz, ax3_xz, ax4_xz, ax_tke_xz = axs['u_xz'], axs['co2_xz'], axs['zoom_u_xz'], axs['zoom_co2_xz'], axs['tke_xz']

# Create the colorbars
u_min, u_max = np.percentile(u_xz_data[START_IDX:, :, 0, :], [1, 99])
u_xz_norm = Normalize(vmin=u_min, vmax=u_max)

co2_min, co2_max = np.percentile(co2_xz_data[START_IDX:, :, 0, :], [1, 99.9])
co2_xz_norm = Normalize(vmin=co2_min, vmax=co2_max)

tke_min, tke_max = np.percentile(tke_data[START_IDX:, :].sum(axis=1), [0, 100])
tke_min -= 0.1 * (tke_max - tke_min)
tke_max += 0.1 * (tke_max - tke_min)
tke_xz_norm = Normalize(vmin=tke_min, vmax=tke_max)

cbar1_xz = fig.colorbar(ax1_xz.imshow(u_xz_data[START_IDX, :, 0, :], cmap='turbo', extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]), origin='lower'), ax=ax1_xz, aspect=15, pad=0.03, shrink=0.8)
cbar1_xz.set_label(r'$u$ (m/s)')
cbar2_xz = fig.colorbar(ax2_xz.imshow(co2_xz_data[START_IDX, :, 0, :], cmap='turbo', extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]), origin='lower'), ax=ax2_xz, aspect=15, pad=0.03, shrink=0.8)
cbar2_xz.set_label(r'CO$_2$ (a.u.)')

# Function to update the plots
def update_xz(frame):
    frame *= FACTOR
    print(f'Frame: {frame+1:4d}/{len(time_xz)}', end='\r')
    ax1_xz.clear()
    ax2_xz.clear()
    ax3_xz.clear()
    ax4_xz.clear()
    ax_tke_xz.clear()

    # Plot the u component and CO2 concentration in the full domain
    im1 = ax1_xz.imshow(u_xz_data[frame+START_IDX, :, 0, :], cmap='turbo', extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]), norm=u_xz_norm, origin='lower')
    ax1_xz.set_title(r'$u$')
    ax1_xz.set_xlabel('x (m)')
    ax1_xz.set_ylabel('z (m)')

    im2 = ax2_xz.imshow(co2_xz_data[frame+START_IDX, :, 0, :], cmap='turbo', extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]), norm=co2_xz_norm, origin='lower')
    ax2_xz.set_title(r'CO$_2$ concentration')
    ax2_xz.set_xlabel('x (m)')
    ax2_xz.set_ylabel('z (m)')

    # Plot the zoomed-in views
    z_start = 0
    z_end = int(0.4 * u_xz_data.shape[1])
    x_start = 0
    x_end = int(0.4 * u_xz_data.shape[3])

    ax3_xz.imshow(u_xz_data[frame+START_IDX, z_start:z_end, 0, x_start:x_end], cmap='turbo', extent=(x_xz[x_start], x_xz[x_end-1], z_xz[z_start], z_xz[z_end-1]), norm=u_xz_norm, origin='lower')
    ax3_xz.set_title(r'$u$ (zoomed)')
    ax3_xz.set_xlabel('x (m)')
    ax3_xz.set_ylabel('z (m)')

    ax4_xz.imshow(co2_xz_data[frame+START_IDX, z_start:z_end, 0, x_start:x_end], cmap='turbo', extent=(x_xz[x_start], x_xz[x_end-1], z_xz[z_start], z_xz[z_end-1]), norm=co2_xz_norm, origin='lower')
    ax4_xz.set_title(r'CO$_2$ concentration (zoomed)')
    ax4_xz.set_xlabel('x (m)')
    ax4_xz.set_ylabel('z (m)')

    # Plot the TKE
    ax_tke_xz.plot(time_xz, tke_data[START_IDX:, :].sum(axis=1))
    ax_tke_xz.scatter(time_xz[frame], tke_data[frame+START_IDX, :].sum(), color='r', s=10)
    ax_tke_xz.vlines(time_xz[frame], 0, 1000, color='r', linestyle='--', linewidth=1)
    ax_tke_xz.set_xlabel('Time (s)')
    ax_tke_xz.set_ylabel(r'TKE (m$^2$/s$^2$)')
    ax_tke_xz.set_title(r'Total Turbulent Kinetic Energy')
    ax_tke_xz.set_ylim(tke_xz_norm.vmin, tke_xz_norm.vmax)

    # Draw a rectangle in the full domain plots to indicate the zoomed region
    ax1_xz.add_patch(plt.Rectangle((x_xz[x_start], z_xz[z_start]), x_xz[x_end-1]-x_xz[x_start], z_xz[z_end-1]-z_xz[z_start], fill=False, color='r', linewidth=2))
    ax2_xz.add_patch(plt.Rectangle((x_xz[x_start], z_xz[z_start]), x_xz[x_end-1]-x_xz[x_start], z_xz[z_end-1]-z_xz[z_start], fill=False, color='r', linewidth=2))

    fig.suptitle(f'Time step: {time_xz[frame]:.2f} s', fontsize=24)

    return [ax1_xz, ax2_xz, ax3_xz, ax4_xz, ax_tke_xz, cbar1_xz, cbar2_xz]

# Create the animation
ani_xz = FuncAnimation(fig, update_xz, frames=len(time_xz)//FACTOR, interval=50, blit=False)

# Save the animation
ani_xz.save('wind_and_co2_xz_zoom_animation.mp4', writer='ffmpeg')
plt.close()

