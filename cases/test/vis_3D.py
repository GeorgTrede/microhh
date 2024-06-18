import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
import sys
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# parse command line arguments and:
# - look for -f FACTOR in arguments (default 3)
# - look for -s START_IDX in arguments (default 1000)
# - look for -z Z_LEVEL in arguments (default 0)

FACTOR = 3
START_IDX = 500
Z_LEVEL = 0
for i, arg in enumerate(sys.argv):
    if arg == "-f":
        FACTOR = int(sys.argv[i + 1])
    if arg == "-s":
        START_IDX = int(sys.argv[i + 1])
    if arg == "-z":
        Z_LEVEL = int(sys.argv[i + 1])

print(f"FACTOR: {FACTOR}, START_IDX: {START_IDX}, Z_LEVEL: {Z_LEVEL}")
print()
print("Creating xy plane and 3D animation")


def update_projection(ax, name, projection="3d", fig=None):
    axi = ax[name]
    if fig is None:
        fig = plt.gcf()
    rows, cols, start, stop = axi.get_subplotspec().get_geometry()
    ax[name].remove()
    ax[name] = fig.add_subplot(rows, cols, start + 1, projection=projection)


# Open the relevant data files
u_xy_file = nc.Dataset("u.xy.nc", "r")  # type: ignore
v_xy_file = nc.Dataset("v.xy.nc", "r")  # type: ignore
w_xy_file = nc.Dataset("w.xy.nc", "r")  # type: ignore
co2_xy_file = nc.Dataset("co2.xy.nc", "r")  # type: ignore
tke_file = nc.Dataset("test.default.0000000.nc", "r")  # type: ignore

# Extract the relevant data
u_xy_data = u_xy_file.variables["u"]
v_xy_data = v_xy_file.variables["v"]
w_xy_data = w_xy_file.variables["w"]
co2_xy_data = co2_xy_file.variables["co2"]
tke_data = tke_file["default"]["tke"]
z_xy = co2_xy_file.variables["z"][:]

# Get the x, y, and time dimensions
x_xy = u_xy_file.variables["xh"][:]
y_xy = u_xy_file.variables["y"][:]
time_xy = u_xy_file.variables["time"][START_IDX:]

# Create the figure and subplots
fig = plt.figure(figsize=(12, 9))
# make mosaic of 2x2 plots and one plot for TKE
axs = fig.subplot_mosaic(
    [["u_xy", "v_xy"], ["w_xy", "co2_plume"], ["tke_xy", "tke_xy"]], width_ratios=[1, 1]
)

update_projection(axs, "co2_plume", "3d", fig=fig)

ax1_xy, ax2_xy, ax3_xy, ax4_xy, ax_tke_xy = (
    axs["u_xy"],
    axs["v_xy"],
    axs["w_xy"],
    axs["co2_plume"],
    axs["tke_xy"],
)

# Create the colorbars
u_min, u_max = np.percentile(u_xy_data[START_IDX:, Z_LEVEL, :, :], [1, 99])
u_xy_norm = Normalize(vmin=u_min, vmax=u_max)

v_min, v_max = np.percentile(v_xy_data[START_IDX:, Z_LEVEL, :, :], [1, 99])
v_xy_norm = Normalize(vmin=v_min, vmax=v_max)

w_min, w_max = np.percentile(w_xy_data[START_IDX:, Z_LEVEL, :, :], [1, 99])
w_xy_norm = Normalize(vmin=w_min, vmax=w_max)

co2_min, co2_max = np.percentile(co2_xy_data[START_IDX:, :, :, :], [10, 100])

tke_min, tke_max = np.percentile(
    tke_data[START_IDX : START_IDX + len(time_xy), 0], [0, 100]
)
# extend tke range by 10% of the range
tke_min -= 0.1 * (tke_max - tke_min)
tke_max += 0.1 * (tke_max - tke_min)
tke_xy_norm = Normalize(vmin=tke_min, vmax=tke_max)

cbar1_xy = fig.colorbar(
    ax1_xy.imshow(
        u_xy_data[START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        origin="lower",
    ),
    ax=ax1_xy,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar1_xy.set_label(r"$u$ (m/s)")
cbar2_xy = fig.colorbar(
    ax2_xy.imshow(
        v_xy_data[START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        origin="lower",
    ),
    ax=ax2_xy,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar2_xy.set_label(r"$v$ (m/s)")
cbar3_xy = fig.colorbar(
    ax3_xy.imshow(
        w_xy_data[START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        origin="lower",
    ),
    ax=ax3_xy,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar3_xy.set_label(r"$w$ (m/s)")


# Function to update the plots
def update_xy(frame):
    frame *= FACTOR
    print(f"Frame: {frame+1:4d}/{len(time_xy)}", end="\r")
    ax1_xy.clear()
    ax2_xy.clear()
    ax3_xy.clear()
    ax4_xy.clear()
    ax_tke_xy.clear()

    # Plot the wind components in the xy plane
    ax1_xy.imshow(
        u_xy_data[frame + START_IDX, Z_LEVEL, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=u_xy_norm,
        origin="lower",
    )
    ax1_xy.set_title(r"$u$")
    ax1_xy.set_xlabel("x (m)")
    ax1_xy.set_ylabel("y (m)")

    ax2_xy.imshow(
        v_xy_data[frame + START_IDX, Z_LEVEL, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=v_xy_norm,
        origin="lower",
    )
    ax2_xy.set_title(r"$v$")
    ax2_xy.set_xlabel("x (m)")
    ax2_xy.set_ylabel("y (m)")

    ax3_xy.imshow(
        w_xy_data[frame + START_IDX, Z_LEVEL, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=w_xy_norm,
        origin="lower",
    )
    ax3_xy.set_title(r"$w$")
    ax3_xy.set_xlabel("x (m)")
    ax3_xy.set_ylabel("y (m)")

    # Plot the 3D CO2 plume
    co2_data = co2_xy_data[frame + START_IDX, :, :, :]
    ax4_xy.clear()
    ax4_xy.set_xlabel("x (m)")
    ax4_xy.set_ylabel("y (m)")
    ax4_xy.set_zlabel("z (m)")  # type: ignore
    for z in range(co2_data.shape[0]):
        co2_slice = co2_data[z, :, :]
        mask = co2_slice > co2_min
        x_coords, y_coords = np.meshgrid(x_xy, y_xy)
        ax4_xy.scatter(
            x_coords[mask],
            y_coords[mask],
            z_xy[z] * np.ones_like(x_coords[mask]),
            c=co2_slice[mask],
            alpha=0.025,
            cmap="gray_r",
            vmin=co2_min,
            vmax=co2_max,
        )
    ax4_xy.set_title(r"CO$_2$ concentration")
    # set axis limits
    ax4_xy.set_xlim(x_xy[0], x_xy[-1])
    ax4_xy.set_ylim(y_xy[0], y_xy[-1])
    ax4_xy.set_zlim(0, 5000)  # type: ignore
    # set aspect ratio
    ax4_xy.set_box_aspect([12, 3.2, 5.0])  # type: ignore

    # Plot the TKE
    ax_tke_xy.plot(time_xy, tke_data[START_IDX : START_IDX + len(time_xy), 0])
    ax_tke_xy.scatter(
        time_xy[frame], tke_data[frame + START_IDX, Z_LEVEL], color="r", s=10
    )
    ax_tke_xy.vlines(time_xy[frame], 0, 1000, color="r", linestyle="--", linewidth=1)
    ax_tke_xy.set_xlabel("Time (s)")
    ax_tke_xy.set_ylabel(r"TKE (m$^2$/s$^2$)")
    ax_tke_xy.set_title(r"Total Turbulent Kinetic Energy")
    ax_tke_xy.set_ylim(tke_xy_norm.vmin, tke_xy_norm.vmax)

    fig.suptitle(f"Time step: {time_xy[frame]:.2f} s", fontsize=24)

    return [
        ax1_xy,
        ax2_xy,
        ax3_xy,
        ax4_xy,
        ax_tke_xy,
        cbar1_xy,
        cbar2_xy,
        cbar3_xy,
    ]


# Create the animation
ani_xy = FuncAnimation(
    fig, update_xy, frames=len(time_xy) // FACTOR, interval=50, blit=False
)

# Save the animation
ani_xy.save("wind_and_co2_xy_3D_animation.mp4", writer="ffmpeg")
plt.close()

print()
print("Saved wind_and_co2_xy_3D_animation.mp4")
print()
