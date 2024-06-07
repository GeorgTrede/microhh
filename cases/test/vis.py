# %%
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
import sys
import warnings


# parse command line arguments and:
# - look for -f FACTOR in arguments (default 3)
# - look for -s START_IDX in arguments (default 1000)

FACTOR = 3
START_IDX = 500
for i, arg in enumerate(sys.argv):
    if arg == "-f":
        FACTOR = int(sys.argv[i + 1])
    if arg == "-s":
        START_IDX = int(sys.argv[i + 1])

print(f"FACTOR: {FACTOR}, START_IDX: {START_IDX}")

# Open the relevant data files
u_xy_file = nc.Dataset("u.xy.nc", "r")  # type: ignore
v_xy_file = nc.Dataset("v.xy.nc", "r")  # type: ignore
w_xy_file = nc.Dataset("w.xy.nc", "r")  # type: ignore
co2_xy_file = nc.Dataset("co2_path.xy.nc", "r")  # type: ignore
tke_file = nc.Dataset("test.default.0000000.nc", "r")  # type: ignore

# Extract the relevant data
u_xy_data = u_xy_file.variables["u"]
v_xy_data = v_xy_file.variables["v"]
w_xy_data = w_xy_file.variables["w"]
co2_xy_data = co2_xy_file.variables["co2_path"]
tke_data = tke_file["default"]["tke"]
z_xy = u_xy_file.variables["z"][:]

# Get the x, y, and time dimensions
x_xy = u_xy_file.variables["xh"][:]
y_xy = u_xy_file.variables["y"][:]
time_xy = u_xy_file.variables["time"][START_IDX:]

# Create the figure and subplots
fig = plt.figure(figsize=(12, 9))
# make mosaic of 2x2 plots and one plot for TKE
axs = fig.subplot_mosaic([["u_xy", "v_xy"], ["w_xy", "co2_xy"], ["tke_xy", "tke_xy"]])
ax1_xy, ax2_xy, ax3_xy, ax4_xy, ax_tke_xy = (
    axs["u_xy"],
    axs["v_xy"],
    axs["w_xy"],
    axs["co2_xy"],
    axs["tke_xy"],
)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    # Create the colorbars
    u_min, u_max = np.percentile(u_xy_data[START_IDX:, 0, :, :], [1, 99])
    u_xy_norm = Normalize(vmin=u_min, vmax=u_max)

    v_min, v_max = np.percentile(v_xy_data[START_IDX:, 0, :, :], [1, 99])
    v_xy_norm = Normalize(vmin=v_min, vmax=v_max)

    w_min, w_max = np.percentile(w_xy_data[START_IDX:, 0, :, :], [1, 99])
    w_xy_norm = Normalize(vmin=w_min, vmax=w_max)

    co2_min, co2_max = np.percentile(co2_xy_data[START_IDX:, :, :], [1, 99.9])
    co2_xy_norm = Normalize(vmin=co2_min, vmax=co2_max)

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
cbar4_xy = fig.colorbar(
    ax4_xy.imshow(
        co2_xy_data[START_IDX, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        origin="lower",
    ),
    ax=ax4_xy,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar4_xy.set_label(r"CO$_2$ (a.u.)")


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
        u_xy_data[frame + START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=u_xy_norm,
        origin="lower",
    )
    ax1_xy.set_title(r"$u$")
    ax1_xy.set_xlabel("x (m)")
    ax1_xy.set_ylabel("y (m)")

    ax2_xy.imshow(
        v_xy_data[frame + START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=v_xy_norm,
        origin="lower",
    )
    ax2_xy.set_title(r"$v$")
    ax2_xy.set_xlabel("x (m)")
    ax2_xy.set_ylabel("y (m)")

    ax3_xy.imshow(
        w_xy_data[frame + START_IDX, 0, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=w_xy_norm,
        origin="lower",
    )
    ax3_xy.set_title(r"$w$")
    ax3_xy.set_xlabel("x (m)")
    ax3_xy.set_ylabel("y (m)")

    ax4_xy.imshow(
        co2_xy_data[frame + START_IDX, :, :],
        cmap="turbo",
        extent=(x_xy[0], x_xy[-1], y_xy[0], y_xy[-1]),
        norm=co2_xy_norm,
        origin="lower",
    )
    ax4_xy.set_title(r"CO$_2$ concentration")
    ax4_xy.set_xlabel("x (m)")
    ax4_xy.set_ylabel("y (m)")

    # Plot the TKE
    ax_tke_xy.plot(time_xy, tke_data[START_IDX : START_IDX + len(time_xy), 0])
    ax_tke_xy.scatter(time_xy[frame], tke_data[frame + START_IDX, 0], color="r", s=10)
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
        cbar4_xy,
    ]


# Create the animation
ani_xy = FuncAnimation(
    fig, update_xy, frames=len(time_xy) // FACTOR, interval=50, blit=False
)

# Save the animation
ani_xy.save("wind_and_co2_xy_animation.mp4", writer="ffmpeg")
plt.close()


# Open the relevant data files
u_xz_file = nc.Dataset("u.xz.nc", "r")  # type: ignore
v_xz_file = nc.Dataset("v.xz.nc", "r")  # type: ignore
w_xz_file = nc.Dataset("w.xz.nc", "r")  # type: ignore
co2_xz_file = nc.Dataset("co2.xz.nc", "r")  # type: ignore
tke_file = nc.Dataset("test.default.0000000.nc", "r")  # type: ignore

# Extract the relevant data
u_xz_data = u_xz_file.variables["u"]
v_xz_data = v_xz_file.variables["v"]
w_xz_data = w_xz_file.variables["w"]
co2_xz_data = co2_xz_file.variables["co2"]
tke_data = tke_file["default"]["tke"]
z_xz = u_xz_file.variables["z"][:]

# Get the x, z, and time dimensions
x_xz = u_xz_file.variables["xh"][:]
z_xz = u_xz_file.variables["z"][:]
time_xz = u_xz_file.variables["time"][START_IDX:]

# Create the figure and subplots
fig = plt.figure(figsize=(12, 9))
# make mosaic of 2x2 plots and one plot for TKE
axs = fig.subplot_mosaic([["u_xz", "v_xz"], ["w_xz", "co2_xz"], ["tke_xz", "tke_xz"]])
ax1_xz, ax2_xz, ax3_xz, ax4_xz, ax_tke_xz = (
    axs["u_xz"],
    axs["v_xz"],
    axs["w_xz"],
    axs["co2_xz"],
    axs["tke_xz"],
)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    # Create the colorbars
    u_min, u_max = np.percentile(u_xz_data[START_IDX:, :, 0, :], [1, 99])
    u_xz_norm = Normalize(vmin=u_min, vmax=u_max)

    v_min, v_max = np.percentile(v_xz_data[START_IDX:, :, 0, :], [1, 99])
    v_xz_norm = Normalize(vmin=v_min, vmax=v_max)

    w_min, w_max = np.percentile(w_xz_data[START_IDX:, :, 0, :], [1, 99])
    w_xz_norm = Normalize(vmin=w_min, vmax=w_max)

    co2_min, co2_max = np.percentile(co2_xz_data[START_IDX:, :, 0, :], [1, 99.9])
    co2_xz_norm = Normalize(vmin=co2_min, vmax=co2_max)

    tke_min, tke_max = np.percentile(
        tke_data[START_IDX : START_IDX + len(time_xz), :].sum(axis=1), [0, 100]
    )
# extend tke range by 10% of the range
tke_min -= 0.1 * (tke_max - tke_min)
tke_max += 0.1 * (tke_max - tke_min)
tke_xz_norm = Normalize(vmin=tke_min, vmax=tke_max)

cbar1_xz = fig.colorbar(
    ax1_xz.imshow(
        u_xz_data[START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        origin="lower",
    ),
    ax=ax1_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar1_xz.set_label(r"$u$ (m/s)")
cbar2_xz = fig.colorbar(
    ax2_xz.imshow(
        v_xz_data[START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        origin="lower",
    ),
    ax=ax2_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar2_xz.set_label(r"$v$ (m/s)")
cbar3_xz = fig.colorbar(
    ax3_xz.imshow(
        w_xz_data[START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        origin="lower",
    ),
    ax=ax3_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar3_xz.set_label(r"$w$ (m/s)")
cbar4_xz = fig.colorbar(
    ax4_xz.imshow(
        co2_xz_data[START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        origin="lower",
    ),
    ax=ax4_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar4_xz.set_label(r"CO$_2$ (a.u.)")


# Function to update the plots
def update_xz(frame):
    frame *= FACTOR
    print(f"Frame: {frame+1:4d}/{len(time_xz)}", end="\r")
    ax1_xz.clear()
    ax2_xz.clear()
    ax3_xz.clear()
    ax4_xz.clear()
    ax_tke_xz.clear()

    # Plot the wind components in the xz plane
    ax1_xz.imshow(
        u_xz_data[frame + START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        norm=u_xz_norm,
        origin="lower",
    )
    ax1_xz.set_title(r"$u$")
    ax1_xz.set_xlabel("x (m)")
    ax1_xz.set_ylabel("z (m)")

    ax2_xz.imshow(
        v_xz_data[frame + START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        norm=v_xz_norm,
        origin="lower",
    )
    ax2_xz.set_title(r"$v$")
    ax2_xz.set_xlabel("x (m)")
    ax2_xz.set_ylabel("z (m)")

    ax3_xz.imshow(
        w_xz_data[frame + START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        norm=w_xz_norm,
        origin="lower",
    )
    ax3_xz.set_title(r"$w$")
    ax3_xz.set_xlabel("x (m)")
    ax3_xz.set_ylabel("z (m)")

    ax4_xz.imshow(
        co2_xz_data[frame + START_IDX, :, 0, :],
        cmap="turbo",
        extent=(x_xz[0], x_xz[-1], z_xz[0], z_xz[-1]),
        norm=co2_xz_norm,
        origin="lower",
    )
    ax4_xz.set_title(r"CO$_2$ concentration")
    ax4_xz.set_xlabel("x (m)")
    ax4_xz.set_ylabel("z (m)")

    # Plot the TKE
    ax_tke_xz.plot(
        time_xz, tke_data[START_IDX : START_IDX + len(time_xz), :].sum(axis=1)
    )
    ax_tke_xz.scatter(
        time_xz[frame], tke_data[frame + START_IDX, :].sum(), color="r", s=10
    )
    ax_tke_xz.vlines(time_xz[frame], 0, 1000, color="r", linestyle="--", linewidth=1)
    ax_tke_xz.set_xlabel("Time (s)")
    ax_tke_xz.set_ylabel(r"TKE (m$^2$/s$^2$)")
    ax_tke_xz.set_title(r"Total Turbulent Kinetic Energy")
    ax_tke_xz.set_ylim(tke_xz_norm.vmin, tke_xz_norm.vmax)

    fig.suptitle(f"Time step: {time_xz[frame]:.2f} s", fontsize=24)

    return [
        ax1_xz,
        ax2_xz,
        ax3_xz,
        ax4_xz,
        ax_tke_xz,
        cbar1_xz,
        cbar2_xz,
        cbar3_xz,
        cbar4_xz,
    ]


# Create the animation
ani_xz = FuncAnimation(
    fig, update_xz, frames=len(time_xz) // FACTOR, interval=50, blit=False
)

# Save the animation
ani_xz.save("wind_and_co2_xz_animation.mp4", writer="ffmpeg")
plt.close()
