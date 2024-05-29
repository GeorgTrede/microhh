# %%
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize

# Open the relevant data files
u_file = nc.Dataset("u.xy.nc", "r")
v_file = nc.Dataset("v.xy.nc", "r")
w_file = nc.Dataset("w.xy.nc", "r")
co2_file = nc.Dataset("co2_path.xy.nc", "r")
tke_file = nc.Dataset("test.default.0000000.nc", "r")

# Extract the relevant data
u_data = u_file.variables["u"]
v_data = v_file.variables["v"]
w_data = w_file.variables["w"]
co2_data = co2_file.variables["co2_path"]
tke_data = tke_file["default"]["tke"]
z = u_file.variables["z"][:]

# Get the x, y, and time dimensions
x = u_file.variables["xh"][:]
y = u_file.variables["y"][:]
time = u_file.variables["time"][:]  # Use only the first 25 time steps

# Create the figure and subplots
fig = plt.figure(figsize=(12, 9))
# make mosaic of 2x2 plots and one plot for TKE
axs = fig.subplot_mosaic([["u", "v"], ["w", "co2"], ["tke", "tke"]])
ax1, ax2, ax3, ax4, ax_tke = axs["u"], axs["v"], axs["w"], axs["co2"], axs["tke"]

# Create the colorbars
u_norm = Normalize(vmin=u_data[:, 0, :, :].min(), vmax=u_data[:, 0, :, :].max())
v_norm = Normalize(vmin=v_data[:, 0, :, :].min(), vmax=v_data[:, 0, :, :].max())
w_norm = Normalize(vmin=w_data[:, 0, :, :].min(), vmax=w_data[:, 0, :, :].max())
co2_norm = Normalize(vmin=co2_data[:].min(), vmax=co2_data[:].max())
tke_norm = Normalize(
    vmin=tke_data[: len(time), 2].min(), vmax=tke_data[: len(time), 2].max()
)

cbar1 = fig.colorbar(
    ax1.imshow(
        u_data[0, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=u_norm,
    ),
    ax=ax1,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar2 = fig.colorbar(
    ax2.imshow(
        v_data[0, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=v_norm,
    ),
    ax=ax2,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar3 = fig.colorbar(
    ax3.imshow(
        w_data[0, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=w_norm,
    ),
    ax=ax3,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar4 = fig.colorbar(
    ax4.imshow(
        co2_data[0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=co2_norm,
    ),
    ax=ax4,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)


# Function to update the plots
def update_xy(frame):
    print(f"Frame: {frame:4d}/{len(time)}", end="\r")
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax_tke.clear()

    # Plot the wind components at the lowest height
    ax1.imshow(
        u_data[frame, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=u_norm,
    )
    ax1.set_title(f"$u$ at $z={z[0]:.2f}$ m")
    cbar1.set_alpha(1 if ax1.images else 0)
    cbar1.set_label(r"$u$ (m/s)")
    ax2.imshow(
        v_data[frame, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=v_norm,
    )
    ax2.set_title(f"$v$ at $z={z[0]:.2f}$ m")
    cbar2.set_alpha(1 if ax2.images else 0)
    cbar2.set_label(r"$v$ (m/s)")
    ax3.imshow(
        w_data[frame, 0, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=w_norm,
    )
    ax3.set_title(f"$w$ at $z={z[0]:.2f}$ m")
    cbar3.set_alpha(1 if ax3.images else 0)
    cbar3.set_label(r"$w$ (m/s)")

    # Plot the CO2 path concentration
    ax4.imshow(
        co2_data[frame, :, :],
        cmap="viridis",
        extent=(x[0], x[-1], y[0], y[-1]),
        norm=co2_norm,
    )
    ax4.set_title(r"CO$_2$ concentration (integrated)")
    cbar4.set_alpha(1 if ax4.images else 0)
    cbar4.set_label(r"CO$_2$ (a.u.)")

    # Plot the TKE
    ax_tke.plot(time, tke_data[: len(time), 2])
    ax_tke.scatter(time[frame], tke_data[frame, 2], color="r")
    ax_tke.vlines(time[frame], 0, tke_data[frame, 2], color="r", linestyle="--")
    ax_tke.set_xlabel("Time (s)")
    ax_tke.set_ylabel(r"TKE (m$^2$/s$^2$)")
    ax_tke.set_title(f'Turbulent Kinetic Energy at $z={tke_file["z"][2]:.2f}$ m')
    ax_tke.set_ylim(tke_norm.vmin, tke_norm.vmax)

    fig.suptitle(f"Time step: {time[frame]:.2f} s", fontsize=24)

    return [ax1, ax2, ax3, ax4, ax_tke, cbar1, cbar2, cbar3, cbar4]


# Create the animation
ani = FuncAnimation(fig, update_xy, frames=len(time), interval=50, blit=False)

# Save the animation
ani.save("wind_and_co2_animation.mp4", writer="ffmpeg")
plt.close()

# Open the relevant data files
u_xz_file = nc.Dataset("u.xz.nc", "r")
v_xz_file = nc.Dataset("v.xz.nc", "r")
w_xz_file = nc.Dataset("w.xz.nc", "r")
co2_xz_file = nc.Dataset("co2.xz.nc", "r")
tke_file = nc.Dataset("test.default.0000000.nc", "r")

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
# Use only the first 25 time steps
time_xz = u_xz_file.variables["time"][500:]

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

# Create the colorbars
u_xz_norm = Normalize(
    vmin=u_xz_data[500:, :, 0, :].min(), vmax=u_xz_data[500:, :, 0, :].max()
)
v_xz_norm = Normalize(
    vmin=v_xz_data[500:, :, 0, :].min(), vmax=v_xz_data[500:, :, 0, :].max()
)
w_xz_norm = Normalize(
    vmin=w_xz_data[500:, :, 0, :].min(), vmax=w_xz_data[500:, :, 0, :].max()
)
co2_xz_norm = Normalize(
    vmin=co2_xz_data[500:, :, 0, :].min(), vmax=co2_xz_data[500:, :, 0, :].max()
)
tke_xz_norm = Normalize(
    vmin=tke_data[500 : 500 + len(time_xz), :].sum(axis=1).min(),
    vmax=tke_data[500 : 500 + len(time_xz), :].sum(axis=1).max(),
)

cbar1_xz = fig.colorbar(
    ax1_xz.imshow(
        u_xz_data[500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        origin="lower",
    ),
    ax=ax1_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar2_xz = fig.colorbar(
    ax2_xz.imshow(
        v_xz_data[500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        origin="lower",
    ),
    ax=ax2_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar3_xz = fig.colorbar(
    ax3_xz.imshow(
        w_xz_data[500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        origin="lower",
    ),
    ax=ax3_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)
cbar4_xz = fig.colorbar(
    ax4_xz.imshow(
        co2_xz_data[500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        origin="lower",
    ),
    ax=ax4_xz,
    aspect=15,
    pad=0.03,
    shrink=0.8,
)


# Function to update the plots
def update_xz(frame):
    print(f"Frame: {frame+1:4d}/{len(time_xz)}", end="\r")
    ax1_xz.clear()
    ax2_xz.clear()
    ax3_xz.clear()
    ax4_xz.clear()
    ax_tke_xz.clear()

    # Plot the wind components in the xz plane
    ax1_xz.imshow(
        u_xz_data[frame + 500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        norm=u_xz_norm,
        origin="lower",
    )
    ax1_xz.set_title(r"$u$")
    cbar1_xz.set_alpha(1 if ax1_xz.images else 0)
    cbar1_xz.set_label(r"$u$ (m/s)")

    ax2_xz.imshow(
        v_xz_data[frame + 500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        norm=v_xz_norm,
        origin="lower",
    )
    ax2_xz.set_title(r"$v$")
    cbar2_xz.set_alpha(1 if ax2_xz.images else 0)
    cbar2_xz.set_label(r"$v$ (m/s)")

    ax3_xz.imshow(
        w_xz_data[frame + 500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        norm=w_xz_norm,
        origin="lower",
    )
    ax3_xz.set_title(r"$w$")
    cbar3_xz.set_alpha(1 if ax3_xz.images else 0)
    cbar3_xz.set_label(r"$w$ (m/s)")

    ax4_xz.imshow(
        co2_xz_data[frame + 500, :, 0, :],
        cmap="viridis",
        extent=(x_xz[0], x_xz[-1], z_xz[-1], z_xz[0]),
        norm=co2_xz_norm,
        origin="lower",
    )
    ax4_xz.set_title(r"CO$_2$ concentration")
    cbar4_xz.set_alpha(1 if ax4_xz.images else 0)
    cbar4_xz.set_label(r"CO$_2$ (a.u.)")

    # Plot the TKE
    ax_tke_xz.plot(time_xz, tke_data[500 : 500 + len(time_xz), :].sum(axis=1))
    ax_tke_xz.scatter(time_xz[frame], tke_data[frame + 500, :].sum(), color="r", s=10)
    ax_tke_xz.vlines(
        time_xz[frame],
        0,
        tke_data[frame + 500, :].sum(),
        color="r",
        linestyle="--",
        linewidth=1,
    )
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
ani_xz = FuncAnimation(fig, update_xz, frames=len(time_xz), interval=50, blit=False)

# Save the animation
ani_xz.save("wind_and_co2_xz_animation.mp4", writer="ffmpeg")
plt.close()
