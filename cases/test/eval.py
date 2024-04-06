# %%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import sys
import multiprocessing

# set plt font size
plt.rcParams.update({'font.size': 16})

show = plt.close
# show = plt.show

if len(sys.argv) > 1:
    folder = sys.argv[1]
else:
    folder = "jaenschwalde"

# %%
# res = "256_64_192"
# res = "128_32_96"
res = "64_16_48"
uflux = "X"
itot = 128
jtot = 32
ktot = 96
# read res and uflux from .ini file
with open(f"{folder}.ini") as f:
    for line in f:
        if line.split("=")[0] == "itot":
            itot = line.split("=")[1].strip()
        if line.split("=")[0] == "jtot":
            jtot = line.split("=")[1].strip()
        if line.split("=")[0] == "ktot":
            ktot = line.split("=")[1].strip()
        if line.split("=")[0] == "uflux":
            uflux = str(round(float(line.split("=")[1].strip())))
    res = f"{itot}_{jtot}_{ktot}"

print(f"res={res}, uflux={uflux}")


if "--with-frames" in sys.argv:
    no_frames = False
else:
    no_frames = True

os.chdir("/home/georg.trede/MasterThesis/env/microhh/cases/{}/snaps_{}_uflux{}".format(folder, res, uflux))

# %%
ds_co2_path_xy = nc.Dataset("nc_files/co2_path.xy.nc") # type: ignore
ds_co2_xy = nc.Dataset("nc_files/co2.xy.nc") # type: ignore
ds_u_xy = nc.Dataset("nc_files/u.xy.nc") # type: ignore
ds_v_xy = nc.Dataset("nc_files/v.xy.nc") # type: ignore
ds_w_xy = nc.Dataset("nc_files/w.xy.nc") # type: ignore

# %%
co2_path_xy = np.array(ds_co2_path_xy.variables["co2_path"][:])
co2_xy = np.array(ds_co2_xy.variables["co2"][:])
u_xy = np.array(ds_u_xy.variables["u"][:])
v_xy = np.array(ds_v_xy.variables["v"][:])
w_xy = np.array(ds_w_xy.variables["w"][:])

z_xy = np.array(ds_co2_xy.variables["z"][:])

# %%
# save co2_path, co2, u, v, w
np.save("npy_files/co2_path.npy", co2_path_xy)
for i in range(co2_xy.shape[1]):  # co2.shape[1] should give the number of heights
    np.save(f'npy_files/co2_h{i+1}_xy.npy', co2_xy[:, i, :, :])
    np.save(f'npy_files/u_h{i+1}_xy.npy', u_xy[:, i, :, :])
    np.save(f'npy_files/v_h{i+1}_xy.npy', v_xy[:, i, :, :])
    np.save(f'npy_files/w_h{i+1}_xy.npy', w_xy[:, i, :, :])

# %%
# load co2_path, co2, u, v, w from npy files
# co2_path = np.load("npy_files/co2_path.npy")
# co2 = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
# u = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
# v = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
# w = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))

# for i in range(co2.shape[1]):  # co2.shape[1] should give the number of heights
#     co2[:, i, :, :] = np.load(f'npy_files/co2_h{i+1}_xy.npy')
#     u[:, i, :, :] = np.load(f'npy_files/u_h{i+1}_xy.npy')
#     v[:, i, :, :] = np.load(f'npy_files/v_h{i+1}_xy.npy')
#     w[:, i, :, :] = np.load(f'npy_files/w_h{i+1}_xy.npy')

# %%
if not no_frames:
    z_level = 0

    co2_min, co2_max = np.min(co2_path_xy), np.max(co2_path_xy)
    u_min, u_max = np.min(u_xy[:, z_level]), np.max(u_xy[:, z_level])
    u_range = u_max - u_min
    u_min, u_max = u_min + 0.1 * u_range, u_max - 0.1 * u_range
    w_min, w_max = np.min(w_xy[:, z_level]), np.max(w_xy[:, z_level])
    w_range = w_max - w_min
    w_min, w_max = w_min + 0.1 * w_range, w_max - 0.1 * w_range

    # %%
    def process_frame(i, co2_slice, u_slice, w_slice, z_level, z, uflux, co2_max, u_min, u_max, w_min, w_max):
        print(f"\rFrame {i+1:3d}/{len(co2_path_xy):3d}", end="")
        plt.figure(figsize=(15, 12))
        plt.suptitle(f"Simulation with uflux={uflux} m/s\nt={i*100} s")
        plt.subplot(3, 1, 1)
        plt.imshow(co2_slice, vmin=0, vmax=co2_max*0.8, extent=[0, 12.8, 0, 3.2])
        plt.title("CO2 (integrated)")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        plt.colorbar(shrink=0.7)
        plt.subplot(3, 1, 2)
        plt.imshow(u_slice, vmin=u_min, vmax=u_max, extent=[0, 12.8, 0, 3.2])
        plt.title(f"u (at z={int(z[z_level])}m, in m/s)")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        plt.colorbar(shrink=0.7)
        plt.subplot(3, 1, 3)
        plt.imshow(w_slice, vmin=w_min, vmax=w_max, extent=[0, 12.8, 0, 3.2])
        plt.title(f"w (at z={int(z[z_level])}m, in m/s)")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        plt.colorbar(shrink=0.7)
        plt.tight_layout()
        plt.savefig(f"../frames/xy_{i:04d}.png")
        plt.close()

    with multiprocessing.Pool() as pool:
        pool.starmap(process_frame, [(i, co2_path_xy[i], u_xy[i, z_level], w_xy[i, z_level], z_level, z_xy, uflux, co2_max, u_min, u_max, w_min, w_max) for i in range(len(co2_path_xy))])

    # %%
    # make video with all frames using ffmpeg with 10 fps
    print("\nMaking video...", end="")
    # !ffmpeg -y -r 7 -i ../frames/xy_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xy_{res}_uflux{uflux}.mp4
    os.system(f"ffmpeg -y -r 7 -i ../frames/xy_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xy_{res}_uflux{uflux}.mp4 > /dev/null 2>&1")
    # os.system(f"ffmpeg -y -r 3 -i ../frames/xy_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xy_{res}_uflux{uflux}_slow.mp4")
    print("done")


# %%
u_xy_mean = np.mean(u_xy[100:], axis=0)
v_xy_mean = np.mean(v_xy[100:], axis=0)
w_xy_mean = np.mean(w_xy[100:], axis=0)

u_prime = u_xy - u_xy_mean
v_prime = v_xy - v_xy_mean
w_prime = w_xy - w_xy_mean

tke = 0.5 * (u_prime**2 + v_prime**2 + w_prime**2)
tke_time_series = np.mean(tke, axis=(2, 3))

# %%
plt.figure(figsize=(15, 6))
plt.title("Time series of TKE")
plt.plot(tke_time_series)
plt.xlabel("Time [s]")
plt.xticks(np.arange(0, 340, 50), np.arange(0, 340*300, 50*300))
plt.ylabel("avg'd TKE [m²/s²]")
plt.legend([f"h{i+1}={int(z_xy[i])}m" for i in range(len(z_xy))])
plt.savefig(f"../tke_{res}_uflux{uflux}.png")

# %%
ds_co2_xz = nc.Dataset("nc_files/co2.xz.nc") # type: ignore
ds_u_xz = nc.Dataset("nc_files/u.xz.nc") # type: ignore
ds_v_xz = nc.Dataset("nc_files/v.xz.nc") # type: ignore
ds_w_xz = nc.Dataset("nc_files/w.xz.nc") # type: ignore

# %%
co2_xz = np.array(ds_co2_xz.variables["co2"][:])
u_xz = np.array(ds_u_xz.variables["u"][:])
v_xz = np.array(ds_v_xz.variables["v"][:])
w_xz = np.array(ds_w_xz.variables["w"][:])

y_xz = np.array(ds_co2_xz.variables["y"][:])

# %%
# save co2, u, v, w
np.save("npy_files/co2_xz.npy", co2_xz)
np.save("npy_files/u_xz.npy", u_xz)
np.save("npy_files/v_xz.npy", v_xz)
np.save("npy_files/w_xz.npy", w_xz)

# %%
# load co2, u, v, w from npy files
# co2 = np.load("npy_files/co2_xz.npy")
# u = np.load("npy_files/u_xz.npy")
# v = np.load("npy_files/v_xz.npy")
# w = np.load("npy_files/w_xz.npy")

if not no_frames:
    # %%
    co2_min, co2_max = np.min(co2_xz), np.max(co2_xz)
    u_min, u_max = np.min(u_xz), np.max(u_xz)
    u_range = u_max - u_min
    u_min, u_max = u_min + 0.1*u_range, u_max - 0.1*u_range
    w_min, w_max = np.min(w_xz), np.max(w_xz)
    w_range = w_max - w_min
    w_min, w_max = w_min + 0.1*w_range, w_max - 0.1*w_range

    def process_frame_xz(i, co2_slice, u_slice, w_slice, y_level, y, uflux, co2_max, u_min, u_max, w_min, w_max):
        print(f"\rFrame {i+1:3d}/{len(co2_xz):3d}", end="")
        plt.figure(figsize=(15, 12))
        plt.suptitle(f"Simulation with uflux={uflux} m/s\nt={i*100} s")
        plt.subplot(3, 1, 1)
        plt.imshow(co2_slice, vmin=0, vmax=co2_max*0.8, extent=[0, 12.8, 0, 3.2], origin="lower")
        plt.title("CO2 (integrated)")
        plt.xlabel("x [km]")
        plt.ylabel("z [m]")
        plt.colorbar(shrink=0.7)
        plt.subplot(3, 1, 2)
        plt.imshow(u_slice, vmin=u_min, vmax=u_max, extent=[0, 12.8, 0, 3.2], origin="lower")
        plt.title(f"u (at y={int(y[y_level])}km, in m/s)")
        plt.xlabel("x [km]")
        plt.ylabel("z [m]")
        plt.colorbar(shrink=0.7)
        plt.subplot(3, 1, 3)
        plt.imshow(w_slice, vmin=w_min, vmax=w_max, extent=[0, 12.8, 0, 3.2], origin="lower")
        plt.title(f"w (at y={int(y[y_level])}km, in m/s)")
        plt.xlabel("x [km]")
        plt.ylabel("z [m]")
        plt.colorbar(shrink=0.7)
        plt.tight_layout()
        plt.savefig(f"../frames/xz_{i:04d}.png")
        plt.close()

    y_level = 0
    with multiprocessing.Pool() as pool:
        pool.starmap(process_frame_xz, [(i, co2_xz[i, :, y_level], u_xz[i, :, y_level], w_xz[i, :, y_level], y_level, y_xz, uflux, co2_max, u_min, u_max, w_min, w_max) for i in range(len(co2_xz))])

    # %%
    # make video with all frames using ffmpeg with 10 fps
    print("\nMaking video...", end="")
    # !ffmpeg -y -r 7 -i ../frames/xz_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xz_{res}_uflux{uflux}.mp4
    os.system(f"ffmpeg -y -r 7 -i ../frames/xz_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xz_{res}_uflux{uflux}.mp4 > /dev/null 2>&1")
    # os.system(f"ffmpeg -y -r 3 -i ../frames/xz_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xz_{res}_uflux{uflux}_slow.mp4")
    print("done")

# %%
u_xz_mean = np.mean(u_xz[:], axis=0)
v_xz_mean = np.mean(v_xz[:], axis=0)
w_xz_mean = np.mean(w_xz[:], axis=0)

u_prime = u_xz - u_xz_mean
v_prime = v_xz - v_xz_mean
w_prime = w_xz - w_xz_mean

tke = 0.5 * (u_prime**2 + v_prime**2 + w_prime**2)
tke_time_series = np.mean(tke, axis=(1, 3))

# %%
plt.figure(figsize=(15, 6))
plt.title("Time series of TKE")
plt.plot(tke_time_series)
plt.xlabel("Time [s]")
plt.xticks(np.arange(0, 340, 50), np.arange(0, 340*300, 50*300))
plt.ylabel("avg'd TKE [m²/s²]")
plt.ylim(bottom=0)
plt.savefig(f"../tke_{res}_uflux{uflux}_xz.png")
