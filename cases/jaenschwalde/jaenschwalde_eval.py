# %%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
from IPython.display import clear_output

# set plt font size
plt.rcParams.update({'font.size': 16})

# %%
# res = "256_64_192"
# res = "128_32_96"
res = "64_16_48"
uflux = "2"

os.chdir("/home/georg.trede/MasterThesis/env/microhh/cases/jaenschwalde/snaps_{}_uflux{}".format(res, uflux))

# %%
ds_co2_path = nc.Dataset("nc_files/co2_path.xy.nc")
ds_co2 = nc.Dataset("nc_files/co2.xy.nc")
ds_u = nc.Dataset("nc_files/u.xy.nc")
ds_v = nc.Dataset("nc_files/v.xy.nc")
ds_w = nc.Dataset("nc_files/w.xy.nc")

# %%
co2_path = np.array(ds_co2_path.variables["co2_path"][:])
co2 = np.array(ds_co2.variables["co2"][:])
u = np.array(ds_u.variables["u"][:])
v = np.array(ds_v.variables["v"][:])
w = np.array(ds_w.variables["w"][:])

z = np.array(ds_co2.variables["z"][:])

# %%
# save co2_path, co2, u, v, w
np.save("npy_files/co2_path.npy", co2_path)
for i in range(co2.shape[1]):  # co2.shape[1] should give the number of heights
    np.save(f'npy_files/co2_h{i+1}_xy.npy', co2[:, i, :, :])
    np.save(f'npy_files/u_h{i+1}_xy.npy', u[:, i, :, :])
    np.save(f'npy_files/v_h{i+1}_xy.npy', v[:, i, :, :])
    np.save(f'npy_files/w_h{i+1}_xy.npy', w[:, i, :, :])

# %%
# load co2_path, co2, u, v, w from npy files
co2_path = np.load("npy_files/co2_path.npy")
co2 = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
u = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
v = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))
w = np.zeros((co2_path.shape[0], 3, co2_path.shape[1], co2_path.shape[2]))

for i in range(co2.shape[1]):  # co2.shape[1] should give the number of heights
    co2[:, i, :, :] = np.load(f'npy_files/co2_h{i+1}_xy.npy')
    u[:, i, :, :] = np.load(f'npy_files/u_h{i+1}_xy.npy')
    v[:, i, :, :] = np.load(f'npy_files/v_h{i+1}_xy.npy')
    w[:, i, :, :] = np.load(f'npy_files/w_h{i+1}_xy.npy')

# %%
z_level = 2

co2_min, co2_max = np.min(co2_path), np.max(co2_path)
u_min, u_max = np.min(u[:, z_level]), np.max(u[:, z_level])
u_range = u_max - u_min
u_min, u_max = u_min + 0.1 * u_range, u_max - 0.1 * u_range
w_min, w_max = np.min(w[:, z_level]), np.max(w[:, z_level])
w_range = w_max - w_min
w_min, w_max = w_min + 0.1 * w_range, w_max - 0.1 * w_range

# %%
for i in range(len(co2_path)):
    print(f"\rFrame {i+1:3d}/{len(co2_path):3d}", end="")
    plt.figure(figsize=(15, 12))
    plt.suptitle(f"Simulation with uflux={uflux} m/s\nt={i*300} s")
    plt.subplot(3, 1, 1)
    plt.imshow(co2_path[i], vmin=0, vmax=co2_max*0.8, extent=[0, 12.8, 0, 3.2])
    plt.title("CO2 (integrated)")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    plt.colorbar(shrink=0.7)
    plt.subplot(3, 1, 2)
    plt.imshow(u[i, z_level], vmin=u_min, vmax=u_max, extent=[0, 12.8, 0, 3.2])
    plt.title(f"u (at z={int(z[2])}m, in m/s)")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    plt.colorbar(shrink=0.7)
    plt.subplot(3, 1, 3)
    plt.imshow(w[i, z_level], vmin=w_min, vmax=w_max, extent=[0, 12.8, 0, 3.2])
    plt.title(f"w (at z={int(z[2])}m, in m/s)")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    plt.colorbar(shrink=0.7)
    plt.tight_layout()
    plt.savefig(f"../frames/xy_{i:04d}.png")
    plt.close()
    # plt.pause(0.01)
    # clear_output(wait=True)
print()

# %%
# make video with all frames using ffmpeg with 10 fps
!ffmpeg -y -r 7 -i ../frames/xy_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xy_{res}_uflux{uflux}.mp4
!ffmpeg -y -r 3 -i ../frames/xy_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xy_{res}_uflux{uflux}_slow.mp4

# %%
plt.figure(figsize=(20, 12))
plt.imshow(w[200, z_level], extent=[0, 12.8, 0, 3.2])
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.colorbar(shrink=0.2)
plt.show()

# %%
u_mean = np.mean(u[100:], axis=0)
v_mean = np.mean(v[100:], axis=0)
w_mean = np.mean(w[100:], axis=0)

u_prime = u - u_mean
v_prime = v - v_mean
w_prime = w - w_mean

tke = 0.5 * (u_prime**2 + v_prime**2 + w_prime**2)
tke_time_series = np.mean(tke, axis=(2, 3))

# %%
plt.figure(figsize=(20, 12))
plt.imshow(w_mean[z_level], extent=[0, 12.8, 0, 3.2], vmin=u_min, vmax=u_max)
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.colorbar(shrink=0.2)
plt.show()

# %%
plt.figure(figsize=(15, 6))
plt.title("Time series of TKE")
plt.plot(tke_time_series)
plt.xlabel("Time [s]")
plt.xticks(np.arange(0, 340, 50), np.arange(0, 340*300, 50*300))
plt.ylabel("avg'd TKE [m²/s²]")
plt.legend([f"h{i+1}={int(z[i])}m" for i in range(len(z))])
plt.savefig(f"../tke_{res}_uflux{uflux}.png")

# %%
np.mean(tke_time_series[-100:], axis=0)

# %%
ds_co2 = nc.Dataset("nc_files/co2.xz.nc")
ds_u = nc.Dataset("nc_files/u.xz.nc")
ds_v = nc.Dataset("nc_files/v.xz.nc")
ds_w = nc.Dataset("nc_files/w.xz.nc")

# %%
co2 = np.array(ds_co2.variables["co2"][:])
u = np.array(ds_u.variables["u"][:])
v = np.array(ds_v.variables["v"][:])
w = np.array(ds_w.variables["w"][:])

# %%
# save co2, u, v, w
np.save("npy_files/co2_xz.npy", co2)
np.save("npy_files/u_xz.npy", u)
np.save("npy_files/v_xz.npy", v)
np.save("npy_files/w_xz.npy", w)

# %%
# load co2, u, v, w from npy files
co2 = np.load("npy_files/co2_xz.npy")
u = np.load("npy_files/u_xz.npy")
v = np.load("npy_files/v_xz.npy")
w = np.load("npy_files/w_xz.npy")

# %%
co2_min, co2_max = np.min(co2), np.max(co2)
u_min, u_max = np.min(u), np.max(u)
u_range = u_max - u_min
u_min, u_max = u_min + 0.1*u_range, u_max - 0.1*u_range
w_min, w_max = np.min(w), np.max(w)
w_range = w_max - w_min
w_min, w_max = w_min + 0.1*w_range, w_max - 0.1*w_range

for i in range(len(co2)):
    print(f"\rFrame {i+1:3d}/{len(co2):3d}", end="")
    plt.figure(figsize=(10, 10))
    plt.suptitle(f"Simulation with uflux={uflux} m/s\nt={i*300} s")
    plt.subplot(3, 1, 1)
    plt.imshow(co2[i, :, 0, :], origin="lower", vmin=0, vmax=co2_max*0.8, extent=[0, 12.8, 0, 5.0],
    aspect=0.75)
    plt.title("CO2 (at y=1.6km)")
    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.colorbar(shrink=0.7)
    plt.subplot(3, 1, 2)
    plt.imshow(u[i, :, 0, :], origin="lower", vmin=u_min, vmax=u_max, extent=[0, 12.8, 0, 5.0],
    aspect=0.75)
    plt.title("u (at y=1.6km, in m/s)")
    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.colorbar(shrink=0.7)
    plt.subplot(3, 1, 3)
    plt.imshow(w[i, :, 0, :], origin="lower", vmin=w_min, vmax=w_max, extent=[0, 12.8, 0, 5.0],
    aspect=0.75)
    plt.title("w (at y=1.6km, in m/s)")
    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.colorbar(shrink=0.7)
    plt.tight_layout()
    plt.savefig(f"../frames/xz_{i:04d}.png")
    plt.close()
    # plt.pause(0.01)
    # clear_output(wait=True)
print()

# %%
# make video with all frames using ffmpeg with 10 fps
!ffmpeg -y -r 7 -i ../frames/xz_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xz_{res}_uflux{uflux}.mp4
!ffmpeg -y -r 3 -i ../frames/xz_%04d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p ../xz_{res}_uflux{uflux}_slow.mp4

# %%
u_mean = np.mean(u[:], axis=0)
v_mean = np.mean(v[:], axis=0)
w_mean = np.mean(w[:], axis=0)

u_prime = u - u_mean
v_prime = v - v_mean
w_prime = w - w_mean

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

# %%
plt.figure(figsize=(20, 12))
plt.imshow(u_mean[:,0,:], origin="lower", extent=[0, 12.8, 0, 3.2])#, vmin=u_min, vmax=u_max)
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.colorbar(shrink=0.2)
plt.show()

# %%
np.mean(tke_time_series[-100:], axis=0)

# %% [markdown]
# # Summary of results

# %%
ufluxes = [4, 3, 1, 6, 8]

TKE_h1 = [2.575, 2.595, 2.427, 2.635, 2.783]
TKE_h2 = [2.433, 2.466, 2.304, 2.595, 2.842]
TKE_h3 = [2.500, 2.503, 2.375, 2.610, 2.838]
TKE_xy = [1.861, 1.769, 1.982, 1.870, 1.930]

# sort by ufluxes
TKE_h1 = [TKE_h1[i] for i in np.argsort(ufluxes)]
TKE_h2 = [TKE_h2[i] for i in np.argsort(ufluxes)]
TKE_h3 = [TKE_h3[i] for i in np.argsort(ufluxes)]
TKE_xy = [TKE_xy[i] for i in np.argsort(ufluxes)]
ufluxes = sorted(ufluxes)

plt.plot(ufluxes, TKE_h1, label="h1", marker="o")
plt.plot(ufluxes, TKE_h2, label="h2", marker="o")
plt.plot(ufluxes, TKE_h3, label="h3", marker="o")
plt.show()

plt.plot(ufluxes, TKE_xy, label="xy", marker="o")
plt.show()

# %% [markdown]
# # One cell for execution on server

# %%
import os
import numpy as np
import netCDF4 as nc

# res = "256_64_192"
# res = "128_32_96"
res = "64_16_48"
uflux = "2"

os.chdir("/Users/georg/Studium/Master/Thesis/env/microhh/cases/jaenschwalde/snaps_{}_uflux{}".format(res, uflux))

ds_co2_path = nc.Dataset("nc_files/co2_path.xy.nc")
ds_co2 = nc.Dataset("nc_files/co2.xy.nc")
ds_u = nc.Dataset("nc_files/u.xy.nc")
ds_v = nc.Dataset("nc_files/v.xy.nc")
ds_w = nc.Dataset("nc_files/w.xy.nc")

co2_path = np.array(ds_co2_path.variables["co2_path"][:])
co2 = np.array(ds_co2.variables["co2"][:])
u = np.array(ds_u.variables["u"][:])
v = np.array(ds_v.variables["v"][:])
w = np.array(ds_w.variables["w"][:])

z = np.array(ds_co2.variables["z"][:])

# save co2_path, co2, u, v, w
np.save("npy_files/co2_path.npy", co2_path)
for i in range(co2.shape[1]):  # co2.shape[1] should give the number of heights
np.save(f'npy_files/co2_h{i+1}_xy.npy', co2[:, i, :, :])
np.save(f'npy_files/u_h{i+1}_xy.npy', u[:, i, :, :])
np.save(f'npy_files/v_h{i+1}_xy.npy', v[:, i, :, :])
np.save(f'npy_files/w_h{i+1}_xy.npy', w[:, i, :, :])

ds_co2 = nc.Dataset("nc_files/co2.xz.nc")
ds_u = nc.Dataset("nc_files/u.xz.nc")
ds_v = nc.Dataset("nc_files/v.xz.nc")
ds_w = nc.Dataset("nc_files/w.xz.nc")

co2 = np.array(ds_co2.variables["co2"][:])
u = np.array(ds_u.variables["u"][:])
v = np.array(ds_v.variables["v"][:])
w = np.array(ds_w.variables["w"][:])

# save co2, u, v, w
np.save("npy_files/co2_xz.npy", co2)
np.save("npy_files/u_xz.npy", u)
np.save("npy_files/v_xz.npy", v)
np.save("npy_files/w_xz.npy", w)


