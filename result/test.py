import numpy as np
import matplotlib.pyplot as plt
import os

snap_folder = "snapshots"
timesteps = [0, 200, 400, 600]

def read_minmax(path="stats_minmax.txt"):
    with open(path) as f:
        lines = f.readlines()
    vmin = float(lines[0].split()[1])
    vmax = float(lines[1].split()[1])
    return vmin, vmax

vmin, vmax = read_minmax()
print("Colormap scale from stats_minmax.txt:", vmin, vmax)

def load_snapshot(timestep):
    filename = os.path.join(snap_folder, f"snapshot_{timestep:06d}.dat")
    print(f"Loading {filename}...")
    return np.loadtxt(filename)


data0 = load_snapshot(timesteps[0])

x = data0[:,0]
y = data0[:,1]
z = data0[:,2]

x_vals = np.unique(x)
y_vals = np.unique(y)
z_vals = np.unique(z)

nx = len(x_vals)
ny = len(y_vals)
nz = len(z_vals)
print("Grid:", nx, ny, nz)

z0 = z_vals[nz // 2]
eps = 1e-6

n = len(timesteps)
cols = 2
rows = (n + 1) // 2

plt.figure(figsize=(12, 6))

for k, ts in enumerate(timesteps):
    data = load_snapshot(ts)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    p = data[:,3]

    mask = np.abs(z - z0) < eps
    xs = x[mask]
    ys = y[mask]
    ps = p[mask]

    order = np.lexsort((xs, ys))
    xs = xs[order]
    ys = ys[order]
    ps = ps[order]

    P = ps.reshape(ny, nx)

    plt.subplot(rows, cols, k+1)
    plt.imshow(
        P,
        extent=[x_vals.min(), x_vals.max(),
                y_vals.min(), y_vals.max()],
        origin="lower",
        aspect="equal",
        vmin=vmin,
        vmax=vmax,
        cmap="seismic"
    )
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title(f"Timestep {ts}\n(z = {z0:.2f})")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

plt.tight_layout()
plt.show()