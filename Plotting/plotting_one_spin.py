import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def load_traj(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    t = data["t"]
    mx = data["mx"]
    my = data["my"]
    mz = data["mz"]
    return t, mx, my, mz

def plot_components(t, mx, my, mz):
    plt.figure()
    plt.plot(t, mx, label="m_x")
    plt.plot(t, my, label="m_y")
    plt.plot(t, mz, label="m_z")
    plt.plot(t, np.sqrt(mx**2 + my**2 + mz**2), label="|m|", linestyle="--", color="black")
    plt.xlabel("t [s]")
    plt.ylabel("m")
    plt.title("Spin components over time")
    plt.legend()
    plt.grid(True)

def plot_3d_trajectory(t, mx, my, mz, stride=5, cmap="plasma",
                           show_sphere=True, elev=25, azim=45):
    t  = np.asarray(t)[::stride]
    mx = np.asarray(mx)[::stride]
    my = np.asarray(my)[::stride]
    mz = np.asarray(mz)[::stride]

    pts = np.column_stack([mx, my, mz]).reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")

    # color by time
    norm = Normalize(vmin=t.min(), vmax=t.max())
    lc = Line3DCollection(segs, cmap=cmap, norm=norm, linewidth=2.5, alpha=1.0)
    lc.set_array(t[:-1])
    ax.add_collection3d(lc)

    # unit sphere
    if show_sphere:
        u = np.linspace(0, 2*np.pi, 80)
        v = np.linspace(0, np.pi, 40)
        xs = np.outer(np.cos(u), np.sin(v))
        ys = np.outer(np.sin(u), np.sin(v))
        zs = np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(xs, ys, zs, alpha=0.08, linewidth=0, antialiased=True)

    # axis limits
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 1)) 

    # labels and view angle
    ax.set_xlabel("m_x"); ax.set_ylabel("m_y"); ax.set_zlabel("m_z")
    ax.set_title("Spin trajectory in 3D")
    ax.view_init(elev=elev, azim=azim)

    cbar = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("time [s]")

    plt.tight_layout()

if __name__ == "__main__":
    t, mx, my, mz = load_traj("LLG_solver/traj_single_spin.csv")

    plot_components(t, mx, my, mz)

    plot_3d_trajectory(t, mx, my, mz, stride=5)

    plt.show()