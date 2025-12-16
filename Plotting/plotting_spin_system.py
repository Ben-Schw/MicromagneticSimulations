import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def load_traj(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    t = data["t"]
    mx = []
    my = []
    mz = []
    for i in range(1, len(data.dtype.names), 3):
        mx.append(data[data.dtype.names[i]])
        my.append(data[data.dtype.names[i+1]])
        mz.append(data[data.dtype.names[i+2]])
    
    return t, np.array(mx), np.array(my), np.array(mz)