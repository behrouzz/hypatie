import numpy as np
import matplotlib.pyplot as plt

def plot_xyz(x, y, z, color, size):
    """cartesian plot"""
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.2)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='z', labelsize=8)
    ax.scatter(x,y,z, c=color, s=size)
    return ax

def plot_radec(ra, dec, color, size):
    """polar plot"""
    ra  = [(i/180)*np.pi for i in ra] #convert to radians
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8], polar=True)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(90, 0)
    ax.set_yticks(np.arange(0,91,15))
    ax.scatter(ra, dec, c=color, s=size)
    ax.grid(True, alpha=0.7)
    return ax


