'''
Created on Jul 5, 2022

@author: Stephen
'''
import matplotlib.pyplot as plt
import os

def plot_stuff(i, particles, path, vol_bounds):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.axes.set_xlim3d(-vol_bounds[0]-.01, vol_bounds[0]+.01)
    ax.axes.set_ylim3d(-vol_bounds[1]-.01, vol_bounds[1]+.01)
    ax.axes.set_zlim3d(-vol_bounds[2]-.01, vol_bounds[2]+.01)
    ax.axes.set_xlabel("X")
    ax.axes.set_ylabel("Y")
    ax.axes.set_zlabel("Z")
    
    for particle in particles:
        x,y,z = particle.Position[i]
        ax.scatter(x,y,z, s = particle.Radius)
    filenumber = i
    filenumber = format(filenumber, "05")
    filename = "image{}.png".format(filenumber)
    fig.savefig(os.path.join(path,filename))
    plt.close(fig)
    return True

def plot_stuff_2D(i, particles, path, vol_bounds):
    fig, ax = plt.subplots(1,1)
    #ax = fig.add_subplot()
    ax.set_xlim(-vol_bounds[0]-.01, vol_bounds[0]+.01)
    ax.set_ylim(-vol_bounds[1]-.01, vol_bounds[1]+.01)
    #ax.axes.set_zlim3d(-vol_bounds[2]-.01, vol_bounds[2]+.01)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    #ax.axes.set_zlabel("Z")
    
    for particle in particles:
        x,y,_ = particle.Position[i]
        ax.scatter(x,y, s = particle.Mass/2)
    filenumber = i
    filenumber = format(filenumber, "05")
    filename = "image_2d{}.png".format(filenumber)
    fig.savefig(os.path.join(path,filename))
    plt.close(fig)
    return True