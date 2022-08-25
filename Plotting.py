'''
Created on Jul 5, 2022

@author: Stephen
'''
import matplotlib.pyplot as plt
import numpy as np
import os

def large_plot(time, t_step_number, particles):
    fig, ax = plt.subplots(2, 3)#+len(particles))
    ax = ax.flatten()
    t = np.linspace(time[0], time[1], len(particles[0].KineticEnergy))#t_step_number+1)
    fig.tight_layout()
    particle = particles[0]
    
    ke = np.zeros_like(particle.KineticEnergy)
    pe = np.zeros_like(particle.PotentialEnergy)
    
    colors = ["r", "b"]#plt.cm.rainbow(np.linspace(0,1, len(particles)))
    
    for particle in particles:
        vel = np.asarray(particle.Velocity)
        pos = np.asarray(particle.Position)
        color = colors[particle.Index%len(colors)]
        ax[0].plot(pos[:,0], pos[:,1], color=color, label = "Pos" + str(particle.Index))
        ax[0].plot(pos[0,0], pos[0,1], "o", color=color)#, label = "Pos" + str(particle.Index))
        ax[0].plot(pos[-1,0], pos[-1,1], "*", color=color)
        
        ax[1].plot(vel[:,0], vel[:,1], color=color, label = "Vel" + str(particle.Index))
        ax[2].plot(pos[:,0], pos[:,2], color=color, label = "Pos" + str(particle.Index))
        ke += particle.KineticEnergy
        ax[3].plot(t, particle.KineticEnergy, color=color, label = str(particle.Index))
        
        pe += particle.PotentialEnergy
        ax[4].plot(t, particle.PotentialEnergy, color=color, label = str(particle.Index))
        
        ax[5].plot(t, particle.TotalEnergy, color=color, label = str(particle.Index))
    ax[3].plot(t, ke, label = "Total")
    ax[4].plot(t, pe, label = "Total")
    
    totalE = ke + pe
    ax[5].plot(t, totalE, label = "Total")
    
    '''i = particle.Index
        ax[3+i].plot(t, vel[:,0], label = str(particle.Index) + "x")
        ax[3+i].plot(t, vel[:,1], label = str(particle.Index) + "y")
        ax[3+i].plot(t, vel[:,2], label = str(particle.Index) + "z")
        ax[3+i].set_xlabel("Time")
        ax[3+i].set_ylabel("Velocity Component")
        ax[3+i].legend()'''
        
    ax[0].set_title("Position")
    ax[0].set_xlabel("X component")
    ax[0].set_ylabel("Y component")
    ax[0].legend()
    
    ax[2].set_title("Position")
    ax[2].set_xlabel("X component")
    ax[2].set_ylabel("Z component")
    ax[2].legend()
    
    ax[1].set_title("Velocity")
    ax[1].set_xlabel("X component")
    ax[1].set_ylabel("Y component")
    ax[1].legend()
    
    ax[3].set_title("Kinetic Energy")
    ax[3].set_xlabel("Time")
    ax[3].set_ylabel("Energy")
    ax[3].legend()
    
    ax[4].set_title("Potential Energy")
    ax[4].set_xlabel("Time")
    ax[4].set_ylabel("Energy")
    ax[4].legend()
    
    ax[5].set_title("Total Energy")
    ax[5].set_xlabel("Time")
    ax[5].set_ylabel("Energy")
    ax[5].legend()
    
    #plt.legend()
    plt.show()
    
    return True
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