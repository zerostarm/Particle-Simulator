'''
Created on May 5, 2020

@author: Stephen
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time as tyme
import gc
import multiprocessing

from Particle import *
from Pendulum_Function_Library import *
from Plotting import *
from Mesh import *

C = 299792458                   #speed of light in vacuum
G = 6.674*10**-11               #Newtonian gravity constant
epsilon0 = 8.854*10**-12        #Epsilon naught
mu0 = 1/(C**2*epsilon0)         #Mu naught
e_C = 1.60217663e-19            #Electron Charge [Coulombs]
nucleon_mass = 1.66054e-27      #Kg per Amu

def generateParticles(number, volume, mass=4, charge=2, manager=None): #Generates Particles in the defined rectangular volume [Length-x, width-y, height-z]
    """
    Generates a *number* of particles in a box defined by *volume*
    returns a list of particle objects
    """
    rand = np.random
    particles = []
    
    for i in range(number): 
        position = [0.7,0,0]
        
        t = 1e-9
        vx, vy, vz = 0,0,0
        try:
            vx = np.sqrt(rand.uniform(0, t, dtype=np.float64))
        except:
            pass
        try:
            vy = np.sqrt(rand.uniform(0, t - vx**2, dtype=np.float64))
        except:
            pass
        try:
            vz = np.sqrt(rand.uniform(0, t - vx**2 - vy**2, dtype=np.float64))
        except:
            pass
        velocity =   [vx, vy, vz] #Three Velocity 
        
        '''
        print(type(i))
        print(type(mass))
        print(type(charge))
        print(type(radius))
        print(type(position))
        print(type(velocity))
        print(type(energy))
        print(type(acceleration))'''
        
        p = Particle(i, mass*nucleon_mass, charge*e_C, position, velocity)
        particles.append(p)
    return particles

def update_particle(particle1, mesh):
    """
    Todo: Add and validate angular momentum as a force
    """
    particles = mesh.particles
    dt = mesh.dt
    volume = mesh.volume_bounds
    bounce_factor = mesh.bounce_factor
    
    Scalar_field = np.zeros(3)
    temp_position = particle1.getPosition()
    E_field, B_field, Scalar_field = mesh.get_fields_at_point_scipy(particles, particle1)
    
    #x,y,z = temp_position
    #R = 10
    #B_field += 100*np.asarray([-y/R, x/R, 0])
        
    #E_field += 100*asarray([0, 0, np.sin(0.1*temp_position[0] - 10*len(particle1.Position)*dt)])
    #B_field += 100*asarray([0, np.sin(0.1*temp_position[0] - 10*len(particle1.Position)*dt), 0])
    
    """
    Boris particle mover
    """
    #print(dt)
    #print(particle1.getCharge())
    #print(particle1.getMass())
    q_prime = dt*particle1.getCharge()/(2*particle1.getMass())
    H = q_prime*B_field
    S = 2*H/(1+np.linalg.norm(H)**2)
    U = particle1.getVelocity() + q_prime*E_field #+ dt/(2*particle1.getMass())*Scalar_field
    U_prime = U + np.cross((U + (np.cross(U, H))), S)
    V = U_prime + q_prime*E_field + dt*Scalar_field
    pos = temp_position
    posn = pos + V*dt
    
    """
    Leapfrogging method that works the most simple way with static B-field drift error
    """"""
    a = F/particle1.getMass()
    V = particle1.getVelocity() + a*dt
    """
    
    xn = posn[0]
    yn = posn[1]
    zn = posn[2]
    
    
    vx = V[0] 
    vy = V[1]
    vz = V[2]
    """
    Wall Collisions
        Sets the position to the wall position and reverses the velocity away from the wall
        Bounce Factor allows not perfect energy conserving collisions with the wall
        v= v-intowall *bounce factor 
        position = wall*sign(x) 
    """    
    #xp, yp, zp = xn, yn, zn
    xf, yf, zf = xn, yn, zn
    if np.abs(xn) > volume[0]:
        vx = -vx * bounce_factor
        xf = np.sign(xn)*volume[0]
    else:
        pass
    
    if np.abs(yn) > volume[1]:
        vy = -vy * bounce_factor
        yf = np.sign(yn)*volume[0]
    else:
        pass
    
    if np.abs(zn) > volume[2]:
        vz = -vz * bounce_factor
        zf = np.sign(zn)*volume[0]
    else:
        pass
    
    pfinal = np.asarray([xf,yf,zf])
    
    #particle1.setNewAcceleration(a) #Uncomment if using Leapfrogging methods
    particle1.setNewVelocity([vx, vy, vz])
    particle1.setNewPosition(pfinal)
    return True

def run_simulation(particles, volume_bounds, dt, mesh_fineness = 1.0):
    print("Starting Simulation")
    global mesh
    
    print("Making Mesh")
    mesh = Mesh(mesh_fineness, volume_bounds, dt, particles,c=C, epsilon0=epsilon0)
    print("Mesh Done")
    
    print("Starting For Loop")
    times = []
    total_start_time = tyme.time()
    for i in range(t_step_number): #range(0,1):
        print("Timestep:", i)
        start_time = tyme.time()
        for particle in particles:
            update_particle(particle, mesh)
        
        for particle in particles:
            particle.setAcceleration()
            particle.setVelocity()
            particle.setPosition()
            #particle.setKe()
            #particle.setPe()
            #particle.setTe()
        times.append(tyme.time()-start_time)
    total_end_time = tyme.time()
    
    print("Finished For Loop")
    
    """
    Raw Printing out the particle start and finish positions
    Also printing the standard particle info - mass, charge, radius, energy 
    """
    print("Particle Stats")
    for i in range(len(particles)):
        print(particles[i].Position[0], particles[i].Position[-1])
        print(particles[i])
    
    print("Simulation Done")
    return particles, times, total_start_time, total_end_time

    
if __name__ == "__main__": 
    """
    Defining some paths that are used for debug stuff
    """
    path = r'C:\Users\Stephen\Desktop\Eclipse\Workspace\Particle-Simulator\pictures\\'
    path_2d = r'C:\Users\Stephen\Desktop\Eclipse\Workspace\Particle-Simulator\pictures_2d\\'
    
    """
    Defining simulation constants
    """
    volume_bounds = [10, 10, 10]                                          #Assumes a box centered on [0,0,0] with walls positioned at [+-10,0,0], [0,+-10,0], and [0,0,+-10] 
    number_of_particles = 1                                             #Number of particles
    number_xyz_steps = 100
    dx = volume_bounds[0]/number_xyz_steps
    
    time_bounds = [0,10]                                                       #Time bounds. These aren't necessarily needed if you define dt but it's useful for plotting stuff
    t_step_number = 1000                                                #Number of time steps. Always gotta define this or else things will break.
    
    dt = (time_bounds[-1])/t_step_number                                       #The delta of time. I usually define it as (time[-1] - time[0])/t_step_number 
    
    bounce_factor = 0.1                                                 #This is how bouncy the walls are. Can be any number but physically real values are between 0 and 1. With 0 being particles instantly stop at walls and 1 being the are perfectly reflected with no energy loss.
    
    print("Generating Particles")
    particles = generateParticles(number_of_particles, volume_bounds, mass=4, charge=2)   #This generates particles, at random positions in the volume box, and low random velocities. The charge alternates between +1, and -1
    #print(particles)
    
    print(dx)
    particles, times, total_start_time, total_end_time = run_simulation(particles, volume_bounds, dt, mesh_fineness=dx) 
    
    #large_plot(time, t_step_number, particles)
    
    
    
    n=1000
    theta = np.linspace(-np.pi, np.pi, n)
    time = np.linspace(time_bounds[0], time_bounds[1], t_step_number)
    c, a = 0, 1.3
    Toroid_x    = (c + a*np.cos(theta))
    toroid_y    = a * np.sin(theta)
    
    colors = ["b"]
    r_thetas = np.zeros((len(particles), 2))
    for particle in particles:
        color = colors[particle.Index%len(colors)]
        pos = np.asarray(particle.Position)
        r_thetas[particle.Index] = np.asarray([np.sqrt(pos[-1,0]**2 + pos[-1,1]**2), pos[-1,2]])
    plt.scatter(r_thetas[:,0], r_thetas[:,1], color=color, label="Ions")
    plt.plot(Toroid_x, toroid_y, "r-", label="Toroid")
    plt.legend()
    plt.xlim(0)
    #ply.ylim(0, )
    plt.title("He2+ Poloidal Distribution")
    plt.xlabel("Major Radius [m]")
    plt.ylabel("Vertical Position [m]")
    plt.show()
    
    
    r_thetas = []
    for particle in particles[:10]:
        color = colors[particle.Index%len(colors)]
        pos = np.asarray(particle.Position)
        plt.plot(np.sqrt(pos[:,0]**2 + pos[:,1]**2), pos[:,2], label=str(particle.Index))
    plt.plot(Toroid_x, toroid_y, "r-", label="Toroid")
    plt.legend()
    plt.xlim(0)
    #ply.ylim(0, )
    plt.title("He2+ Poloidal Paths")
    plt.xlabel("Major Radius [m]")
    plt.ylabel("Vertical Position [m]")
    plt.show()
    
    for particle in particles[:10]:
        color = colors[particle.Index%len(colors)]
        pos = np.asarray(particle.Position)
        plt.plot(time, np.sqrt(pos[:,0]**2 + pos[:,1]**2)[:-1], label=str(particle.Index))
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("R Position [m]")
    plt.show()
    
    
    #3D Path Plots
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #ax.axes.set_xlim3d(-volume_bounds[0]-.01, volume_bounds[0]+.01)
    #ax.axes.set_ylim3d(-volume_bounds[1]-.01, volume_bounds[1]+.01)
    #ax.axes.set_zlim3d(-volume_bounds[2]-.01, volume_bounds[2]+.01)
    ax.axes.set_xlabel("X")
    ax.axes.set_ylabel("Y")
    ax.axes.set_zlabel("Z")
    for particle in particles:
        #plt.plot(t, particle.Energy, label = str(particle.Index))
        vel = np.asarray(particle.Velocity)
        pos = np.asarray(particle.Position)
        #ax.plot(vel[:,0], vel[:,1], vel[:,2], label = "Vel" + str(particle.Index))
        ax.plot(pos[:,0], pos[:,1], pos[:,2], label = "Pos" + str(particle.Index))
    #plt.legend()
    plt.show()
    
    
    print("Average Execution time:", np.average(times))
    print("Standard Deviation:", np.std(times))
    print("Total Execution time:", total_end_time-total_start_time)
    '''
    #Execution time plots
    plt.hist(times, bins="auto", color="r")
    #plt.hist(times_mesh, bins="auto", color="b")
    ylims = plt.ylim()
    plt.vlines(np.average(times), ylims[0], ylims[1], colors="k", label="Average"); ylims = plt.ylim();
    plt.vlines(np.average(times) - np.std(times), ylims[0], ylims[1], colors="k", linestyle="--"); ylims = plt.ylim(); 
    plt.vlines(np.average(times) + np.std(times), ylims[0], ylims[1], colors="k", linestyle="--")
    """
    #Mesh
    plt.vlines(np.average(times_mesh), ylims[0], ylims[1], colors="k", label="Average_mesh"); ylims = plt.ylim();
    plt.vlines(np.average(times_mesh) - np.std(times_mesh), ylims[0], ylims[1], colors="k", linestyle="--"); ylims = plt.ylim(); 
    plt.vlines(np.average(times_mesh) + np.std(times_mesh), ylims[0], ylims[1], colors="k", linestyle="--")
    """
    plt.xlabel("Times")
    plt.ylabel("Number")
    plt.legend()
    plt.show()
    '''