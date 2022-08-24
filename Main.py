'''
Created on May 5, 2020

@author: Stephen
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time as tyme
import gc

from Particle import *
from Pendulum_Function_Library import *
from Plotting import *
from Mesh import *

C = 299792458           #speed of light in vacuum
G = .1                  #6.674*10**-11 #Newtonian gravity constant
epsilon0 = .001         #8.854*10**-12 #Epsilon naught
mu0 = 1/(C**2*epsilon0) #Mu naught

def generateParticles(number, volume): #Generates Particles in the defined rectangular volume [Length-x, width-y, height-z]
    """
    Generates a *number* of particles in a box defined by *volume*
    returns a list of particle objects
    """
    rand = np.random
    particles = []
    for i in range(number):
        mass = 20#rand.choice(range(1,2), 1)#integer
        charge =  [1,-1][i%2] #rand.choice([-1,0,1], 1)[0]#+-1, 0 integer 
        position = [asarray([1,0,0]), asarray([-1,0,0]),
                    asarray([0,1,0]), asarray([0,-1,0]), 
                    asarray([0,0,1]), asarray([0,0,-1])][i%6]
        """[[10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)],
                    [-10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)]][i%2]""" 
                #[rand.randint(-volume[0]+1, volume[0]-1), rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)] #Three Position
                # [asarray([1,0,0]), asarray([-1,0,0]), asarray([0,0,1]), asarray([0,0,-1])][i%3]
        #print(position)
        radius = 0.01 #mass/2 #float
        
        t = 2.0
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
        velocity =  [asarray([0,0,0]), asarray([0,0,0]), asarray([0,0,0])][i%3] #[vx, vy, vz] #Three Velocity [asarray([0,1,0]), asarray([0,-1,0]), asarray([1,0,0])][i%3]
        
        linearE = 1/2*mass*np.linalg.norm(velocity)**2
        Ke = linearE    #mass*C**2 +                                #For now just KE, so mass*||v||**2
        acceleration = [0.0, 0.0, 0.0]                                  #Three Acceleration
        
        '''
        print(type(i))
        print(type(mass))
        print(type(charge))
        print(type(radius))
        print(type(position))
        print(type(velocity))
        print(type(energy))
        print(type(acceleration))'''
        
        p = Particle(i, mass, charge, radius, position, velocity, Ke, acceleration)
        particles.append(p)
        
    for particle1 in particles:
        Pe = 0
        for particle2 in particles:
            if particle1.Index != particle2.Index:
                r = particle2.getPosition() - particle1.getPosition()
                r_norm = np.linalg.norm(r, 2)
                Pe += 1/(4*np.pi*epsilon0)*particle2.getCharge()/r_norm #E-field
                Pe += -G*particle1.getMass()*particle2.getMass()/r_norm #Gravity Potential
            else:
                pass
        particle1.setPe_directly(Pe, 0)
    return particles
def gamma(v):
    return 1/np.sqrt(1- v*v/C**2)

def get_forces(Scalar_field, B_field, E_field, bounce_factor, particle1, particle2, temp_position):
    Scalar_additions = np.zeros(3)
    B_additions = np.zeros(3)
    E_additions = np.zeros(3)
    pe = 0
    
    if particle1.Index == particle2.Index:
            pass
    else:
        r = particle2.getPosition() - particle1.getPosition()
        r *=-1
        r_norm = np.linalg.norm(r, 2)
          
        r_min = particle1.getRadius() + particle2.getRadius()
        temp_position = particle1.getPosition()
        if r_norm < r_min:
            #F_additions *= 0#-bounce_factor
            #Scalar_additions *= 0#-bounce_factor
            #B_additions *= 0#-bounce_factor
            pass
        else:
            r_hat = r/r_norm
            
            Scalar_additions += G*particle1.getMass()*particle2.getMass()/r_norm**2 * -r #Gravity
            E_additions += 1/(4*np.pi*epsilon0)*particle1.getCharge()*particle2.getCharge()/r_norm**2 * r #Coulomb
            
            B_additions += mu0/(4*np.pi)*particle2.getCharge()*np.cross(particle2.getVelocity(), r_hat)/r_norm**3 #B field from moving charged particle - Biot-Savart law implementation for a discrete particle
            #F_additions += particle1.getCharge()*np.cross(particle1.getVelocity(), B_vec)
            
            pe += 1/(4*np.pi*epsilon0)*particle2.getCharge()/r_norm #E-fiel Potential
            pe += -g*particle1.getMass()*particle2.getMass()/r_norm #Gravity Potential
    Scalar_field += Scalar_additions
    B_field += B_additions
    E_field += E_additions
    
    return Scalar_field, B_field, E_field, temp_position, pe

def update_particle(particle1, particles, mesh, dt, volume, bounce_factor=0.9, MESH=False):
    """
    Todo: Add and validate angular momentum as a force
    """
    Scalar_field = np.zeros(3)
    B_field = np.zeros(3)
    E_field = np.zeros(3)
    Pe = 0
    temp_position = particle1.getPosition()
    if MESH:
        E_field, B_field = mesh.get_fields_at_point(temp_position[0], temp_position[1], temp_position[2])
    else:
        for particle2 in particles:
            Scalar_field, B_field, E_field, temp_position, pe = get_forces(Scalar_field, B_field, E_field, bounce_factor, particle1, particle2, temp_position)
            Pe += pe
    """
    Fields that are not dependent on other particles must go at this level or in single particle test cases they won't be applied
    """         
    #B_field += np.asarray([0,0,10]) #Static B-Field Check
    #E_field += np.asarray([0,10,0]) #Static E-field Check
    #Scalar_field += -40/np.linalg.norm(particle1.getPosition())**3 *particle1.getPosition()
    #Scalar_field += np.asarray([0,10,0]) #Non E&M field check
    #Scalar_field += -0.5*particle1.getVelocity()**2
    
    #E_field += 100*asarray([0, 0, np.sin(0.1*temp_position[0] - 10*len(particle1.Position)*dt)])
    B_field += 100*asarray([0, np.sin(0.1*temp_position[0] - 10*len(particle1.Position)*dt), 0])
    
    #Pe += -40/np.linalg.norm(particle1.getPosition())
    
    """
    Boris particle mover
    """
    q_prime = dt*particle1.getCharge()/(2*particle1.getMass())
    H = q_prime*B_field
    S = 2*H/(1+np.linalg.norm(H)**2)
    U = particle1.getVelocity() + q_prime*E_field + dt/(2*particle1.getMass())*Scalar_field
    U_prime = U + np.cross((U + (np.cross(U, H))), S)
    V = U_prime + q_prime*E_field + dt/(2*particle1.getMass())*Scalar_field 
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
    
    Ke = 1/2*particle1.getMass()*np.linalg.norm(particle1.getVelocity())**2
    particle1.setNewKe(Ke) #particle1.getMass()*C**2 +
    particle1.setNewPe(Pe)
    particle1.setNewTe(Ke + Pe)
    return True

def update_mesh(mesh, particles):
    mesh.get_new_fields(particles)
    return True

def run_simulation(number_of_particles, volume_bounds, dt, C, epsilon0, mesh_fineness = 5.0, MESH=False):
    particles = generateParticles(number_of_particles, volume_bounds)   #This generates particles, at random positions in the volume box, and low random velocities. The charge alternates between +1, and -1
    
    mesh = Mesh(mesh_fineness, volume_bounds, dt, c=C, epsilon0=epsilon0)
    
    
    #particles_for_multi = []
    #for particle in particles:
    #        particles_for_multi.append([particle, particles, dt, volume_bounds])
    
    times = []
    total_start_time = tyme.time()
    for i in range(t_step_number):
        print(i)
        #MultiprocessingStandard(update_particle, particles_for_multi, prints=False)
        start_time = tyme.time()
        if MESH:
            update_mesh(mesh, particles)
        for particle in particles:
            update_particle(particle, particles, mesh, dt, volume_bounds, MESH=MESH)
        for particle in particles:
            particle.setAcceleration()
            particle.setVelocity()
            particle.setPosition()
            particle.setKe()
            particle.setPe()
            particle.setTe()
        times.append(tyme.time()-start_time)
        #gc.collect()
    total_end_time = tyme.time()
    """
    Raw Printing out the particle start and finish positions
    Also printing the standard particle info - mass, charge, radius, energy 
    """
    for i in range(len(particles)):
        print(particles[i].Position[0], particles[i].Position[-1])
        print(particles[i])
    
    values = []
    values_2d = []
    for i in range(t_step_number):
        values.append((i, particles, path, volume_bounds))
        values_2d.append((i, particles, path_2d, volume_bounds))
    
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
    volume_bounds = [10,10,10]                                          #Assumes a box centered on [0,0,0] with walls positioned at [+-10,0,0], [0,+-10,0], and [0,0,+-10] 
    number_of_particles = 1                                             #Number of particles
    
    time = [0,100]                                                       #Time bounds. These aren't necessarily needed if you define dt but it's useful for plotting stuff
    t_step_number = 10000                                                #Number of time steps. Always gotta define this or else things will break.
    
    dt = (time[-1])/t_step_number                                       #The delta of time. I usually define it as (time[-1] - time[0])/t_step_number 
    
    bounce_factor = .5                                                  #This is how bouncy the walls are. Can be any number but physically real values are between 0 and 1. With 0 being particles instantly stop at walls and 1 being the are perfectly reflected with no energy loss.
    
    particles, times, total_start_time, total_end_time = run_simulation(number_of_particles, volume_bounds, dt, C, epsilon0) 
    
    #particles_mesh, times_mesh, total_start_time_mesh, total_end_time_mesh = run_simulation(number_of_particles, volume_bounds, dt, C, epsilon0, mesh_fineness=5.0, MESH=True)
    
    #MultiprocessingStandard(plot_stuff, values)
    #createVideo(path, "Particles_in_a_box", time, t_step_number)
    #MultiprocessingStandard(plot_stuff_2D, values_2d)
    #createVideo(path_2d, "Particles_in_a_box_2d", time, t_step_number)
    
    large_plot(time, t_step_number, particles)
     
    '''   
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.axes.set_xlim3d(-volume_bounds[0]-.01, volume_bounds[0]+.01)
    ax.axes.set_ylim3d(-volume_bounds[1]-.01, volume_bounds[1]+.01)
    ax.axes.set_zlim3d(-volume_bounds[2]-.01, volume_bounds[2]+.01)
    ax.axes.set_xlabel("X")
    ax.axes.set_ylabel("Y")
    ax.axes.set_zlabel("Z")
    for particle in particles:
        #plt.plot(t, particle.Energy, label = str(particle.Index))
        vel = np.asarray(particle.Velocity)
        pos = np.asarray(particle.Position)
        #ax.plot(vel[:,0], vel[:,1], vel[:,2], label = "Vel" + str(particle.Index))
        ax.plot(pos[:,0], pos[:,1], pos[:,2], label = "Pos" + str(particle.Index))
    plt.legend()
    plt.show()
    '''
    
    print("Average Execution time:", np.average(times))
    print("Standard Deviation:", np.std(times))
    print("Total Execution time:", total_end_time-total_start_time)
    
    plt.hist(times, bins="auto", color="r")
    #plt.hist(times_mesh, bins="auto", color="b")
    ylims = plt.ylim()
    plt.vlines(np.average(times), ylims[0], ylims[1], colors="k", label="Average"); ylims = plt.ylim();
    plt.vlines(np.average(times) - np.std(times), ylims[0], ylims[1], colors="k", linestyle="--"); ylims = plt.ylim(); 
    plt.vlines(np.average(times) + np.std(times), ylims[0], ylims[1], colors="k", linestyle="--")
    '''
    #Mesh
    plt.vlines(np.average(times_mesh), ylims[0], ylims[1], colors="k", label="Average_mesh"); ylims = plt.ylim();
    plt.vlines(np.average(times_mesh) - np.std(times_mesh), ylims[0], ylims[1], colors="k", linestyle="--"); ylims = plt.ylim(); 
    plt.vlines(np.average(times_mesh) + np.std(times_mesh), ylims[0], ylims[1], colors="k", linestyle="--")
    '''
    plt.xlabel("Times")
    plt.ylabel("Number")
    plt.legend()
    plt.show()
    
    
    