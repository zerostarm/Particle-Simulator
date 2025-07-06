'''
Created on Feb 27, 2024

@author: Stephen
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time as tyme
#import multiprocess
from multiprocess import Manager, Pool, cpu_count, set_start_method

g = 9.81 #m/s^2
C = 299792458           #speed of light in vacuum
G = .1                  #6.674*10**-11 #Newtonian gravity constant
epsilon0 = .001         #8.854*10**-12 #Epsilon naught
mu0 = 1/(C**2*epsilon0) #Mu naught

def generateParticles(number, volume):
    """
    Generates a *number* of particles in a box defined by *volume*
    returns list for each of the required vectors
    """
    rand = np.random
    
    mass_arr = np.zeros(number)
    charge_arr = np.zeros(number)
    radius_arr = np.zeros(number)
    
    position_arr = np.zeros((number, 3))
    velocity_arr = np.zeros((number, 3))
    acceleration_arr = np.zeros((number, 3))
    
    Ke_arr = np.zeros(number)
    Pe_arr = np.zeros(number)
    Te_arr = np.zeros(number)
    
    
    for i in range(number):
        mass = 20#rand.choice(range(1,2), 1)#integer
        charge =  [1,-1][i%2] #rand.choice([-1,0,1], 1)[0]#+-1, 0 integer 
        position = [[10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)],
                    [-10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)]][i%2] 
                #[rand.randint(-volume[0]+1, volume[0]-1), rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)] #Three Position
                # [asarray([1,0,0]), asarray([-1,0,0]), asarray([0,0,1])][i%3]
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
        velocity =  [np.asarray([0,0,0]), np.asarray([0,0,0]), np.asarray([0,0,0])][i%3] #[vx, vy, vz] #Three Velocity [asarray([0,1,0]), asarray([0,-1,0]), asarray([1,0,0])][i%3]
        
        linearE = 1/2*mass*np.linalg.norm(velocity)**2
        Ke = linearE    #mass*C**2 +                                #For now just KE, so mass*||v||**2
        acceleration = [0.0, 0.0, 0.0] 
        
        
        mass_arr[i] = mass
        charge_arr[i] = charge
        radius_arr[i] = radius
        
        position_arr[i] = np.asarray(position)
        velocity_arr[i] = np.asarray(velocity)
        acceleration_arr[i] = np.asarray(acceleration)
        
        Ke_arr[i] = Ke
    
    for i in range(number):
        Pe = 0
        for j in range(number):
            if i != j:
                r = position_arr[j] - position_arr[i]
                r_norm = np.linalg.norm(r, 2)
                Pe += 1/(4*np.pi*epsilon0)*charge_arr[j]/r_norm #E-field
                Pe += -G*mass_arr[i]*mass_arr[j]/r_norm #Gravity Potential
            else:
                pass
            
        Pe_arr[i] = Pe
    Te_arr = Ke_arr + Pe_arr
    
    return mass_arr, charge_arr, radius_arr, position_arr, velocity_arr, acceleration_arr, Ke_arr, Pe_arr, Te_arr 

def update_particle(
        index,
        number_of_particles,
        dt, 
        volume, 
        
        mass_arr, 
        charge_arr, 
        radius_arr, 
        position_arr, 
        velocity_arr, 
        acceleration_arr, 
        Ke_arr, 
        Pe_arr, 
        Te_arr,
        
        r_mass_arr, 
        r_charge_arr, 
        r_radius_arr, 
        r_position_arr, 
        r_velocity_arr, 
        r_acceleration_arr, 
        r_Ke_arr, 
        r_Pe_arr,
        r_Te_arr,
        
        bounce_factor=0.9,
        ):
    '''
    All of these should be multiprocessing.manager objects 
    EXCEPT index and number_of_particles which should be ints. 
    '''
    print("index:", index, "inital ma:", mass_arr)
    print("index:", index, "inital r_ma:", r_mass_arr)
    
    
    particle1_mass = mass_arr[index]
    particle1_charge = charge_arr[index]
    particle1_radius = radius_arr[index]
    
    particle1_current_pos = position_arr[index]
    particle1_current_velocity = velocity_arr[index]
    
    
    F = np.zeros(3)
    Scalar_field = np.zeros_like(F)
    B_field = np.zeros_like(F)
    E_field = np.zeros_like(F)
    Pe = 0
    for j in range(number_of_particles):
        if index == j:
            pass
        else:
            particle2_mass = mass_arr[j] 
            particle2_charge = charge_arr[j]
            particle2_radius = radius_arr[j]
            
            particle2_current_pos = position_arr[j]
            
            r = particle1_current_pos - particle2_current_pos #vector pointing to the main particle
            r_norm = np.linalg.norm(r, 2)
            r_min = particle1_radius + particle2_radius
            
            """
            if r_norm < r_min:
                '''
                Something was supposed to go here to make the particles have different collisions 
                like hard sphere or something but it never happened so idk
                if i will actually add this functionality
                '''
                pass
            else:
                pass
            """
            
            r_hat = r/r_norm
            Scalar_field += G*particle1_mass*particle2_mass/r_norm**2 * -r #Gravity
            E_field += 1/(4*np.pi*epsilon0)*particle1_charge*particle2_charge/r_norm**2 * r #Coulomb
            
            
            '''
            There was some B-field stuff here but am ignoring for now.
            just gonna write the wrapper code for it around this function
            '''
            
            Pe += -g*particle1_mass*particle2_mass/r_norm #Gravity Potential
            Pe += 1/(4*np.pi*epsilon0)*particle2_charge/r_norm #E-fiel Potential
            Pe *= -1 * np.sign( particle1_charge*particle2_charge) ##Pretty sure I forgor the negative sign in the potential energy expression and effed up the signs
    """
    Boris particle mover
    """
    q_prime = dt*particle1_charge/(2*particle1_mass)
    H = q_prime*B_field
    S = 2*H/(1+np.linalg.norm(H)**2)
    U = particle1_current_velocity + q_prime*E_field + dt/(2*particle1_mass)*Scalar_field
    U_prime = U + np.cross((U + (np.cross(U, H))), S)
    V = U_prime + q_prime*E_field + dt/(2*particle1_mass)*Scalar_field 
    
    posn = particle1_current_pos + V*dt
    
    """
    Boundary Conditions Stuff
    """
    xn = posn[0]
    yn = posn[1]
    zn = posn[2]
    
    
    vx = V[0] 
    vy = V[1]
    vz = V[2]
    
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
    velfinal = np.asarray([vx, vy, vz])
    
    Ke = 1/2*particle1_mass*np.linalg.norm(velfinal)**2
    
    
    mass_arr[index] = particle1_mass
    charge_arr[index] = particle1_charge
    radius_arr[index] = particle1_radius
    #position_arr[index] = pfinal 
    #velocity_arr[index] = velfinal
    #acceleration_arr[index] = acceleration_arr[index]
    #Ke_arr[index] = Ke
    #Pe_arr[index] = Pe
    #Te_arr[index] = Ke+Pe
    
    r_mass_arr = mass_arr
    r_charge_arr = charge_arr
    r_radius_arr = radius_arr
    
    r_position_arr[index] = pfinal#position_arr 
    r_velocity_arr[index] = velfinal#velocity_arr
    #r_acceleration_arr = acceleration_arr
    r_Ke_arr = Ke_arr
    r_Pe_arr = Pe_arr
    r_Te_arr = Te_arr
    
    print("index:", index, "final ma:", mass_arr)
    print("index:", index, "final r_ma:", r_mass_arr)
    
    return r_mass_arr, r_charge_arr, r_radius_arr, r_position_arr, r_velocity_arr, r_acceleration_arr, r_Ke_arr, r_Pe_arr,r_Te_arr
    #return True

if __name__ == '__main__':
    
    """
    Defining simulation constants
    """
    volume_bounds = [100,100,100]                                          #Assumes a box centered on [0,0,0] with walls positioned at [+-10,0,0], [0,+-10,0], and [0,0,+-10] 
    number_of_particles = 2                                             #Number of particles
    
    time = [0,100]                                                       #Time bounds. These aren't necessarily needed if you define dt but it's useful for plotting stuff
    t_step_number = 10                                                #Number of time steps. Always gotta define this or else things will break.
    
    dt = (time[-1])/t_step_number                                       #The delta of time. I usually define it as (time[-1] - time[0])/t_step_number 
    
    bounce_factor = .5                                                  #This is how bouncy the walls are. Can be any number but physically real values are between 0 and 1. With 0 being particles instantly stop at walls and 1 being the are perfectly reflected with no energy loss.
    
    
    args = generateParticles(number_of_particles, volume_bounds)
    mass_arr0, charge_arr0, radius_arr0, position_arr0, velocity_arr0, acceleration_arr0, Ke_arr0, Pe_arr0, Te_arr0 = args
    
    manager = Manager()
    
    t_mass_arr = [manager.list(mass_arr0)]
    t_charge_arr = [manager.list(charge_arr0)]
    t_radius_arr = [manager.list(radius_arr0)]  
    t_position_arr = [manager.list(position_arr0)]
    t_velocity_arr = [manager.list(velocity_arr0)]
    t_acceleration_arr = [manager.list(acceleration_arr0)]
    t_Ke_arr = [manager.list(Ke_arr0)]
    t_Pe_arr = [manager.list(Pe_arr0)]
    t_Te_arr = [manager.list(Te_arr0)]
    
    times = []
    total_start_time = tyme.time()
    for i in range(t_step_number):
        print(i)
        start_time = tyme.time()
        
        '''
        Make the multiprocess array
        '''
        multiprocess_arr = []
        for j in range(number_of_particles):
            index = j
            print("indices:", index)
            number_of_particles = number_of_particles
            dt = dt
            volume = volume_bounds
            mass_arr = t_mass_arr[-1]
            charge_arr = t_charge_arr[-1]
            radius_arr = t_radius_arr[-1]
            position_arr = t_position_arr[-1]
            velocity_arr = t_velocity_arr[-1]
            acceleration_arr = t_acceleration_arr[-1]
            Ke_arr = t_Ke_arr[-1]
            Pe_arr = t_Pe_arr[-1]
            Te_arr = t_Te_arr[-1]
            
            r_mass_arr = manager.list(np.zeros_like(mass_arr0))
            r_charge_arr = manager.list(np.zeros_like(charge_arr0))
            r_radius_arr = manager.list(np.zeros_like(radius_arr0))
            r_position_arr = manager.list(np.zeros_like(position_arr0))
            r_velocity_arr = manager.list(np.zeros_like(velocity_arr0))
            r_acceleration_arr = manager.list(np.zeros_like(acceleration_arr0))
            r_Ke_arr = manager.list(np.zeros_like(Ke_arr0))
            r_Pe_arr = manager.list(np.zeros_like(Pe_arr0))
            r_Te_arr = manager.list(np.zeros_like(Te_arr0))
            
            bounce_factor = bounce_factor
            
            temp = [index,
                number_of_particles,
                dt, 
                volume, 
                
                mass_arr, 
                charge_arr, 
                radius_arr, 
                position_arr, 
                velocity_arr, 
                acceleration_arr, 
                Ke_arr, 
                Pe_arr, 
                Te_arr,
                
                r_mass_arr, 
                r_charge_arr, 
                r_radius_arr, 
                r_position_arr, 
                r_velocity_arr, 
                r_acceleration_arr, 
                r_Ke_arr, 
                r_Pe_arr,
                r_Te_arr,
                
                bounce_factor
                ]
            multiprocess_arr.append(temp)
        print("Params Len:", len(multiprocess_arr))
        print(f'starting computations on {cpu_count()} cores')
        with Pool() as pool:
            res = pool.starmap(update_particle, multiprocess_arr)
            r_mass_arr, r_charge_arr, r_radius_arr, r_position_arr, r_velocity_arr, r_acceleration_arr, r_Ke_arr, r_Pe_arr,r_Te_arr = res[1]
        #print(manager2.join())
        print("r_mass:", r_mass_arr)
        print("r_charge:", r_charge_arr)
        #pool.join()
            
        t_mass_arr.append(np.asarray(list(r_mass_arr)))
        t_charge_arr.append(r_charge_arr)
        t_radius_arr.append(r_radius_arr)
        t_position_arr.append(np.asarray(r_position_arr))
        t_velocity_arr.append(r_velocity_arr)
        t_acceleration_arr.append(r_acceleration_arr)
        t_Ke_arr.append(r_Ke_arr)
        t_Pe_arr.append(r_Pe_arr)
        t_Te_arr.append(r_Te_arr)
        
        times.append(tyme.time()-start_time)
    total_end_time = tyme.time()
    
    t_mass_arr = np.asarray(t_mass_arr)
    t_charge_arr = np.asarray(t_charge_arr)
    t_radius_arr = np.asarray(t_radius_arr)
    
    t_position_arr = np.asarray(list(t_position_arr))
    t_velocity_arr = np.asarray(t_velocity_arr)
    t_acceleration_arr = np.asarray(t_acceleration_arr)
    t_Ke_arr = np.asarray(t_Ke_arr)
    t_Pe_arr = np.asarray(t_Pe_arr)
    t_Te_arr = np.asarray(t_Te_arr)
    
    #print("len t_mass_arr:", len(t_mass_arr))
    print("t_mass_arr:", str(t_mass_arr))
    print("t_charge_arr:", str(t_charge_arr))
    print("t_radius_arr:", str(t_radius_arr))
    print("t_position_arr:", str(t_position_arr))
    print("t_velocity_arr:", str(t_velocity_arr))
    print("t_acceleration_arr:", str(t_acceleration_arr))
    print("t_Ke_arr:", str(t_Ke_arr))
    print("t_Pe_arr:", str(t_Pe_arr))
    print("t_Te_arr:", str(t_Te_arr))
    
    colors = ["r", "b"]#plt.cm.rainbow(np.linspace(0,1, len(particles)))
    for i in range(number_of_particles):
        color = colors[i%len(colors)]
        plt.plot(t_position_arr[:][i][0], t_position_arr[:][i][1], color=color, label=str(i))
        plt.plot(t_position_arr[0][i][0], t_position_arr[0][i][1], "o", color=color)
        plt.plot(t_position_arr[-1][i][0], t_position_arr[-1][i][1], "*", color=color)
    plt.show()
        