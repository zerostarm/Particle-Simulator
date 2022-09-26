'''
Created on Sept 14, 2022

@author: Stephen
'''
import numpy as np
from scipy.ndimage import laplace #This function returns the laplacian of things correctly
import Multiprocessing_Library as mp
import time

class Mesh:
    '''
    classdocs
    '''
    def __init__(self, dv, volume_bounds, dt, particles, c=299792458, epsilon0=.001, bounce_factor=0.9):
        '''
        Constructor
        '''
        
        if type(dv) == type(1.0):
            dx, dy, dz = dv, dv, dv
        elif len(dv) == 2:
            dx, dy, dz = dv[0], dv[0], dv[1]
        else:
            dx, dy, dz = dv[0], dv[1], dv[2]
        self.dx, self.dy, self.dz = dx, dy, dz
        
        len_x = int(2*volume_bounds[0]/dx)
        len_y = int(2*volume_bounds[1]/dy)
        len_z = int(2*volume_bounds[2]/dz)
        
        self.dt = dt
        self.c = c
        self.mu0 = 1/(c**2*epsilon0)
        self.epsilon0 = epsilon0
        
        self.scalar_mesh_0 = np.zeros((len_x, len_y, len_z))
        self.vector_mesh_0 = np.zeros((3, len_x, len_y, len_z))
        
        self.psi = np.zeros_like(self.scalar_mesh_0)
        
        self.E_field = np.zeros_like(self.vector_mesh_0)
        self.B_field = np.zeros_like(self.vector_mesh_0)
        self.Scalar_field = np.zeros_like(self.vector_mesh_0)
        
        
        self.x = np.linspace(-volume_bounds[0], volume_bounds[0], len_x)
        self.y = np.linspace(-volume_bounds[1], volume_bounds[1], len_y)
        self.z = np.linspace(-volume_bounds[2], volume_bounds[2], len_z)
        
        self.particles = particles
        self.bounce_factor = bounce_factor
        self.volume_bounds = volume_bounds
        
        self.Standard_Test_fields()
        #self.test_psi()
        #self.test_plots()
    
    def Standard_Test_fields(self):
        """
        Fields that are not dependent on other particles must go at this level or in single particle test cases they won't be applied
        """
        import itertools
        for i, j, k in itertools.combinations(range(len(self.x)), 3):
            self.B_field[:, i, j, k] += np.asarray([0,0,100]) #Static B-Field Check
            #self.E_field[:, i, j, k] += np.asarray([0,0,1]) #Static E-field Check
            #self.Scalar_field[:, i, j, k] += -40/np.linalg.norm(particle1.getPosition())**3 *particle1.getPosition()
            #self.Scalar_field[:, i, j, k] += np.asarray([0,10,0]) #Non E&M field check
            #self.Scalar_field[:, i, j, k] += -0.5*particle1.getVelocity()**2
        
        
    
    def test_psi(self):
        self.psi = np.zeros_like(self.scalar_mesh_0)
        length = len(self.x)-1
        
        xx, yy, zz = np.meshgrid(self.x, self.y, self.z)
        
        self.psi = np.sqrt(xx**2 + yy**2 + zz**2)
        
        print("Avg Psi", np.average(self.psi))
        
        self.E_field = -1*np.asarray(np.gradient(self.psi))
        #self.B_field = 1*np.asarray(self.curl(np.gradient(self.psi)))
        
        return True
    
    def find_nearest_neighbors_better(self, x0, y0, z0):
        i_0 = np.absolute(np.floor(x0/self.dx))
        j_0 = np.absolute(np.floor(y0/self.dy))
        k_0 = np.absolute(np.floor(z0/self.dz))
        
        i_plus = np.absolute(np.ceil(x0/self.dx)) 
        j_plus = np.absolute(np.ceil(y0/self.dy))
        k_plus = np.absolute(np.ceil(z0/self.dz))
        
        if i_0 == i_plus:
            i_plus = i_0 + 1
        if j_0 == j_plus:
            j_plus = j_0 + 1
        if k_0 == k_plus:
            k_plus = k_0 + 1
        
        xi = self.x[0] + i_0*self.dx
        yj = self.y[0] + j_0*self.dy
        zk = self.z[0] + k_0*self.dz
        
        dx = self.dx
        dy = self.dy
        dz = self.dz
        
        A0 = (xi + dx - x0)*(yj + dy - y0)*(zk + dz - z0)
        A1 = (x0 - xi)*(yj + dy - y0)*(zk + dz - z0)
        A2 = (x0 - xi)*(y0 - yj)*(zk + dz - z0)
        A3 = (x0 - xi)*(y0 - yj)*(z0 - zk)
        A4 = (xi + dx - x0)*(y0 - yj)*(zk + dz - z0)
        A5 = (xi + dx - x0)*(y0 - yj)*(z0 - zk)
        A6 = (xi + dx - x0)*(yj + dy - y0)*(z0 - zk)
        A7 = (x0 - xi)*(yj + dy - y0)*(z0 - zk)
        At = dx*dy*dz
        
        weights = np.asarray([A0, A1, A2, A3, A4, A5, A6, A7])/At
        
        indices = np.asarray([[i_0, j_0, k_0],
                              [i_plus, j_0, k_0],
                              [i_plus, j_plus, k_0],
                              [i_plus, j_plus, k_plus],
                              [i_0, j_plus, k_0],
                              [i_0, j_plus, k_plus],
                              [i_0, j_0, k_plus],
                              [i_plus, j_0, k_plus]], dtype=int)
        
        
        return indices, weights
    
    def get_fields_at_point(self, particles, particle):
        temp_position = particle.getPosition()
        x0, y0, z0 = temp_position
        #print("Temp_position:", temp_position)
        indices, weights = self.find_nearest_neighbors_better(x0, y0, z0)
    
        point_Efield = np.zeros((3))
        point_Bfield = np.zeros((3))
        point_Scalarfield = np.zeros((3))
        for i in range(len(indices)):
            point_Efield += self.E_field[:, indices[i,0], indices[i,1], indices[i,2]]*weights[i] 
            point_Bfield += self.B_field[:, indices[i,0], indices[i,1], indices[i,2]]*weights[i]
            point_Scalarfield += self.Scalar_field[:, indices[i,0], indices[i,1], indices[i,2]]*weights[i]
        return point_Efield, point_Bfield, point_Scalarfield
    
    def getE(self):
        return self.E_field
    def getB(self):
        return self.B_field
    
    def curl(self, vector_field):
        dx = self.dx
        dy = self.dy
        dz = self.dz
        
        u, v, w = np.asarray(vector_field)
        
        dummy, dFx_dy, dFx_dz = np.gradient (u, dx, dy, dz, axis=[1,0,2])
        dFy_dx, dummy, dFy_dz = np.gradient (v, dx, dy, dz, axis=[1,0,2])
        dFz_dx, dFz_dy, dummy = np.gradient (w, dx, dy, dz, axis=[1,0,2])
    
        rot_x = dFz_dy - dFy_dz
        rot_y = dFx_dz - dFz_dx
        rot_z = dFy_dx - dFx_dy
    
        l = np.sqrt(np.power(u,2.0) + np.power(v,2.0) + np.power(w,2.0));
    
        m1 = np.multiply(rot_x,u)
        m2 = np.multiply(rot_y,v)
        m3 = np.multiply(rot_z,w)
    
        tmp1 = (m1 + m2 + m3)
        tmp2 = np.multiply(l,2.0)
    
        av = np.divide(tmp1, tmp2)
    
        return rot_x, rot_y, rot_z#, av
    
    def test_plots(self):
        import matplotlib.pyplot as plt
        """
        A bunch of test plots
        Uncomment the relevant ones you want to look at.
        """
        
        '''
        plt.plot(self.x, self.psi[:, 0, 0], label="Psi_x")
        plt.plot(self.x, self.psi[0, :, 0], label="Psi_y")
        plt.plot(self.x, self.psi[0, 0, :], label="Psi_z")
        plt.legend()
        plt.show()
        '''
        
        xx, yy, zz = np.meshgrid(self.x, self.y, self.z)
        
        import itertools
        pos = itertools.permutations(self.x, 2)
        pos = np.asarray([item for item in pos])
        print(pos)
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.axes.set_xlabel("X")
        ax.axes.set_ylabel("Y")
        ax.axes.set_zlabel("Z")
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        #ax.plot(np.ravel(xx[:,:,0]), np.ravel(yy[:,:,0]),  np.ravel(self.psi[:, :, 0]), label="psi")
        ax.quiver(np.ravel(xx[:,:,0]), np.ravel(yy[:,:,0]), np.ravel(zz[:,:,0]), np.ravel(self.E_field[0,:,:,0]), np.ravel(self.E_field[1,:,:,0]), np.ravel(self.E_field[2,:,:,0]), label="E_field")
        
        plt.show()
        
        '''
        plt.plot(self.x, self.Ti[:, 0, 0], label="Ti_x")
        plt.plot(self.x, self.Ti[0, :, 0], label="Ti_y")
        plt.plot(self.x, self.Ti[0, 0, :], label="Ti_z")
        plt.legend()
        plt.show()
        
        plt.plot(self.x, self.Ni[:, :, 10], label="ni_x")
        #plt.plot(self.x, self.Ni[0, :, 0], label="ni_y")
        #plt.plot(self.x, self.Ni[0, 0, :], label="ni_z")
        plt.legend()
        plt.show()
        
        plt.plot(self.x, self.B_field[0, :, 0, 0], label="B_x_x")
        plt.plot(self.x, self.B_field[1, :, 0, 0], label="B_x_y")
        plt.plot(self.x, self.B_field[2, :, 0, 0], label="B_x_z")
        plt.legend()
        plt.show()
        
        plt.plot(self.x, self.E_field[0, :, 0, 0], label="E_x_x")
        plt.plot(self.x, self.E_field[1, :, 0, 0], label="E_x_y")
        plt.plot(self.x, self.E_field[2, :, 0, 0], label="E_x_z")
        plt.plot(self.x, np.linalg.norm(self.E_field, axis=0)[:, 0, 0], label="E_tot")
        plt.legend()
        plt.show()
        '''

def generateParticles(number, volume): #Generates Particles in the defined rectangular volume [Length-x, width-y, height-z]
    """
    Generates a *number* of particles in a box defined by *volume*
    returns a list of particle objects
    """
    rand = np.random
    particles = []
    for i in range(number):
        mass = 2#rand.choice(range(1,2), 1)#integer
        charge =  2 #rand.choice([-1,0,1], 1)[0]#+-1, 0 integer 
        position = [0,0,0]
        '''[asarray([1,0,0]), asarray([-1,0,0]),
                    asarray([0,1,0]), asarray([0,-1,0]), 
                    asarray([0,0,1]), asarray([0,0,-1])][i%6]'''
        """[[10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)],
                    [-10, rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)]][i%2]""" 
                #[rand.randint(-volume[0]+1, volume[0]-1), rand.randint(-volume[1]+1, volume[1]-1), rand.randint(-volume[2]+1, volume[2]-1)] #Three Position
                # [asarray([1,0,0]), asarray([-1,0,0]), asarray([0,0,1]), asarray([0,0,-1])][i%3]
        #print(position)
        radius = 0.01 #mass/2 #float
        
        t = 1.0
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
        velocity =   [vx, vy, vz] #Three Velocity [asarray([0,1,0]), asarray([0,-1,0]), asarray([1,0,0])][i%3]
            #[asarray([1,0,0]), asarray([0,0,0]), asarray([0,0,0])][i%3]
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
    '''
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
        particle1.setPe_directly(Pe, 0)'''
    return particles

if __name__ == "__main__":
    from Particle import *
#    from Main import generateParticles
    volume = [10, 10, 10]
    mesh = Mesh(0.5, volume, 0.1)
    particles = generateParticles(1, volume)
    #mesh.generate_rho(particles)
    #mesh.get_new_fields(particles[0], particles)
    print(mesh.get_fields_at_point(particles, particles[0]))
    #print(mesh.rho)

#last best ~ 4.2s avg per step 
