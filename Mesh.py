'''
Created on Jul 6, 2022

@author: Stephen
'''
import numpy as np
from scipy.ndimage import laplace #This function returns the laplacian of things correctly

class Mesh:
    '''
    classdocs
    '''
    def __init__(self, dv, volume_bounds, dt, c=299792458, epsilon0=.001):
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
        
        self.rho = np.zeros_like(self.scalar_mesh_0)
        
        self.E_field = np.zeros_like(self.vector_mesh_0)
        self.B_field = np.zeros_like(self.vector_mesh_0)
        
        self.A_field = np.zeros_like(self.vector_mesh_0)
        self.V_field = np.zeros_like(self.scalar_mesh_0)
        
        self.new_E = np.zeros_like(self.vector_mesh_0)
        self.new_B = np.zeros_like(self.vector_mesh_0)
        self.new_A = np.zeros_like(self.vector_mesh_0)
        self.new_V = np.zeros_like(self.scalar_mesh_0)
        
        self.A_list = [self.A_field, self.new_A]
        self.V_list = [self.V_field, self.new_V]
        
        self.E_list = [self.E_field, self.new_E]
        self.B_list = [self.B_field, self.new_B]
        
        self.J_list = [self.vector_mesh_0]
        self.rho_list = []
        
        self.x = np.linspace(-volume_bounds[0], volume_bounds[0], len_x)
        self.y = np.linspace(-volume_bounds[1], volume_bounds[1], len_y)
        self.z = np.linspace(-volume_bounds[2], volume_bounds[2], len_z)
    
    def find_nearest_neighbors(self, x0, y0, z0):
        xm = np.abs(self.x-x0).argmin()
        ym = np.abs(self.y-y0).argmin()
        zm = np.abs(self.z-z0).argmin()
        
        #print(xm, ym, zm)
        
        temp_x = np.delete(self.x, xm)
        temp_y = np.delete(self.y, ym)
        temp_z = np.delete(self.z, zm)
        
        xmp, ymp, zmp = self.x[xm], self.y[ym], self.z[zm]
        #print(xmp, ymp, zmp)
        
        xp = np.abs(temp_x-x0).argmin()
        yp = np.abs(temp_y-y0).argmin()
        zp = np.abs(temp_z-z0).argmin()
        
        #print(xp, yp, zp)
        
        xpp, ypp, zpp = temp_x[xp], temp_y[yp], temp_z[zp]
        
        #print(xpp, ypp, zpp)
        
        points = [[xmp, ymp, zmp],
                  [xmp, ymp, zpp],
                  [xmp, ypp, zmp],
                  [xmp, ypp, zpp],
                  [xpp, ymp, zmp],
                  [xpp, ymp, zpp],
                  [xpp, ypp, zmp],
                  [xpp, ypp, zpp]]
        
        points = np.asarray(points)
        
        #print(points)
        
        xi = np.where(self.x == xpp)[0][0]
        yi = np.where(self.y == ypp)[0][0]
        zi = np.where(self.z == zpp)[0][0]
        
        indices = [[xm, ym, zm],
                   [xm, ym, zi],
                   [xm, yi, zm],
                   [xm, yi, zi],
                   [xi, ym, zm],
                   [xi, ym, zi],
                   [xi, yi, zm],
                   [xi, yi, zi]]
        indices = np.asarray(indices)
        #print(indices)
        
        weights = np.sqrt((points[:,0]-x0)**2 + (points[:,1]-y0)**2 + (points[:,2]-z0)**2)
        weights = weights/np.sum(weights)
        
        #print(weights)
        
        return indices, weights
    
    def generate_rho(self, particles):
        """
        Will take in particles and generate a rho map 
        """
        
        #points, indices, weights = find_nearest_neighbors(0.5, 0.5, 0.5)
        
        rho = np.zeros_like(self.scalar_mesh_0)
        current = np.zeros_like(self.vector_mesh_0)
        for particle in particles:
            x, y, z = particle.getPosition()
            indices, weights = self.find_nearest_neighbors(x, y, z)
            
            for i in range(len(indices)):
                rho[indices[i,0], indices[i,1], indices[i,2]] += weights[i]*particle.getCharge()
                for j in range(3):
                    current[j,indices[i,0], indices[i,1], indices[i,2]] += weights[i]*particle.getCharge()*particle.getVelocity()[j]
        #self.rho = rho
        return rho, current
    
    def get_fields_at_point(self, x0, y0, z0):
        indices, weights = self.find_nearest_neighbors(x0, y0, z0)
        
        point_Efield = np.zeros((3))
        point_Bfield = np.zeros((3))
        for i in range(len(indices)):
            point_Efield += self.E_field[:, indices[i,0], indices[i,1], indices[i,2]]*weights[i]
            point_Bfield += self.B_field[:, indices[i,0], indices[i,1], indices[i,2]]*weights[i]
        return point_Efield, point_Bfield
        
    def get_new_fields(self, particles):
        rho_new, J_new = self.generate_rho(particles)
        
        last_A = self.A_field
        last_V = self.V_field
        
        current_A = self.new_A
        current_V = self.new_V
        
        grad2A = self.gradient_squared_A(current_A, edge_type="nearest", cval=0)/self.c**2
        grad2V = self.gradient_squared_V(current_V, edge_type="nearest", cval=0)/self.c**2
        
        #print("c gradV", grad2V.shape)
        #print(grad2V)
        
        A_next = grad2A*self.dt**2*self.c**2 + 2*current_A - last_A + self.mu0*J_new*self.dt**2*self.c**2
        V_next = grad2V*self.dt**2*self.c**2 + 2*current_V - last_V + self.epsilon0**-1*rho_new*self.dt**2*self.c**2
        
        #print("next A", A_next.shape)
        #print("Next V", V_next.shape)
        
        self.A_field = np.copy(self.new_A)
        self.V_field = np.copy(self.new_V)
        
        self.new_A = A_next
        self.new_V = V_next
        
        #self.A_list.append(A_next)
        #self.V_list.append(V_next)
        
        E = self.get_E(V_next, self.A_list, self.dt)
        B = self.get_B(A_next)
        
        #self.E_list.append(E)
        #self.B_list.append(B)
        
        self.rho_list.append(rho_new)
        self.J_list.append(J_new)
        
        self.E_field = E
        self.B_field = B
        
        return E, B
        
    def get_E(self, V, A_list, dt):
    
        #print(np.gradient(V))
        #print(np.diff(A, axis=1))
        
        grad_V = np.asarray(np.gradient(V))
        #print(grad_V.shape)
        
        #der_A = np.asarray(np.gradient(A, axis=0))
        der_A = (self.new_A - self.A_field)/dt
        
        #print(der_A.shape)
        
        return -grad_V - der_A
    
    def get_B(self, A):
        if np.shape(A)[0] == 2:
            der_A = np.gradient(A, axis=0)
        else:
            x_hat = np.asarray(np.gradient(A[2, :, :])[1] - np.gradient(A[1, :, :])[2]) #dAz/dy - dAy/dz
            y_hat = np.asarray(np.gradient(A[2, :, :])[0] - np.gradient(A[0, :, :])[2]) #dAz/dx - dAx/dz
            z_hat = np.asarray(np.gradient(A[1, :, :])[0] - np.gradient(A[0, :, :])[1]) #dAy/dx - dAx/dy
        
            der_A = np.asarray([x_hat, -y_hat, z_hat])
        #print("All_grad_A", np.shape(All_grad_A))
        #print("Der_A", der_A.shape)
        return der_A
    
    
    def gradient_squared_A(self, A, edge_type="nearest", cval=0):
        temp = []
        for i in range(np.shape(A)[0]):
            temp.append(laplace(A[i], mode=edge_type, cval=cval))
        return np.asarray(temp)
    
    def gradient_squared_V(self, V, edge_type="nearest", cval=0):
        return laplace(V, mode=edge_type, cval=cval)
    
    def find_grad_dot(self, E):
        if np.shape(E)[0] == 2:
            temparr = np.zeros_like(E[0])
            for i in range(np.shape(E)[0]):
                temparr += np.gradient(E[i, :, :])[i]
        elif np.shape(E)[0] == 3:
            temparr = np.zeros_like(E[0])
            for i in range(np.shape(E)[0]):
                temparr += np.gradient(E[i, :, :, :])[i]
        else:
            temparr = np.zeros_like(E[0])
        
        return np.asarray(temparr) 
    
    def find_grad_cross(self, A):
        if np.shape(A)[0] == 2:
            der_A = np.gradient(A, axis=0)
        else:
            x_hat = np.asarray(np.gradient(A[2, :, :])[1] - np.gradient(A[1, :, :])[2]) #dAz/dy - dAy/dz
            y_hat = np.asarray(np.gradient(A[2, :, :])[0] - np.gradient(A[0, :, :])[2]) #dAz/dx - dAx/dz
            z_hat = np.asarray(np.gradient(A[1, :, :])[0] - np.gradient(A[0, :, :])[1]) #dAy/dx - dAx/dy
            der_A = np.asarray([x_hat, -y_hat, z_hat]) 
        return der_A
    
    def getE(self):
        return self.E_field
    def getB(self):
        return self.B_field
    
if __name__ == "__main__":
    from Particle import *
    from Main import generateParticles
    volume = [10, 10, 10]
    mesh = Mesh(0.5, volume, 0.1)
    particles = generateParticles(1, volume)
    #mesh.generate_rho(particles)
    mesh.get_new_fields(particles)
    print(mesh.getE())
    #print(mesh.rho)
