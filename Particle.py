'''
Created on Sept 14, 2020

@author: Stephen
'''

from numpy import asarray
import numpy as np

class Particle:
    def __init__(self, index, mass, charge, position, velocity, radius=0, Ke=0, acceleration=[0,0,0], Pe=0):
        self.Index = index
        self.Mass = mass
        self.Charge = charge
        self.Position = [asarray(position, dtype=np.float64)]
        self.Velocity = [asarray(velocity, dtype=np.float64)]
        self.KineticEnergy = [Ke]
        self.PotentialEnergy = [Pe]
        self.TotalEnergy = [Ke + Pe]
        self.Acceleration = [asarray(acceleration, dtype=np.float64)]
        self.Radius = radius 
        
        self.NewPosition = asarray(position)
        self.NewVelocity = asarray(velocity)
        self.NewEnergy = Ke + Pe
        self.NewAcceleration = asarray(acceleration)
        
        self.indices = np.asarray([[0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0]])
        
        self.weights = np.asarray([[0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0],
                                   [0,0,0]])
        
        
    def setMass(self, mass):
        self.Mass = mass
        return True
    def getMass(self):
        return self.Mass
    
    def setCharge(self, charge):
        self.Charge = charge
        return True
    def getCharge(self):
        return self.Charge
    
    def setPosition(self):
        self.Position.append(self.NewPosition)
        return True
    def setNewPosition(self, position):
        self.NewPosition = position
        return True
    def setPosition_Directly(self, position):
        self.Position.append(position)
        return True
    def getPosition(self, i=-1):
        return self.Position[i]
    
    def setVelocity(self):
        self.Velocity.append(self.NewVelocity)
        return True
    def setNewVelocity(self, velocity):
        self.NewVelocity = asarray(velocity)
        return True
    def setVelocity_Directly(self, velocity):
        self.Velocity.append(asarray(velocity))
        return True
    def getVelocity(self, i =-1):
        return self.Velocity[i]
    
    def setKe(self):
        self.KineticEnergy.append(self.NewKineticEnergy)
        return True
    def setNewKe(self, energy):
        self.NewKineticEnergy = energy
        return True
    def setPe(self):
        self.PotentialEnergy.append(self.NewPotentialEnergy)
        return True
    
    def setPe_directly(self, Pe, i):
        self.PotentialEnergy[i] = Pe
        return True
    
    def setNewPe(self, Pe):
        self.NewPotentialEnergy = Pe
        return True
    def setNewTe(self, Te):
        self.NewTotalEnergy = Te
        return True
    def setTe(self):
        self.TotalEnergy.append(self.NewTotalEnergy)
        return True
    
        
        
    def setAcceleration(self):
        self.Acceleration.append(self.NewAcceleration)
        return True
    def setNewAcceleration(self, acceleration):
        self.NewAcceleration = acceleration
        return True
    def setAcceleration_Directly(self, acceleration):
        self.Acceleration.append(asarray(acceleration))
        return True
    def getAcceleration(self, i=-1):
        return self.Acceleration[i]
    
    def setRadius(self, radius):
        self.Radius = radius
        return True
    def getRadius(self):
        return self.Radius
    
    def setIndex(self, index):
        self.Index = index
        return True
    def getIndex(self):
        return self.Index
    
    def __str__(self):
        return str(self.Mass) + " " + str(self.Charge) + " " + str(self.Radius) + " " + str(self.TotalEnergy[-1]) 
        
    