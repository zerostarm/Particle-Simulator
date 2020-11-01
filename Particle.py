'''
Created on May 5, 2020

@author: Stephen
'''
from numpy import asarray


class Particle:
    def __init__(self, mass, position, velocity, energy, acceleration, radius):
        self.Mass = mass
        self.Position = asarray(position)
        self.Velocity = asarray(velocity)
        self.Energy = asarray(energy)
        self.Acceleration = asarray(acceleration)
        self.Radius = radius 
    
    def setMass(self, mass):
        self.Mass = mass
        return True
    def getMass(self):
        return self.Mass
    
    def setPosition(self, position):
        self.Position = asarray(position)
        return True
    def getPosition(self):
        return self.Position
    
    def setVelocity(self, velocity):
        self.Velocity = asarray(velocity)
        return True
    def getVelocity(self):
        return self.Velocity
    
    def setEnergy(self, energy):
        self.Energy = asarray(energy)
        return True
    def getEnergy(self):
        return self.Energy
    
    def setAcceleratio(self, acceleration):
        self.Acceleration = asarray(acceleration)
        return True
    def getAcceleration(self):
        return self.Acceleration
    
    def setRadius(self, radius):
        self.Radius = radius
        return True
    def getRadius(self):
        return self.Radius
       