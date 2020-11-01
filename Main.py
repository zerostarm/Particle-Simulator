'''
Created on May 5, 2020

@author: Stephen
'''

import Particle
import numpy as np

C = 299792458 #speed of light in vacuum

def generateParticles(number, volume): #Generates Particles in the defined rectangular volume [Length-x, width-y, height-z]
    rand = np.random
    particles = []
    for i in range(number):
        
        mass = rand.choice([1,16], 1)#integer
        position = [rand.randint(0, volume[0]), rand.randint(0, volume[1]), rand.randint(0, volume[2])] #Three Position
        vx = np.sqrt(rand.randint(0, C**2))
        vy = np.sqrt(rand.randint(0, C**2 - vx**2))
        vz = np.sqrt(rand.randint(0, C**2 - vx*2 - vy**2))
        velocity = [vx, vy, vz] #Three Velocity
        linearE = mass*np.sum(velocity**2)
        energy = [mass*C**2, linearE] #Three or more energy idk, for now just 2 energy
        acceleration = [] #Three Acceleration
        radius = mass #integer
        
        particles.append()