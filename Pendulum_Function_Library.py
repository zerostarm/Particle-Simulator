'''
Created on Dec 15, 2020

@author: Stephen
'''
import numpy as np

g = 9.81 #m/s^2

def x(r, theta):
    return r*np.cos(theta-np.pi/2)
    
def y(r, theta):
    return r*np.sin(theta-np.pi/2)


def PendulumOnRotatingPlate(n, Tbounds, dt, alpha, b):
    w = 1 #rads/sec
    
    theta0 = 0 #rads
    thetaDot0 = -.1 #rads/sec
    
    
    def thetaDoubleDotNext(theta, t):
        return (-g*np.sin(theta) + alpha*w**2*np.cos(theta-w*t))/b
    
    def xu(r,theta, u, phi):
        return r*np.cos(theta-np.pi/2) + u*np.cos(phi-np.pi/2)
    
    def yu(r,theta, u, phi):
        return r*np.sin(theta-np.pi/2) + u*np.sin(phi-np.pi/2)
    
    t = np.linspace(Tbounds[0], Tbounds[1], int(n))
    
    
    theta = np.zeros(t.shape)
    theta[0] = theta0
    thetaDot = np.zeros(t.shape)
    thetaDot[0] = thetaDot0
    thetaDoubleDot = np.zeros(t.shape)
    thetaDoubleDot[0] = 0
    
    for i in range(len(t)-1):
        thetaDoubleDot[i+1] = thetaDoubleDotNext(theta[i], t[i])
        thetaDot[i+1] = thetaDot[i]+thetaDoubleDot[i+1]*dt
        theta[i + 1] = theta[i]+ thetaDot[i+1]*dt
    
    return x(alpha, w*t), y(alpha, w*t), xu(b, theta, alpha, w*t), yu(b, theta, alpha, w*t)


def PendulumSpring(n, Tbounds, k, m, l, l0, y0, theta0, thetaDot0, r0, rDot0):
    dt = (Tbounds[1] - Tbounds[0])/n
    t = np.linspace(Tbounds[0], Tbounds[1], int(n))
    
    y0 *= -1
    
    def thetaDoubleDotNext(theta, thetaDot, r, rDot):
        #r = np.abs(r)
        if np.abs(r > l):
            return -g/(l0+r)*np.sin(theta)-2*rDot/(l0+r)*thetaDot
        else:
            return (-g/(l0+r)*np.sin(theta)-2*rDot/(l0+r)*thetaDot)
        return 0
    def rDoubleDotNext(theta, thetaDot, r):
        #r = np.abs(r)
        if np.abs(r > l):
            return (l0+r)*thetaDot**2+g*np.cos(theta)
        else:
            return (l0+r)*thetaDot**2+g*np.cos(theta)-k/m*r
        return 0
    theta = np.zeros(t.shape)
    theta[0] = theta0
    thetaDot = np.zeros(t.shape)
    thetaDot[0] = thetaDot0
    thetaDoubleDot = np.zeros(t.shape)
    thetaDoubleDot[0] = 0
    
    r = np.zeros(t.shape)
    r[0] = r0 + l0
    rDot = np.zeros(t.shape)
    rDot[0] = rDot0
    rDoubleDot = np.zeros(t.shape)
    rDoubleDot[0] = 0
    
    for i in range(len(t)-1):
        thetaDoubleDot[i+1] = thetaDoubleDotNext(theta[i], thetaDot[i], r[i], rDot[i])
        thetaDot[i+1] = thetaDot[i]+thetaDoubleDot[i+1]*dt
        theta[i + 1] = theta[i]+ thetaDot[i+1]*dt
        
        rDoubleDot[i+1] = rDoubleDotNext(theta[i], thetaDot[i], r[i])
        rDot[i+1] = rDot[i] + rDoubleDot[i+1]*dt
        r[i+1] = r[i] + rDot[i+1]*dt

        if y(r[i+1], theta[i+1]) <= y0+1:
            rDoubleDot[i+1] = 2*r[i]*thetaDot[i]**2 - 4*rDot[i]*np.cos(theta[i])/np.sin(theta[i])*thetaDot[i]  
            rDot[i+1] = rDot[i] + rDoubleDot[i+1]*dt
            r[i+1] = r[i] + rDot[i+1]*dt
        
            thetaDoubleDot[i+1] = np.sin(theta[i]+np.pi/2)/(y0)*(rDoubleDot[i] + (r[i] +l0) * thetaDot[i]**2) 
            thetaDot[i+1] = thetaDot[i]+thetaDoubleDot[i+1]*dt
            theta[i + 1] = theta[i]+ thetaDot[i+1]*dt
    return x(r, theta), y(r, theta), r, t

def MultiprocessingStandard(plotStuff, values, prints=True):
    import multiprocessing
    if prints:
        print(f'starting computations on {multiprocessing.cpu_count()} cores')
    results = []
    with multiprocessing.Pool() as pool:
        res = pool.starmap(plotStuff, values)
        #print(res)
        results.append(res)
    #print(len(results))
    if prints:
        print("Multiprocessing Done")
    return True

def createVideo(path, Name, Tbounds, n):
    import cv2
    import os
    print("Start Video Creation")
    if not ".avi" in Name:
        Name += ".avi"
    dt = (Tbounds[1] - Tbounds[0])/n
    images = os.listdir(path)[:n]
    frame = cv2.imread(os.path.join(path,images[0]))
    height, width, layers = frame.shape
    video = cv2.VideoWriter(Name, 0, dt**-1, (width, height))
    for image in images:
        read = cv2.imread(os.path.join(path, str(image)))
        #cv2.imshow("", read)
        video.write(read)
    cv2.destroyAllWindows()
    video.release()
    print("End Video Creation")
    return True