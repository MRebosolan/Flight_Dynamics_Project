import numpy as np
from newinter import torquerloader3000
from math import *
from integration1d import numerical_integrator1d,multiplicator
load,torque = torquerloader3000()
print(torque)
print(sum(load))
xsamp = np.linspace(0,-2.76,100)
def qcomponents(load):
    qy = [element*sin(radians(64)) for element in load]
    qz = [element*cos(radians(64)) for element in load]
    return qy,qz
x1 = -0.153
x2 = -1.281
x3 = -2.681
xa = 0.28
d1 = 0.01103
d3 = 0.01642
P = 91.7
T = numerical_integrator1d(torque,100,-2.76)
print(T)
def multitron(load,xpoint):
    newlist = []
    for c,element in enumerate(load):
        lol = (element*xpoint[c])
        newlist.append(lol)
    return newlist



def xloadloc(xsamp):
    qy,qz = qcomponents(load)
    xz = (numerical_integrator1d(multitron(qz,xsamp),100,-2.76))/(numerical_integrator1d(qz,100,-2.76))
    xy = (numerical_integrator1d(multitron(qy,xsamp),100,-2.76))/(numerical_integrator1d(qy,100,-2.76))
    return xz,xy
xz,xy = xloadloc(xsamp)
def momenty(x):

def momentz(x):

def torque(x):
