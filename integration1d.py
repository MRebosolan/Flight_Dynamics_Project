import numpy as np
from math import *
import matplotlib.pyplot as plt
from interpolation1d import interpolation
def multiplicator(loadfunction):
    xprime1 = 0
    xprime2 = -0.5464859422096491
    xprimes = np.linspace(xprime1, xprime2, pointnum)
    xhinge = -0.23172705644352412
    newfunction = []
    for c, element in enumerate(loadfunction):
        new = element*( xhinge - xprimes[c])
        newfunction.append(new)
    return newfunction

def numerical_integrator1d(loadfunction,pointnum,xprime2):
    xprime1 = 0
    xprimes = np.linspace(xprime1, xprime2, pointnum)  # array of x coordinates

    tot = 0

    for c,element in enumerate(loadfunction[:-1]):
                #i is index of element in the row
        tot += (loadfunction[c] + loadfunction[c+1]) * (xprimes[c+1] - xprimes[c]) * 0.5
    return tot

pointnum = 400
pointnum2 = pointnum - 1
i = interpolation()
ill = np.linspace(0,-0.5456,pointnum)
magni = []
def resultandtorques():
    resultants = []
    torques = []
    for g in range(41):
        y = i.interpolate(ill,g,'span',magni)
        integral = numerical_integrator1d(y,pointnum,-0.5464859422096491)
        y2 = multiplicator(y)
        twistint = numerical_integrator1d(y2,pointnum,-0.5464859422096491)
        resultants.append(integral)
        torques.append(twistint)
    return resultants,torques
