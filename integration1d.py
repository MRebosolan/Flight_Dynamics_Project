import numpy as np
from interpolation1d import interpolation
def multiplicator(loadfunction):
    xprime1 = 0
    xprime2 = 0.5464859422096491
    xprimes = np.linspace(xprime1, xprime2, pointnum)
    xhinge = (0.225)/2
    newfunction = []
    for c, element in enumerate(loadfunction):
        new = element*(xhinge - xprimes[c])
        newfunction.append(new)
    return newfunction

def numerical_integrator1d(loadfunction,pointnum):
    xprime1 = 0
    xprime2 = 0.5464859422096491
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
y = i.interpolate(ill,20)
print()
integral = numerical_integrator1d(y,pointnum)
print(integral)