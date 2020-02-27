import numpy as np
from newinter import torquerloader3000
import math
from integration1d import numerical_integrator1d,multiplicator
from matplotlib import pyplot as plt
load12,torque = torquerloader3000()

xsamp = np.linspace(0,-2.76,100)
def qcomponents(load):
    qy = [element*math.sin(math.radians(64)) for element in load]
    qz = [element*math.cos(math.radians(64)) for element in load]
    return qy,qz
qy,qz = qcomponents(load12)
L = 2.76
x1 = -0.153
x2 = -1.281
x3 = -2.681
xa = 0.28
ha = 0.225
d1 = 0.01103
d3 = 0.01642
P = 91.7
zsc = -0.11922705644352412 - (ha/2)
cos = math.cos(math.radians(26))
sin = math.sin(math.radians(26))
xp = (x2 + (xa/2))
xI = (x2 - (xa/2))
armact = ((ha/2)-(zsc*cos))
E = 72*10**6
Izz = 1.2807456661786271e-05
G = 28 * 10**6
J = 1.5101498390705797e-05
Iyy = 6.86413741465098e-05
ksi = math.sqrt((ha/2)**+zsc**2)
def multitron(load,xpoint):
    newlist = []
    for c,element in enumerate(load):
        lol = (element*xpoint[c])
        newlist.append(lol)
    return newlist
#doubleqy = numerical_integrator1d(qy,100,-2.77)


def xloadloc(xsamp):
    qy,qz = qcomponents(load)
    xz = (numerical_integrator1d(multitron(qz,xsamp),100,-2.76))/(numerical_integrator1d(qz,100,-2.76))
    xy = (numerical_integrator1d(multitron(qy,xsamp),100,-2.76))/(numerical_integrator1d(qy,100,-2.76))
    return xz,xy
#xz,xy = xloadloc(xsamp)
def xsamp1(num):
    return np.linspace(0,num,100)
def multitron2(load,xpoint,power):
    newlist = []
    for c,element in enumerate(load):
        lol = (element*(xpoint[c]**power))
        newlist.append(lol)
    return newlist
dubly = (1/2)*numerical_integrator1d(multitron(qy,xsamp),100,L)
dublz = (1/2)*numerical_integrator1d(multitron(qz,xsamp),100,L)
quadhinge1 = (1/24)*numerical_integrator1d(multitron2(qy,xsamp,3),100,x1)
dubhinge1 = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x1)
quadhinge2 = (1/24)*numerical_integrator1d(multitron2(qy,xsamp,3),100,x2)
dubhinge2 = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x2)
quadhinge3 = (1/24)*numerical_integrator1d(multitron2(qy,xsamp,3),100,x3)
dubhinge3 = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x3)
quadhinge1z = (1/24)*numerical_integrator1d(multitron2(qz,xsamp,3),100,x1)
dubhinge1z = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x1)
quadhinge2z = (1/24)*numerical_integrator1d(multitron2(qz,xsamp,3),100,x2)
dubhinge2z = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x2)
quadhinge3z = (1/24)*numerical_integrator1d(multitron2(qz,xsamp,3),100,x3)
dubhinge3z = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,x3)
quadp = (1/24)*numerical_integrator1d(multitron2(qz,xsamp,3),100,xp)
dubp = (1/2)*numerical_integrator1d(multitron(torque,xsamp),100,xp)
print(quadhinge1)
print(dublz)
vh1 = 0.01103*cos+1/E/Izz*(quadhinge1)-(dubhinge1)*zsc/G/J
vh2 =(quadhinge2)/E/Izz - zsc/G/J*(dubhinge2)
vh3 = 0.01642*cos+(quadhinge3)/E/Izz-zsc/G/J*(dubhinge3)+P/6/E/Izz*(x3-x2-xa/2)**3*sin+P*zsc/G/J*cos*(ha/2+zsc*sin)*(x3-x2-xa/2)+zsc*P*sin*zsc/G/J*(x3-x2-xa/2)
wh1 =  -0.01103*sin-(quadhinge1z)/E/Iyy
wh2 = -(quadhinge2z)/E/Iyy
wh3 = -0.01642*sin-(quadhinge3z)/E/Iyy+P*cos/E/Iyy/6*(x3-x2-xa/2)**3
wP = -(quadp)/E/Iyy-ksi*sin/G/J*(dubp)
mzl = -dubly - P*(L-xp)*sin
myl = -dublz - P*(L-xp)*cos
szl = - numerical_integrator1d(qz,100,L) - P*cos
syl = - numerical_integrator1d(qy,100,L) - P*sin
tl = - numerical_integrator1d(torque,100,L) + P*cos*(ha/2) + P*sin*zsc

A = np.array([
    [-1*cos*(L-x1),-1*cos*(L-x2),-1*cos*(L-x3),-1*sin*(L-x1),-1*sin*(L-x2),-1*sin*(L-x3),-1*sin*(L-xI),0,0,0,0,0],
    [-1*sin*(L-x1),-1*sin*(L-x2),-1*sin*(L-x3),-1*cos*(L-x1),-1*cos*(L-x2),-1*cos*(L-x3),-1*cos*(L-xI),0,0,0,0,0],
    [-1*sin,-1*sin,-1*sin,-1*cos,-1*cos,-1*cos,-1*cos,0,0,0,0,0],
    [-1*cos,-1*cos,-1*cos,-1*sin,-1*sin,-1*sin,-1*sin,0,0,0,0,0],
    [cos*zsc*(L-x1),cos*zsc*(L-x2),cos*zsc*(L-x3),sin*zsc*(L-x1),sin*zsc*(L-x2),sin*zsc*(L-x3),sin*zsc*xI+cos*(ha/2)*xI,0,0,0,0,0],
    [0, 0, 0, 0, 0, 0, 0, -x1/E/Izz, -1/E/Izz, zsc/(G*J), 0, 0],
    [cos*(x2-x1)**3/6/E/Izz+cos*zsc/G/J*zsc*(x2-x1), 0, 0, sin*(x2-x1)**3/6/E/Izz+zsc*sin/G/J*zsc*(x2-x1), 0, 0, (xa/2)**3*sin/E/Izz/6+zsc*sin/G/J*zsc*xa/2 + zsc*cos/G/J*(ha/2)*xa/2, -x2/E/Izz, -1/E/Izz, zsc/G/J, 0, 0],
    [cos*(x3-x1)**3/6/E/Izz+cos*zsc/G/J*zsc*(x3-x1), cos*(x3-x2)**3/6/E/Izz+cos*zsc/G/J*zsc*(x3-x2), 0, sin*(x3-x1)**3/6/E/Izz+sin*zsc/G/J*zsc*(x3-x1), sin*(x3-x2)**3/6/E/Izz+sin*zsc/G/J*zsc*(x3-x2), 0, (x3-x2+xa/2)**3*sin/6/E/Izz+sin*zsc/G/J*(x3-x2+xa/2)+cos*zsc/G/J*(ha/2)*(x3-x2+xa/2), -x3/E/Izz, -1/E/Izz, zsc/G/J, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x1/E/Iyy, -1/E/Iyy],
    [sin/6/E/Iyy*(x2-x1)**3, 0, 0, cos/6/E/Iyy*(x2-x1)**3, 0, 0, cos/6/E/Iyy*(xa/2)**3, 0, 0, 0, -x2/E/Iyy, -1/E/Iyy],
    [sin/6/E/Iyy*(x3-x1)**3, sin/6/E/Iyy*(x3-x2)**3, 0, cos/6/E/Iyy*(x3-x1)**3, cos/6/E/Iyy*(x3-x2)**3, 0, cos/6/E/Iyy*(x3-x2+xa/2)**3, 0, 0, 0, -x3/E/Iyy, -1/E/Iyy],
    [sin/E/Iyy/6*(xp-x1)**3+ksi*sin*cos*zsc/G/J*(xp-x1), sin/E/Iyy/6*(xa/2)**3+ksi*sin*cos*zsc/G/J*(xp-x2), 0, cos/E/Iyy/6*(xp-x1)**3+ksi*sin/G/J*sin*zsc*(xp-x1), cos/E/Iyy/6*(xa/2)**3+ksi*sin/G/J*sin*zsc*(xp-x2), 0, cos/E/Iyy/6*(xa)**3+ksi*sin/G/J*sin*zsc*(xa)+ksi*sin/G/J*cos*(ha/2)*(xa), 0, 0, sin*ksi/G/J, -xp/E/Iyy, -1/E/Iyy]
])


X = np.array([
    [mzl],
    [myl],
    [szl],
    [syl],
    [tl],
    [vh1],
    [vh2],
    [vh3],
    [wh1],
    [wh2],
    [wh3],
    [wP]
])

solution = np.linalg.solve(A,X)
print(solution)
def stepfunction(x,xhinge):
    if x - xhinge >= 0:
        return x-xhinge
    else:
        return 0
def stepfunction2(x,x2):
    if x - x2 >= 0:
        return 1
    else:
        return 0
def Momentz(x):
    dub = (1/2)*numerical_integrator1d(multitron(qy,xsamp),100,x)
    return dub-1*cos*stepfunction(x,x1)*solution[0]-1*cos*stepfunction(x,x2)*solution[1]-1*cos*stepfunction(x,x3)*solution[2]-1*sin*stepfunction(x,x1)*solution[3]-1*sin*stepfunction(x,x2)*solution[4]-1*sin*stepfunction(x,x3)*solution[5]-1*sin*stepfunction(x,xI)*solution[6] + P*stepfunction(x,xp)*sin
def Momenty(x):
    dub = (1/2)*numerical_integrator1d(multitron(qz,xsamp),100,x)
    return dub-1*sin*stepfunction(x,x1)*solution[0]-1*sin*stepfunction(x,x2)*solution[1]-1*sin*stepfunction(x,x3)*solution[2]-1*cos*stepfunction(x,x1)*solution[3]-1*cos*stepfunction(x,x2)*solution[4]-1*cos*stepfunction(x,x3)*solution[5]-1*cos*stepfunction(x,xI)*solution[6] + P*stepfunction(x,xp)*cos
def Shearz(x):
    return numerical_integrator1d(qz,100,x) -1*sin*solution[0]*stepfunction2(x,x1)-1*sin*solution[1]*stepfunction2(x,x2)-1*sin*solution[2]*stepfunction2(x,x3)-1*cos*solution[3]*stepfunction2(x,x1)-1*cos*solution[4]*stepfunction2(x,x2)-1*cos*solution[5]*stepfunction2(x,x3)-1*cos*solution[6]*stepfunction2(x,xI) + P*cos*stepfunction2(x,xp)
def Sheary(x):
    return numerical_integrator1d(qy,100,x) -1*cos*solution[0]*stepfunction2(x,x1)-1*cos*solution[1]*stepfunction2(x,x2)-1*cos*solution[2]*stepfunction2(x,x3)-1*sin*solution[3]*stepfunction2(x,x1)-1*sin*solution[4]*stepfunction2(x,x2)-1*sin*solution[5]*stepfunction2(x,x3)-1*sin*solution[6]*stepfunction2(x,xI) + P*sin*stepfunction2(x,xp)
def Torque(x):
    return numerical_integrator1d(torque,100,x)+cos*zsc*stepfunction2(x,x1)*solution[0]+cos*zsc*stepfunction2(x,x2)*solution[1]+cos*zsc*stepfunction2(x,x3)*solution[2]+sin*zsc*stepfunction2(x,x1)*solution[3]+sin*zsc*stepfunction2(x,x2)*solution[4]+sin*zsc*stepfunction2(x,x3)*solution[5]+sin*zsc*stepfunction2(x,xI)*solution[6]+cos*(ha/2)*stepfunction2(x,xI)*solution[6] - P*cos*(ha/2)*stepfunction2(x,xp) - P*sin*zsc*stepfunction2(x,xp)
y23 = []
y24 = []
y25 = []
y26 = []
y27 = []
for element in xsamp:
    y23.append(Momentz(element))
    y24.append(Momenty(element))
    y25.append(Shearz(element))
    y26.append(Sheary(element))
    y27.append(Torque(element))

plt.plot(xsamp,y27)
plt.show()