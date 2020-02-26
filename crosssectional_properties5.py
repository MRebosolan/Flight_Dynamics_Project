
"""
Created on Mon Feb 17 12:32:42 2020
@author: lizzy & richelle
"""
import numpy as np
import math as m
import matplotlib.pyplot as plt

aircraft = "B737" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.605  # m
la = 2.661  # m
x1 = 0.172  # m
x2 = 1.211  # m
x3 = 2.591  # m
xa = 0.35   # m
ha = 0.205  # m
tsk = 1.1/1000  # m
tsp = 2.8/1000  # m
tst = 1.2/1000  # m
hst = 16./1000   # m
wst = 19./1000   # m
nst = 15  # -
d1 = 0.01154  # m
d3 = 0.01840  # m
theta = m.radians(28)  # rad
P = 97.4*1000  # N
"airfoil plot"
R=ha/2          #radius 
length=m.pi*R**2

Ab=wst*tst+(hst-tst)*tst #m
l=2*m.sqrt((Ca-R)**2+(R)**2)+m.pi*R
b=l/nst
k=b/R #radians

z_coordinates=[]
y_coordinates=[]
z_spar=[]
y_spar=[]
A_st=(hst)*tst+tst*wst

"Stringer placement"
z_stringer=[] 
y_stringer=[]
stringer_array=[]
z_stringer.append(-R)
y_stringer.append(0)
stringer_array.append([A_st, 0,-R])
h=0
for n in range(0,10001):
    z=-R*m.cos(n*b/(R))
    y=R*m.sin(n*b/(R))
    if z<=0. and z>-R+0.01:
        z_stringer.append(z)
        z_stringer.append(z)
        y_stringer.append(y)
        y_stringer.append(-y)
        stringer_array.append([A_st, y,z])
        stringer_array.append([A_st, -y,z])
        h=h+2
    elif z>0:
        break
for n in range(0,10001):
     last_room=0.5*m.pi*R-(h/2)*b
     z1=0
     z2=(n-2)*(Ca-R)/10000
     y1=R
     y2=-R/(Ca-R)*z2+R
     l_p=m.sqrt((z1-z2)**2+(y1-y2)**2)
     z=round(z2,6)
     y=round(y2,6)
     for i in range (0,nst):
         if abs(l_p-(i*b-last_room))<=0.000025:
             z_stringer.append(z)
             z_stringer.append(z)
             y_stringer.append(y)
             y_stringer.append(-y)
             stringer_array.append([A_st, y,z])
             stringer_array.append([A_st, -y,z])
         elif z>(Ca-R):
                 break
   
stringers=np.array(stringer_array)
print(stringers)   
"airfoil placement"
for i in range(0,101):
    z=-R+R*i/(100)
    y=np.sqrt(R**2-z**2)
    z_coordinates.append(z)
    y_coordinates.append(y)
for i in range(-100,101):
    z=0
    y=i/100*R
    z_spar.append(z)
    y_spar.append(y)
for i in range(0,101):
    z=i*(Ca-R)/100
    y=-R/(Ca-R)*z+R
    z_coordinates.append(z)
    y_coordinates.append(y)

for i in range(0,101):
    z=(100-i)*(Ca-R)/100
    y=-R/(Ca-R)*z+R
    z_coordinates.append(z)
    y_coordinates.append(-y)

for i in range(0,101):
    z=-R+R*(100-i)/100
    y=np.sqrt(R**2-z**2)
    z_coordinates.append(z)
    y_coordinates.append(-y)

"Calculation Centroid"
A_sk1=m.pi*R*tsk
A_sk2=m.sqrt((Ca-R)**2+R**2)*tsk
A_sp=ha*tsp

A_tot=A_sk1+2*A_sk2+A_sp+nst*A_st
print("Total area:", A_tot, "m^2")

Az_st=A_st*sum(z_stringer) #m^2
Az_sk1=A_sk1*2*-R/(m.pi)  #m^2
Az_sk2=A_sk2*1/2*(Ca-R)   #m^2
Ay_st=A_st*sum(y_stringer)
Ay_sk1=A_sk1*0
Ay_sk2=A_sk2*1/2*(R)
Ay_sk3=A_sk2*-1/2*(R)
zc=(Az_st+Az_sk1+2*Az_sk2)/A_tot
yc=(Ay_st+Ay_sk1+Ay_sk2+Ay_sk3)/A_tot          #due to symmetry it should be 0
print("centroid",zc,yc)
print(l)
"airfoil plot"

plt.plot(z_coordinates,y_coordinates, color='black',label="Periphery")
plt.plot(z_spar,y_spar,color='green', label="Spar")
plt.xlim(-0.2,0.6)
plt.ylim(-0.2,0.3)
plt.scatter(z_stringer,y_stringer, s=30, color='red',label="Stringers")
plt.scatter(zc, yc, marker="x", color='orange',label="Centroid")
plt.legend()

"Moment of Inertia results"
#Calculate moment of inertia of straght skin parts
angle = m.atan(R/(Ca-R)) #radians
length_beam = m.sqrt(R**2+(Ca-R)**2) #m
Iyy_straight = ((tsk*length_beam**(3)*m.cos(angle)*m.cos(angle))/12+A_sk2*((Ca-R)/2-zc)**2) #2 beams
#print("Moment of inertia of straight skin parts Iyy is:", Iyy_straight, "m^4")

Izz_straight = ((tsk*length_beam**(3)*m.sin(angle)*m.sin(angle))/12+A_sk2*(R/2)**2) #2 beams
#print("Moment of inertia of straight skin parts Izz is:", Izz_straight, "m^4")

#Calculate moment of inertia of arc skin part
Izz_arc = 1/2*m.pi*tsk*R**3         #
Iyy_arc = 1/2*m.pi*tsk*R**3-4/m.pi*R**3*tsk+A_sk1*(2*R/m.pi+zc)**2
#print("Moment of inertia of arc skin part Iyy is:", Iyy_arc, "m^4")
#print("Moment of inertia of arc skin part Izz is:", Izz_arc, "m^4")

#Calculate moment of inertia of spar
Izz_spar = (tsp*ha**3)/12
Iyy_spar = A_sp*(zc)**2

#print("Moment of inertia of spar Iyy is:", Iyy_spar,"m^4")
#print("Moment of inertia of spar Izz is:", Izz_spar, "m^4")

#Calculate moment of inertia of stiffners
values_Ady=[]
for i in y_stringer:
    dy=i**2
    Ady=A_st*dy
    values_Ady.append(Ady)
Izz_st=sum(values_Ady)

values_Adz=[]
for i in z_stringer:
    dz=(i-zc)**2
    Adz=A_st*dz
    values_Adz.append(Adz)
Iyy_st=sum(values_Adz)

#print("Moment of inertia of stiffners Izz is:", Izz_st, "m^4") 
#print("Moment of inertia of stiffners Iyy is:", Iyy_st, "m^4")
#Calculate total moment of inertia 
Iyy_total = 2*Iyy_straight + Iyy_arc + Iyy_spar + Iyy_st
Izz_total = 2*Izz_straight + Izz_arc + Izz_spar + Izz_st

print("Total moment of inertia Iyy =", Iyy_total, "m^4")
print("Total moment of inertia Izz =", Izz_total, "m^4") 

"torsional constant"
A1 = m.pi * R ** 2 / 2.
A2 = (Ca - R) * R
l2=m.sqrt((Ca-R)**2+R**2)

A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
b = np.array([0., 0., 0.])

### First row
A[0, 0] = 2. * A1
A[0, 1] = 2. * A2
b[0] = 1

### Second row
A[1, 0] = (R * m.pi / tsk + 2 * R / tsp) / (2 * A1)
A[1, 1] = (-2 * R / tsp) / (2 * A1)
A[1, 2] = -1.
b[1] = 0.

### Third row
A[2, 0] = (-2 * R / tsp) / (2 * A2)
A[2, 1] = (2 * l2 / tsk + 2 * R / tsp) / (2 * A2)
A[2, 2] = -1
b[2] = 0.

solution = np.linalg.solve(A, b)
J = 1. / solution[-1]
print(J)
