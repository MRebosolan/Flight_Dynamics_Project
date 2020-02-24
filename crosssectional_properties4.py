"""
Created on Mon Feb 17 12:32:42 2020
@author: lizzy & richelle
"""
import numpy as np
import math as m
import matplotlib.pyplot as plt

Ca = 0.547  # m
la = 2.771  # m
x1 = 0.153  # m
x2 = 1.281  # m
x3 = 2.681  # m
xa = 0.28   # m
ha = 0.225  # m
tsk = 1.1/1000  # m
tsp = 2.9/1000  # m
tst = 1.2/1000  # m
hst = 15./1000   # m
wst = 20./1000   # m
nst = 17 # -
d1 = 0.01103  # m
d3 = 0.01642  # m
theta = m.radians(26)  # rad
P = 91.7*1000  # N

"airfoil plot"
R=ha/2          #radius 
length=m.pi*R**2

Ab=wst*tst+(hst-tst)*tst #m
l=2*m.sqrt((Ca-R)**2+(R)**2)+m.pi*R
b=l/17
k=b/R #radians

z_coordinates=[]
y_coordinates=[]
A_st=(hst)*tst+tst*wst

"Stringer placement"
z_stringer=[] 
y_stringer=[]
stringer_array=[]
for n in range(0,1000000):
    z=-R*m.cos(n*b/(R))
    y=R*m.sin(n*b/(R))
    if n>3:
        last_room=0.5*m.pi*R-2*b
        z1=0
        z2=(n-2)*(Ca-R)/1000000
        y1=R
        y2=-R/(Ca-R)*z2+R
        l_p=m.sqrt((z1-z2)**2+(y1-y2)**2)
        z=z2
        y=y2
        for i in range (0,10):
            if abs(l_p-(i*b-last_room))<0.0000002:
                z_stringer.append(z)
                z_stringer.append(z)
                y_stringer.append(y)
                y_stringer.append(-y)
                stringer_array.append([A_st, z,y])
                stringer_array.append([A_st, z,-y])
            elif z>(Ca-R):
                break
    elif n<3:
        z_stringer.append(z)
        y_stringer.append(y)
        stringer_array.append([A_st, z,y])      
        if n>0:
            z_stringer.append(z)
            y_stringer.append(-y)
            stringer_array.append([A_st, z,-y])

print(stringer_array)
       
"airfoil placement"
for i in range(0,r):
    z=-R+R*i/(r-1)
    y=np.sqrt(R**2-z**2)
    z_coordinates.append(z)
    y_coordinates.append(y)

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

"airfoil plot"
plt.plot(z_coordinates,y_coordinates, color='black')
plt.xlim(-0.2,0.6)
plt.ylim(-0.2,0.3)
plt.scatter(z_stringer,y_stringer, s=30, color='red')
plt.scatter(zc, yc, marker='x', color='orange')

"Calculation Centroid"
A_sk1=m.pi*R*tsk
A_sk2=m.sqrt((Ca-R)**2+R**2)*tsk
A_sp=ha*tsp

A_tot=A_sk1+2*A_sk2+A_sp+17*A_st
print("Total area:", A_tot, "m^2")

Az_st=A_st*sum(z_stringer)
Az_sk1=A_sk1*2*-R/(m.pi)
Az_sk2=A_sk2*1/2*(Ca-R)
zc=(Az_st+Az_sk1+2*Az_sk2)/A_tot
yc=0                            #due to symmetry
#recalculated again due to big difference with verification model
print("centroid",zc)
"Moment of Inertia results"
#Calculate moment of inertia of straght skin parts
angle = m.atan(R/(Ca-R))

length_beam = m.sqrt(R**2+(Ca-R)**2)


Iyy_straight = ((tsk*length_beam**(3)*m.cos(angle)*m.cos(angle))/12+A_sk2*((Ca-R)/2-zc)**2) #2 beams
print("Moment of inertia of straight skin parts Iyy is:", Iyy_straight, "m^4")

Izz_straight = ((tsk*length_beam**(3)*m.sin(angle)*m.sin(angle))/12+A_sk2*(R/2)**2) #2 beams
print("Moment of inertia of straight skin parts Izz is:", Izz_straight, "m^4")

#Calculate moment of inertia of arc skin part
Izz_arc = 1/2*m.pi*tsk*R**3         #
Iyy_arc = 1/2*m.pi*tsk*R**3-4*R**3*tsk/m.pi+A_sk1*(2*R/m.pi+zc)**2
print("Moment of inertia of arc skin part Iyy is:", Iyy_arc, "m^4")
print("Moment of inertia of arc skin part Izz is:", Izz_arc, "m^4")

#Calculate moment of inertia of spar
Izz_spar = (tsp*ha**3)/12
Iyy_spar = A_sp*(zc)**2

print("Moment of inertia of spar Iyy is:", Iyy_spar,"m^4")
print("Moment of inertia of spar Izz is:", Izz_spar, "m^4")

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

print("Moment of inertia of stiffners Izz is:", Izz_st, "m^4") 
print("Moment of inertia of stiffners Iyy is:", Iyy_st, "m^4")
#Calculate total moment of inertia 
Iyy_total = 2*Iyy_straight + Iyy_arc + Iyy_spar + Iyy_st
Izz_total = 2*Izz_straight + Izz_arc + Izz_spar + Izz_st

print("Total moment of inertia Iyy =", Iyy_total, "m^4")
print("Total moment of inertia Izz =", Izz_total, "m^4") 

"Torsional constant calculation"