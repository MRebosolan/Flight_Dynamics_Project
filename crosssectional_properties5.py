"""
Created on Mon Feb 17 12:32:42 2020
@author: lizzy & richelle
"""
import numpy as np
import math as m
import matplotlib.pyplot as plt

aircraft = "A320"
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
Ra=ha/2          #radius 
length=m.pi*Ra**2

Ab=wst*tst+(hst-tst)*tst #m
l=2*m.sqrt((Ca-Ra)**2+(Ra)**2)+m.pi*Ra
b=l/nst
k=b/Ra #radians

z_coordinates=[]
y_coordinates=[]
z_spar=[]
y_spar=[]
A_st=(hst)*tst+tst*wst

"Stringer placement"
z_stringer=[] 
y_stringer=[]
stringer_array=[]
z_stringer.append(-Ra)
y_stringer.append(0)
stringer_array.append([A_st, 0,-Ra])
h=0
for n in range(0,10001):
    z=-Ra*m.cos(n*b/(Ra))
    y=Ra*m.sin(n*b/(Ra))
    if z<=0. and z>-Ra+0.01:
        z_stringer.append(z)
        z_stringer.append(z)
        y_stringer.append(y)
        y_stringer.append(-y)
        stringer_array.append([A_st, y,z])
        stringer_array.append([A_st, -y,z])
        h=h+2
    elif z>0 or nst==1:
        break
for n in range(0,10001):
     last_room=0.5*m.pi*Ra-(h/2)*b
     z1=0
     z2=(n-2)*(Ca-Ra)/10000
     y1=Ra
     y2=-Ra/(Ca-Ra)*z2+Ra
     l_p=round(m.sqrt((z1-z2)**2+(y1-y2)**2),7)
     z=z2
     y=y2
     for i in range (0,nst):
         if abs(l_p-round((i*b-last_room),7))<=0.000022:
             z_stringer.append(z)
             z_stringer.append(z)
             y_stringer.append(y)
             y_stringer.append(-y)
             stringer_array.append([A_st, y,z])
             stringer_array.append([A_st, -y,z])
         elif z>(Ca-Ra) or nst==1:
                 break
   
stringers=np.array(stringer_array)
print(stringers) 
"airfoil placement"
for i in range(0,101):
    z=-Ra+Ra*i/(100)
    y=np.sqrt(Ra**2-z**2)
    z_coordinates.append(z)
    y_coordinates.append(y)
for i in range(-100,101):
    z=0
    y=i/100*Ra
    z_spar.append(z)
    y_spar.append(y)
for i in range(0,101):
    z=i*(Ca-Ra)/100
    y=-Ra/(Ca-Ra)*z+Ra
    z_coordinates.append(z)
    y_coordinates.append(y)

for i in range(0,101):
    z=(100-i)*(Ca-Ra)/100
    y=-Ra/(Ca-Ra)*z+Ra
    z_coordinates.append(z)
    y_coordinates.append(-y)

for i in range(0,101):
    z=-Ra+Ra*(100-i)/100
    y=np.sqrt(Ra**2-z**2)
    z_coordinates.append(z)
    y_coordinates.append(-y)

"Calculation Centroid"
A_sk1=m.pi*Ra*tsk
A_sk2=m.sqrt((Ca-Ra)**2+Ra**2)*tsk
A_sp=ha*tsp

A_tot=A_sk1+2*A_sk2+A_sp+nst*A_st
print("Total area:", A_tot, "m^2")

Az_st=A_st*sum(z_stringer) #m^2
Az_sk1=A_sk1*2*-Ra/(m.pi)  #m^2
Az_sk2=A_sk2*1/2*(Ca-Ra)   #m^2
Ay_st=A_st*sum(y_stringer)
Ay_sk1=A_sk1*0
Ay_sk2=A_sk2*1/2*(Ra)
Ay_sk3=A_sk2*-1/2*(Ra)
zc=(Az_st+Az_sk1+2*Az_sk2)/A_tot
yc=(Ay_st+Ay_sk1+Ay_sk2+Ay_sk3)/A_tot          #due to symmetry it should be 0
print("centroid",zc,yc)

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
angle = m.atan(Ra/(Ca-Ra)) #radians
length_beam = m.sqrt(Ra**2+(Ca-Ra)**2) #m
Iyy_straight = ((tsk*length_beam**(3)*m.cos(angle)*m.cos(angle))/12+A_sk2*((Ca-Ra)/2-zc)**2) #2 beams
#print("Moment of inertia of straight skin parts Iyy is:", Iyy_straight, "m^4")

Izz_straight = ((tsk*length_beam**(3)*m.sin(angle)*m.sin(angle))/12+A_sk2*(Ra/2)**2) #2 beams
#print("Moment of inertia of straight skin parts Izz is:", Izz_straight, "m^4")

#Calculate moment of inertia of arc skin part
Izz_arc = 1/2*m.pi*tsk*Ra**3         #
Iyy_arc = 1/2*m.pi*tsk*Ra**3-4/m.pi*Ra**3*tsk+A_sk1*(2*Ra/m.pi+zc)**2
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
J1=4*(0.5*m.pi*Ra**2)**2/((m.pi*Ra/tsk))
J2=4*(Ra*(Ca-Ra))**2/((2*m.sqrt(Ra**2+(Ca-Ra)**2)/tsk))

J=J1+J2
print("Torsional constant", J)