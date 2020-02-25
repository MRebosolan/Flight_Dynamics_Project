# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 22:21:33 2020

@author: Mart
"""

import numpy as np

"This first part gives values that have to be inputted an can be changed"

theta =  #input here the angle between coordinate systems in radians
E =  #input here the E-modulus of the material
Izz =  #input here the moment of inertia about the z-prime axis
Iyy =  #input here the moment of inertia about the y-prima axis
eta = #input here the distance between shear centre and hinge line (must be a positive number!!!!)
G =  #input here the G-constant of the material
J =  #input here the torsional stiffnes constant of the structure
x1 =  #input here the x-position of hinge 1
x2 =  #input here the x-position of hinge 2
x3 =  #input here the x-position of hinge 3
xa =  #input what you defined as xa in your moment thingeys
ha =  #input what you defined as ha in your moment thingeys
xp =  #input here the x-position of actuator P
P =  #input the loading P
ksi = #input here the distance between the point where P acts and the shear centre


A = np.array([[R1y, R2y, R3y, R1z, R2z, R3z, RI, C1, C2, C3, C4, C5],\
             [0, 0, 0, 0, 0, 0, 0, -x1/E/Izz, -1/E/Izz, eta/G/J, 0, 0], \
             [np.cos(theta)*(x2-x1)**3/6/E/Izz+np.cos(theta)*eta/G/J*eta*(x2-x1), 0, 0, np.sin(theta)*(x2-x1)**3/6/E/Izz+eta*np.sin(theta)/G/J*eta*(x2-x1), 0, 0, (xa/2)**3*np.sin(theta)/E/Izz/6+eta*np.cos(theta)/G/J*eta*xa/2 + eta*np.cos(theta)/G/J*(ha/2+eta*np.sin(theta))*xa/2, -x2/E/Izz, -1/E/Izz, eta/G/J, 0, 0],\
             [np.cos(theta)*(x3-x1)**3/6/E/Izz+np.cos(theta)*eta/G/J*eta*(x3-x1), np.cos(theta)*(x3-x2)**3/6/E/Izz+np.cos(theta)*eta/G/J*eta*(x3-x2), 0, np.sin(theta)*(x3-x1)**3/6/E/Izz+np.sin(theta)*eta/G/J*eta*(x3-x1), np.sin(theta)*(x3-x2)**3/6/E/Izz+np.sin(theta)*eta/G/J*eta*(x3-x2), 0, (x3-x2+xa/2)**3*np.sin(theta)/6/E/Izz+np.cos(theta)*eta/G/J*(x3-x2+xa/2)+np.cos(theta)*eta/G/J*(ha/2+eta*np.sin(theta))*(x3-x2+xa/2), -x3/E/Izz, -1/E/Izz, eta/G/J, 0, 0],\
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x1/E/Iyy, -1/E/Iyy],\
             [np.sin(theta)/6/E/Iyy*(x2-x1)**3, 0, 0, np.cos(theta)/6/E/Iyy*(x2-x1)**3, 0, 0, np.cos(theta)/6/E/Iyy*(xa/2)**3, 0, 0, 0, -x2/E/Iyy, -1/E/Iyy],\
             [np.sin(theta)/6/E/Iyy*(x3-x1)**3, np.sin(theta)/6/E/Iyy*(x3-x2)**3, 0, np.cos(theta)/6/E/Iyy*(x3-x1)**3, np.cos(theta)/6/E/Iyy*(x3-x2)**3, 0, np.cos(theta)/6/E/Iyy*(x3-x2+xa/2)**3, 0, 0, 0, -x3/E/Iyy, -1/E/Iyy],\
             [np.sin(theta)/E/Iyy/6*(xp-x1)**3+ksi*np.sin(theta)*np.cos(theta)*eta/G/J*(xp-x1), np.sin(theta)/E/Iyy/6*(xa/2)**3++ksi*np.sin(theta)*np.cos(theta)*eta/G/J*(xp-x2), 0, np.cos(theta)/E/Iyy/6*(xp-x1)**3+ksi*np.sin(theta)/G/J*np.sin(theta)*eta*(xp-x1), np.cos(theta)/E/Iyy/6*(xa/2)**3+ksi*np.sin(theta)/G/J*np.sin(theta)*eta*(xp-x2), 0, np.cos(theta)/E/Iyy/6*(xa)**3+ksi*np.sin(theta)/G/J*np.cos(theta)*eta*(xa)+ksi*np.sin(theta)/G/J*np.cos(theta)*(ha/2+eta*np.sin(theta))*(xa), 0, 0, np.sin(theta)*ksi/G/J, -xp/E/Iyy, -1/E/Iyy],\
                                                                                                                               ]                                                                                           )

  
"Constant terms:"
"v(h1)" 0.01103*np.cos(theta)+1/E/Izz*(TRIPLEINTEGRALOFLOADINGQYUPTOHINGE1)-(DOUBLEINTEGRALOFLOADINTAUGUPTOHINGE1)*eta/G/J
"v(h2)" (TRIPLEINTEGRALOFLOADINGQYUPTOHINGE2)/E/Izz - eta/G/J*(DOUBLEINTEGRALOFLOADINGTAUUPTOHINGE2)
"v(h3)" 0.01642*np.cos(theta)+(TRIPLEINTEGRALOFLOADINGQYUPTOHINGE3)/E/Izz-eta?G/J*(DOUBLEINTEGRALOFLOADINGTAUUPTOHINGE3)+P/6/E/Izz*(x3-x2-xa/2)**3*np.sin(theta)+P*eta/G/J*np.cos(theta)*(ha/2+eta*np.sin(theta))*(x3-x2-xa/2)+eta*P*np.sin(theta)*eta/G/J*(x3-x2-xa/2)
"w(h1)" -0.01103*np.sin(theta)-(TRIPLEINTEGRALOFLOADINGQZUPTOHINGE1)/E/Iyy
"w(h2)" -(TRIPLEINTEGRALOFLOADINGQZUPTOHINGE2)/E/Iyy
"w(h3)" -0.01642*np.sin(theta)-(TRIPLEINTEGRALOVERLOADINGQZUPTOHINGE3)/E/Iyy+P*np.cos(theta)/E/Iyy/6*(x3-x2-xa/2)**3
"w(P)" -(TRIPLEINTEGRALOFLOADINGQZUPTOACTUATORP)/E/Iyy-ksi*np.sin(theta)/G/J*(DOUBLEINTEGRALOFLOADINGTAUUPTOACTUATORP)

