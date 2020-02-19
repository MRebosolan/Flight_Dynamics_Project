# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:23:42 2020

@author: Nils
"""
import numpy as np
import math


def delta_q_i_base(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, boom_area_i, y_prime, z_prime):
    delta_q_i= - shear_force_y*boom_area_i*y_prime/MOI_z_prime - shear_force_z*boom_area_i*z_prime/MOI_y_prime
    return delta_q_i

#----------------------------------------------------------------------------------------------------------------
#defintion of areas:
#top1===> the y positive part of the circular section of the aileron
#top2==> the y positive part of the inclined section of the aileron
#bottom1==>the y negative part of the circular part of the aileron
#bottom2==>the y negative part of the inclined section of the aileron
#spar==> spar of the aileron
#sparA==> positive y part of spar
#sparB==> negative y part of spar

def alternative_q_base_top1(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, integral_step, aileron_height, skin_thickness,z_pos,y_pos):
    #first part a: determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    radius=aileron_height/2
    s=np.arctan(y_pos/z_pos)*radius
    s_1=0
    while s_1<=s:
        y=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius*(np.arctan(s_1/radius))
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius*(np.arctan(s_1/radius))
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime - shear_force_z*integral_value2/MOI_y_prime
    
    return q_base_i

def alternative_q_base_bottom1(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, integral_step, aileron_height, skin_thickness,z_pos,y_pos):
    #first part a: determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    radius=aileron_height/2
    s=(np.pi*0.5-np.arctan(y_pos/z_pos))*radius
    s_1=0
    while s_1<=s:
        y=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius*(np.arctan(s_1/radius))
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius*(np.arctan(s_1/radius))
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime - shear_force_z*integral_value2/MOI_y_prime
    
    return q_base_i


def alternative_q_base_top2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,z_pos,y_pos):
    #first part a: determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=math.sqrt((aileron_height*0.5-y_pos)**2+z_pos**2)
    theta=np.arctan(aileron_height/(chord_length-0.5*aileron_height))
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1*np.sin(theta)
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=-1*s_1*np.cos(theta)
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime - shear_force_z*integral_value2/MOI_y_prime
    
    return q_base_i

def alternative_q_base_bottom2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,z_pos,y_pos):
    #first part a: determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=math.sqrt((y_pos)**2+(chord_length-aileron_height*0.5-z_pos)**2)
    theta=np.arctan(aileron_height/(chord_length-0.5*aileron_height))
    s_1=0
    while s_1<=s:
        y=s_1*np.sin(theta)
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=chord_length-aileron_height*0.5-s_1*np.cos(theta)
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime - shear_force_z*integral_value2/MOI_y_prime
    
    return q_base_i

def alternative_q_base_spar(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,y_pos):
    #first part : determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=aileron_height*0.5-y_pos
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime 
    
    return q_base_i

def alternative_q_base_sparA(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,y_pos):
    #first part : determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=aileron_height*0.5-y_pos
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime 
    
    return q_base_i

def alternative_q_base_sparB(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,y_pos):
    #first part : determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=-aileron_height*0.5-y_pos
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_i= - shear_force_y*integral_value1/MOI_z_prime 
    
    return q_base_i

def total_q_i(base_shear,torque_shear, zero_shear):
    total_q_i=base_shear+torque_shear+zero_shear
    return total_q_i

def relation_shear_1_and_2_torque_and_torque(Overall_torque, aileron_height, skin_thickness, spar_thickness, chord_length):
    
    #The first out output is the torque shear distribution, the second is the second torque shear ditribution and the third is the ratio between both shear distributions both for the shear due to torque and the residual shear.
    
    t=True

    length_after_spar=(chord_length-aileron_height*0.5)
    radius=aileron_height/2
    A_1=(np.pi)**2*aileron_height/2 #cross section area 1 
    A_2=aileron_height*(length_after_spar) #cros section area 2
    
    K_1= radius*np.pi/(skin_thickness)+2*aileron_height/spar_thickness
    K_2= (2/spar_thickness)*math.sqrt(length_after_spar**2+radius**2)+aileron_height/spar_thickness
    K_f=K_2/K_1
    
    shear_2=0.0
    while t is True:
        shear_1_A=(Overall_torque-2*A_2*shear_2)/(2*A_1)
        shear_1_B=K_f*shear_2

        shear_2+=0.0001
        
        if round(shear_1_A,4)==round(shear_1_B,4):
            t=False
            
    return shear_1_A/K_f, shear_1_A, K_f

def rate_twist_at_x(shear_torque_1,shear_torque_2, shear_zero_1,shear_zero_2, aileron_height, spar_thickness, skin_thickness):
    rate_twist_torque = shear_torque_1*np.pi*aileron_height*0.5/skin_thickness+(shear_torque_1-shear_torque_2)*aileron_height/spar_thickness
    rate_twist_zero=shear_zero_1*np.pi*aileron_height*0.5/skin_thickness+(shear_zero_1-shear_zero_2)*aileron_height/spar_thickness
    
    rate_twist_x=rate_twist_zero+rate_twist_torque
    
    return rate_twist_x

def twist(x_set_of_positions, rate_twist_at_x):#the x set of positions has to start at zero
    #the x_set_of_positions starts at the origin, that is at hinge 1. thus, hinge 1 is fixed in twist
    x_set_of_twist=[]
    twist=0
    i=1
    while i<=len(x_set_of_positions):
        twist+=(x_set_of_positions[i]-x_set_of_positions[i-1])*rate_twist_at_x #calculates the twist at each point
        x_set_of_twist.append(twist)
        i+=1
    return twist, x_set_of_twist #the first output is the twist at the edge of the aileron
        
    
def deflection_due_to_torque(twist, pos_z, pos_y, shear_center_y, shear_center_z):
    phi=np.arctan(pos_y/pos_z)
    distance=math.sqrt(pos_y**2+pos_z**2)
    delta_deflec=twist*distance
    deflection_due_to_torque_y=delta_deflec*np.sin(phi)
    deflection_due_to_torque_z=delta_deflec*np.cos(phi)
    
    return deflection_due_to_torque_y, deflection_due_to_torque_z

def deflection_y(moment_z_set,x_set_of_positions):
    
    #finding the first step of the integral
    i=0
    while i<=len(x_set_of_positions):
        
        
    
    

    





        
    
    
    
    
    
    
    
    
    
    
    
    