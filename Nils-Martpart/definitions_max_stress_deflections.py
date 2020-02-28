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

def alternative_q_base_top1(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, integral_step, aileron_height, skin_thickness,z_pos,y_pos,boom_list):
    #y_pos and z_pos are the positions of the maximum point of s for this part
    #first part a: determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    radius=aileron_height/2
    if z_pos==0:
        s=np.pi*0.5
    else:
        s=np.arctan(y_pos/z_pos)*radius
    s_1=0
    while s_1<=s:
        y=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius*(np.arctan(s_1/radius))
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            y=math.sqrt(1/(1+(np.arctan(s_1_values[i]/radius))**2))*radius*(np.arctan(s_1_values[i]/radius))
            if y>=boom_list[u][1]:
                extra+=boom_list[u][0]*boom_list[u][1]
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1+extra)
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    integral_values2=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            z=math.sqrt(1/(1+(np.arctan(s_1_values[i]/radius))**2))*radius
            if z<=boom_list[u][2] and boom_list[u][2]>0:
                extra+=boom_list[u][0]*boom_list[u][2]
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values2.append(integral_value2)
        i+=1  

    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values2)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime - shear_force_z*integral_values2[i]/MOI_y_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values

def alternative_q_base_bottom1(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, integral_step, aileron_height, skin_thickness,z_pos,y_pos, boom_list):
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
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)-1):
            y=math.sqrt(1/(1+(np.arctan(s_1_values[i]/radius))**2))*radius*(np.arctan(s_1_values[i]/radius))
            if y>=boom_list[u][1]:
                extra+=boom_list[u][0]*boom_list[u][1]
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1)
        i+=1
        
    #first part b: determining the integral of thickness and position with respect to z
    s_1_values=[]
    z_values=[]
    s_1=0
    while s_1<=s:
        z=math.sqrt(1/(1+(np.arctan(s_1/radius))**2))*radius
        z_values.append(z*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value2=0
    integral_values2=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            z=math.sqrt(1/(1+(np.arctan(s_1_values[i]/radius))**2))*radius
            if z>=boom_list[u][2] and boom_list[u][2]>0:
                extra+=boom_list[u][0]*boom_list[u][2]
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values2.append(integral_value2)
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values2)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime - shear_force_z*integral_values2[i]/MOI_y_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values


def alternative_q_base_top2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,z_pos,y_pos, boom_list):
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
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            y=aileron_height*0.5-s_1_values[i]*np.sin(theta)
            if y<=boom_list[u][1] and boom_list[u][2]<0:
                extra+=boom_list[u][0]*boom_list[u][1]
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1+extra)
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
    integral_values2=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            z=-1*s_1_values[i]*np.cos(theta)
            if z<=boom_list[u][2] and boom_list[u][2]<0:
                extra+=boom_list[u][0]*boom_list[u][2]
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values2.append(integral_value2+extra)
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values2)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime - shear_force_z*integral_values2[i]/MOI_y_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values

def alternative_q_base_bottom2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,z_pos,y_pos, boom_list):
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
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            y=s_1_values[i]*np.sin(theta)
            if y<=boom_list[u][1] and boom_list[u][2]<0 :
                extra+=boom_list[u][0]*boom_list[u][1]
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1+extra)
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
    integral_values2=[]
    i=1
    while i <=(len(s_1_values)-1):
        extra=0
        for u in range(len(boom_list)):
            z=chord_length-aileron_height*0.5-s_1_values[i]*np.cos(theta)
            if z>=boom_list[u][2] and boom_list[u][2]<0 :
                extra+=boom_list[u][0]*boom_list[u][2]
        integral_value2+=(s_1_values[i]-s_1_values[i-1])*z_values[i-1]+(z_values[i]-z_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values2.append(integral_value2+extra)
        i+=1  
    
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values2)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime - shear_force_z*integral_values2[i]/MOI_y_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values

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
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1)
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values1)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values

def alternative_q_base_sparA(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,y_pos):
    #first part : determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=y_pos
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    
    integral_value1=0
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1)
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values1)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime
        q_base_list.append(q_base_i)
        i+=1
    
    return q_base_list, s_1_values

def alternative_q_base_sparB(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness,y_pos):
    #first part : determining the integral of thickness and position with respect to y
    s_1_values=[]
    y_values=[]
    s=aileron_height*0.5+y_pos
    s_1=0
    while s_1<=s:
        y=aileron_height*0.5-s_1
        y_values.append(y*skin_thickness)
        s_1_values.append(s_1)
        s_1+=integral_step
    
    integral_value1=0
    integral_values1=[]
    i=1
    while i <=(len(s_1_values)-1):
        integral_value1+=(s_1_values[i]-s_1_values[i-1])*y_values[i-1]+(y_values[i]-y_values[i-1])*(s_1_values[i]-s_1_values[i-1])
        integral_values1.append(integral_value1)
        i+=1
          
    #second part: determing the actual q_base for the point with position y and z
    q_base_list=[]
    i=0
    while i<=(len(integral_values1)-1):
        q_base_i= - shear_force_y*integral_values1[i]/MOI_z_prime
        q_base_list.append(q_base_i)
        i+=1 
    
    return q_base_list, s_1_values

def total_q_i(base_shear,torque_shear, zero_shear):
    total_q_i=base_shear+torque_shear+zero_shear
    return total_q_i

def relation_shear_1_and_2_torque_and_torque(Overall_torque, aileron_height, skin_thickness, spar_thickness, chord_length):
    
    #The first out output is the torque shear distribution, the second is the second torque shear ditribution and the third is the ratio between both shear distributions both for the shear due to torque and the residual shear.
    
    t=True
    length_after_spar=(chord_length-aileron_height*0.5)
    radius=aileron_height/2
    A_1=((np.pi)**2*(aileron_height/2)**2)/2 #cross section area 1 
    A_2=radius*(length_after_spar) #cros section area 2
    
    K_1= (radius*np.pi/(skin_thickness)+aileron_height/spar_thickness)*A_2/A_1 + aileron_height/spar_thickness
    K_2= (2/skin_thickness)*math.sqrt(length_after_spar**2+radius**2)+(aileron_height/spar_thickness)*(1+A_2/A_1)
    K_f=K_2/K_1
    
    shear_2= Overall_torque/(2*A_1*K_f+2*A_2)
    shear_1 =shear_2*K_f

    return shear_2, shear_1, K_f

def deflection_z_bending_stress(E, I_yy, moment_y_set,x_set_of_positions,MOI_y_prime, E_modulus):
    #deflection=int_numerical+D*x+C
    #finding the first step of the integral
    #1A=-K*dx
    #C4=36.54718612
    #C5=29.79453841
    i=0
    K_list=[]
    while i<=(len(x_set_of_positions)-1):
        K=-moment_y_set[i]/(MOI_y_prime*E_modulus)
        K_list.append(K)
        i+=1

    integral_value1=0
    integral_values1=[]
    i=1
    while i <=(len(x_set_of_positions)-1):
        integral_value1+=(x_set_of_positions[i]-x_set_of_positions[i-1])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i]-x_set_of_positions[i-1])
        integral_values1.append(integral_value1)
        i+=1
        
    integral_values1.insert(0,0)
    
    #x_set_of_positions=x_set_of_positions[1:]
    #finding the second step of the integral
    i=0
    K_list=[]
    while i<=(len(integral_values1)-1):
        K=integral_values1[i]
        K_list.append(K)
        i+=1
        
    integral_value2=0
    i=1
    integral_values2=[]
    while i <=(len(integral_values1)-1):
        if i==1:
            integral_value2+=(x_set_of_positions[i])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i])*0.5
        else:
            integral_value2+=(x_set_of_positions[i]-x_set_of_positions[i-1])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i]- x_set_of_positions[i-1])*0.5
        integral_values2.append(integral_value2)
        i+=1
        
    integral_values2.insert(0,0)   
    return integral_values2

def deflection_y_bending_stress(E, I_zz, moment_z_set,x_set_of_positions,MOI_z_prime, E_modulus):
    #deflection=int_numerical+D*x+C
    #finding the first step of the integral
    #1A=-K*dx
    #C1=6.67973948
    #C2=-0.3474889
    i=0
    K_list=[]
    while i<=(len(x_set_of_positions)-1):
        K=-moment_z_set[i]/(MOI_z_prime*E_modulus)
        K_list.append(K)
        i+=1
    
     
    integral_value1=0
    integral_values1=[]
    i=1
    while i <=(len(x_set_of_positions)-1):
        integral_value1+=(x_set_of_positions[i]-x_set_of_positions[i-1])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i]-x_set_of_positions[i-1])*0.5
        integral_values1.append(integral_value1)
        i+=1
    integral_values1.insert(0,0)
    
    #x_set_of_positions=x_set_of_positions[1:]
    #finding the second step of the integral
    i=0
    K_list=[]
    while i<=(len(integral_values1)-1):
        K=integral_values1[i]
        K_list.append(K)
        i+=1
     
    integral_value2=0
    i=1
    integral_values2=[]
    while i <=(len(integral_values1)-1):
        if i==1:
            integral_value2+=(x_set_of_positions[i])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i])*0.5
        else:
            integral_value2+=(x_set_of_positions[i]-x_set_of_positions[i-1])*K_list[i-1]+(K_list[i]-K_list[i-1])*(x_set_of_positions[i]-x_set_of_positions[i-1])

        integral_values2.append(integral_value2)
        i+=1
    integral_values2.insert(0,0)
    
    return integral_values2


def rate_twist_at_x(shear_torque_1,shear_torque_2, shear_zero_1,shear_zero_2, aileron_height, spar_thickness, skin_thickness, shear_modulus):
    rate_twist_torque = (shear_torque_1*np.pi*aileron_height*0.5/skin_thickness+(shear_torque_1-shear_torque_2)*aileron_height/spar_thickness)/shear_modulus
    rate_twist_zero=(shear_zero_1*np.pi*aileron_height*0.5/skin_thickness+(shear_zero_1-shear_zero_2)*aileron_height/spar_thickness)/shear_modulus
    
    rate_twist_x=rate_twist_zero+rate_twist_torque
    #print(rate_twist_x)
    return rate_twist_x

def twist(G,J, x_set_of_positions, rate_twist_at_x):#the x set of positions has to start at zero
    #the x_set_of_positions starts at the origin, that is at hinge 1. thus, hinge 1 is fixed in twist
    #C3=-15.2202193
    x_set_of_twist=[]
    x_set=[]
    #x_set_of_positions=x_set_of_positions[1:]
    twist=0
    i=1
    while i<=len(x_set_of_positions)-1:
        if i==1:
            twist+=(x_set_of_positions[i]-x_set_of_positions[i-1])*(rate_twist_at_x[i]) #calculates the twist at each point
        else:
            twist+=(x_set_of_positions[i]-x_set_of_positions[i-1])*(rate_twist_at_x[i] )#calculates the twist at each point

        x_set_of_twist.append(twist)
        x_set.append(x_set_of_positions[i])
        i+=1
    x_set_of_twist.insert(0,0)

    '''
    for i in range(len(x_set_of_twist)):
        x_set_of_twist[i]+=-C3/(G*J)

    '''   
    return twist, x_set_of_twist, x_set#the first output is the twist at the edge of the aileron
        
    
def deflection_due_to_torque_and_bending(E, I_zz, I_yy, x_set_of_twist, x_set_of_positions, shear_center_y, shear_center_z, deflection_y_bending_set, deflection_z_bending_set, x_location_hinge1, x_location_hinge3, deflection_hinge_1, deflection_hinge_3 ):
    lst_Deflections=[]

    for i in range(len(x_set_of_twist)):
        distance=math.sqrt(shear_center_y**2+shear_center_z**2)
        delta_deflec_y_torque=np.sin(x_set_of_twist[i])*distance
        delta_deflec_z_torque=distance*(1-np.cos(x_set_of_twist[i]))
        deflection_y=delta_deflec_y_torque+deflection_y_bending_set[i]
        deflection_z=delta_deflec_z_torque+deflection_z_bending_set[i]
        lst_Deflections.append([x_set_of_positions[i], deflection_y, deflection_z])
        
    
    index1=list(x_set_of_positions).index(round(x_location_hinge1,2))
    index3=list(x_set_of_positions).index(round(x_location_hinge3,2))

    diff_x3_x1=round(x_location_hinge3,2)-round(x_location_hinge1,2)
    #using the boundary conditions to find the integration constants, y deflections
    D1=((-1*deflection_hinge_1+lst_Deflections[index1][1]-lst_Deflections[index3][1]+deflection_hinge_3))/(diff_x3_x1)
    C1=(deflection_hinge_1-lst_Deflections[index1][1]-D1*round(x_location_hinge1,2))
    
    #using the boundary conditions to find the integration constants, z deflections
    D_z1=((lst_Deflections[index1][2]-lst_Deflections[index3][2]))/(diff_x3_x1)
    C_z1=(-lst_Deflections[index1][2]-D_z1*round(x_location_hinge1,2))


    #correcting the data using the constants
    i=0
    lst_final_Deflections=[]
    while i<=(len(x_set_of_positions)-1):
        y_deflect=lst_Deflections[i][1]+(D1*x_set_of_positions[i]+C1)
        z_deflect=lst_Deflections[i][2]+(D_z1*x_set_of_positions[i]+C_z1)
                
        lst_final_Deflections.append([x_set_of_positions[i], y_deflect, z_deflect])
        
        i+=1
    
        
    #print(lst_Deflections[index3][1], D, x_set_of_positions[index3], C)
    #print(lst_Deflections[index3][2], D_z, x_set_of_positions[index3], C_z)
    
    return lst_final_Deflections








        


    
        
    
    

    





        
    
    
    
    
    
    
    
    
    
    
    
    
