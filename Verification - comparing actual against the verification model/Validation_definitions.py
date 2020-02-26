# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:25:00 2020

@author: Nils
"""

import numpy as np
import math
import re

#defintion of members:
#top1===> the y positive part of the circular section of the aileron [1]
#top2==> the y positive part of the inclined section of the aileron [2]
#bottom1==>the y negative part of the circular part of the aileron [3]
#bottom2==>the y negative part of the inclined section of the aileron [4]
#spar==> spar of the aileron [5]

def deflection_validation(deflection_set_actual, deflection_set_validation_data):
    error_list=[]
    for i in range(len(deflection_set_actual)):
        if deflection_set_validation_data[i][3]==0:
            y_error=-1
        else:
            y_error=abs(deflection_set_actual[i][1]-deflection_set_validation_data[i][3])
        
        if deflection_set_validation_data[i][4]==0:
            z_error=-1
        else:
            z_error=abs(deflection_set_actual[i][2]-deflection_set_validation_data[i][4])
        error_list.append([deflection_set_actual[i][0], y_error, z_error])
    return error_list

#the original comparison was with percentage but the actual difference is way more telling
def twist_validation(twist_set_actual, twist_set_validation_data):
    error_list=[]
    for i in range(len(twist_set_actual)):
        twist_error=abs(twist_set_actual[i]-twist_set_validation_data[i])
        error_list.append(twist_error)
    
    twists_actual=[]
    twists_validation=[]
    for i in range(len(twist_set_actual)):
        twists_actual.append(twist_set_actual[i])
        twists_validation.append(twist_set_validation_data[i])
        
    maximum_minimum_difference_actual=max(twist_set_actual)-min(twist_set_actual)
    maximum_minimum_difference_validation=max(twist_set_validation_data)-min(twist_set_validation_data)
    
    return error_list, maximum_minimum_difference_actual, maximum_minimum_difference_validation

def max_stress_validation(max_stress_actual, max_stress_validation,  actual_member_max_stress, validation_member_max_stress):
    error_max_stress=abs(max_stress_actual-max_stress_validation)*100/max_stress_validation
    
    if actual_member_max_stress==validation_member_max_stress:
        switch=1
    else:
        switch=0
        
    #diff_x=abs(x_location_actual_percent_span-x_location_validation_percent_span)*100/x_location_validation_percent_span
        
    #the switch which is returned labels wether the maximum stress is located on the same member of the aileron cross section or not
    #these are defined at the beginning of the program of defintions for validation
    return error_max_stress, switch

def reading_inputs(file_name):
    #this is not for Shear Center
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)
        
    lst2=[]
    for i in range(9,6597):
        lst2.append(lst[i])
    
    
    lst3=[]
    for k in range(len(lst2)):
        lst4=[]
        for i in range(len(lst2[k])):
            if lst2[k][i]!='':
                lst4.append(float(lst2[k][i]))
        lst3.append(lst4)
    file.close()
    
    lst4=[]
    for i in range(len(lst3)):
        if lst3[i][2]==0 and lst3[i][3]==0:
            lst4.append(lst3[i])
    
    lst4=sorted(lst4, key = lambda x: x[2])
    
    #lst 4 gives an output list of each line [node number, x, y, z]. The outputs are only the nodes which are on the hinge line
    return lst4, lst3

def reading_outputs_deflections_bent(file_name, interest_nodes):
    
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)

    lst2=[]
    for i in range(26724,33312):
        lst2.append(lst[i])
    
    
    lst3=[]
    for k in range(len(lst2)):
        lst4=[]
        for i in range(len(lst2[k])):
            if lst2[k][i]!='':
                lst4.append(float(lst2[k][i]))
        lst3.append(lst4)
    file.close()
    
    lst4=[]
    for i in range(len(lst3)):
        for j in range(len(interest_nodes)):
            if lst3[i][0]==interest_nodes[j][0]:
                lst4.append([lst3[i][0],lst3[i][1], lst3[i][2]/1000, lst3[i][3]/1000, lst3[i][4]/1000])
    lst4=sorted(lst4, key = lambda x: x[2])
    return len(lst4), lst4

def reading_outputs_deflections_unbent(file_name, interest_nodes):
    
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)

    lst2=[]
    for i in range(33374,39962):
        lst2.append(lst[i])
    
    
    lst3=[]
    for k in range(len(lst2)):
        lst4=[]
        for i in range(len(lst2[k])):
            if lst2[k][i]!='':
                lst4.append(float(lst2[k][i]))
        lst3.append(lst4)
    file.close()
    
    lst4=[]
    for i in range(len(lst3)):
        for j in range(len(interest_nodes)):
            if lst3[i][0]==interest_nodes[j][0]:
                lst4.append([lst3[i][0],lst3[i][1], lst3[i][2]/1000, lst3[i][3]/1000, lst3[i][4]/1000])
    lst4=sorted(lst4, key = lambda x: x[2])           
    return len(lst4), lst4

def reading_outputs_deflections_pure_bending(file_name, interest_nodes):
    
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)

    lst2=[]
    for i in range(20075,26662):
        lst2.append(lst[i])
    
    
    lst3=[]
    for k in range(len(lst2)):
        lst4=[]
        for i in range(len(lst2[k])):
            if lst2[k][i]!='':
                lst4.append(float(lst2[k][i]))
        lst3.append(lst4)
    file.close()
    
    lst4=[]
    for i in range(len(lst3)):
        for j in range(len(interest_nodes)):
            if lst3[i][0]==interest_nodes[j][0]:
                lst4.append([lst3[i][0],lst3[i][1], lst3[i][2]/1000, lst3[i][3]/1000, lst3[i][4]/1000])
    lst4=sorted(lst4, key = lambda x: x[2])           
    return len(lst4), lst4

def average_y_z_deflections_due_to_torque_validated_data(deflection_y_z_unbent, deflection_y_z_bent, deflection_y_z_bending):
    lst=[]
    for i in range(len(deflection_y_z_unbent)):
        D=(deflection_y_z_unbent[i][3]+(deflection_y_z_bent[i][3]-deflection_y_z_bending[i][3]))/2
        C=(deflection_y_z_unbent[i][4]+(deflection_y_z_bent[i][4]-deflection_y_z_bending[i][4]))/2
        E=math.sqrt(C**2+D**2)
        lst.append([i, deflection_y_z_unbent[i][0],D,C,E])
    
    return lst

def computing_outputs_twist(deflections_y_and_z_hinge_lines, shear_center_y, shear_center_z):
    #these deflection in y and z are calculated by take the average of the unbent case and the dif. between the bent and pure bending case
    lst_final_values=[]
    for i in range(len(deflections_y_and_z_hinge_lines)):
        twist=np.arctan((deflections_y_and_z_hinge_lines[i][2]-shear_center_y)/(shear_center_z-deflections_y_and_z_hinge_lines[i][3]))
        lst_final_values.append([deflections_y_and_z_hinge_lines[i][0], deflections_y_and_z_hinge_lines[i][1], twist])
    return lst_final_values    

def reading_max_stress(file_name, node):
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)
    
    lst2=[]
    for i in range(9,6597):
        lst2.append(lst[i])

    
    lst3=[]
    for k in range(len(lst2)):
        lst4=[]
        for i in range(len(lst2[k])):
            if lst2[k][i]!='':
                lst4.append(float(lst2[k][i]))
        lst3.append(lst4)
    file.close()
    
    
    #finding the location of the max stress 

    for i in range(len(lst3)):
        if lst3[i][0]==node:
            if lst3[i][2]>0 and lst3[i][3]>0:
                member=1
            elif lst3[i][2]<0 and lst3[i][3]>0:
                member=3 
            elif lst3[i][2]>0 and lst3[i][3]<0:
                member=2  
            elif lst3[i][2]<0 and lst3[i][3]<0:
                member=4
                
    
    return member


def creating_format_x(number_of_nodes, span):
    delta=span/number_of_nodes
    x_set_points=[]
    x=0
    for i in range(number_of_nodes+1):
        x_set_points.append(x)
        x+=delta
        
    return x_set_points



def converting_deflection_list(new_x_set, deflection_values):
    i=1
    new_list_deflections=[]
    while i<=len(new_x_set)-1:
        sum_y=0
        sum_z=0
        counter=0
        for j in range(len(deflection_values)):
            if new_x_set[i-1]<deflection_values[j][0] and new_x_set[i]>deflection_values[j][0]:
                sum_y+=deflection_values[j][1]
                sum_z+=deflection_values[j][2]
                counter+=1
        if counter==0:
            new_list_deflections.append([new_x_set[i], 0, 0])
        else:
            new_list_deflections.append([new_x_set[i], sum_y/counter, sum_z/counter])
        i+=1
        
    return new_list_deflections

def converting_twist_list(new_x_set, twist_values, x_values):
    i=1
    new_list_twists=[]
    
    while i<=len(new_x_set)-1:
        sum_twist=0
        counter=0
        for j in range(len(x_values)):
            if new_x_set[i-1]<x_values[j] and new_x_set[i]>=x_values[j]:
                sum_twist+=twist_values[j]
                counter+=1
        if counter==0:
            new_list_twists.append([new_x_set[i], 0])
        else:
            new_list_twists.append([new_x_set[i], sum_twist/counter])
        i+=1
        
    return new_list_twists




def testing_space_validation(filename):
    a=reading_inputs(filename)[0]

    a=sorted(a, key = lambda x: x[1])
    
    dif=[]
    i=1
    while i<=len(a)-1:
        dif.append(abs(a[i][1]-a[i-1][1]))
        i+=1

    return dif











