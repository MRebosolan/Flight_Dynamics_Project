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
        y_error=abs(deflection_set_actual[i][1]-deflection_set_validation_data[i][1])*100/(deflection_set_validation_data[i][1])
        z_error=abs(deflection_set_actual[i][2]-deflection_set_validation_data[i][2])*100/(deflection_set_validation_data[i][2])
        error_list.append([deflection_set_actual[i][0], y_error, z_error])
    return error_list

def twist_validation(twist_set_actual, twist_set_validation_data):
    error_list=[]
    for i in range(len(twist_set_actual)):
        twist_error=abs(twist_set_actual[i][1]-twist_set_validation_data[i][1])*100/twist_set_validation_data[i][1]
        error_list.append(twist_error)
    
    twists_actual=[]
    twists_validation=[]
    for i in range(len(twist_set_actual)):
        twists_actual.append(twist_set_actual[i][1])
        twists_validation.append(twist_set_validation_data[i][1])
        
    maximum_minimum_difference_actual=max(twists_actual)-min(twists_actual)
    maximum_minimum_difference_validation=max(twists_validation)-min(twists_validation)
    
    return error_list, maximum_minimum_difference_actual, maximum_minimum_difference_validation

def max_stress_validation(max_stress_actual, max_stress_validation, x_location_actual_percent_span, x_location_validation_percent_span, actual_member_max_stress, validation_member_max_stress):
    error_max_stress=abs(max_stress_actual-max_stress_validation)*100/max_stress_validation
    
    if actual_member_max_stress==validation_member_max_stress:
        switch=1
    else:
        switch=0
        
    diff_x=abs(x_location_actual_percent_span-x_location_validation_percent_span)*100/x_location_validation_percent_span
        
    #the switch which is returned labels wether the maximum stress is located on the same member of the aileron cross section or not
    #these are defined at the beginning of the program of defintions for validation
    return error_max_stress, switch, diff_x

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
    
    #lst 4 gives an output list of each line [node number, x, y, z]. The outputs are only the nodes which are on the hinge line
    return lst4

def reading_inputs_shear_center_valida(file_name):
    #this is for Shear Center
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
        if lst3[i][1]==0 and lst3[i][3]>0:
            lst4.append(lst3[i])
            
    lst5=[]
    for i in range(len(lst3)):
        if lst3[i][1]==0 and lst3[i][3]<0:
            lst5.append(lst3[i])
            
    lst6=[]
    for i in range(len(lst3)):
        if lst3[i][1]==0 and lst3[i][3]==0:
            lst6.append(lst3[i])
            
    lst6=sorted(lst6, key = lambda x: x[2]) 
    lst5=sorted(lst5, key = lambda x: x[2]) 

    return lst4, lst5, lst6

def reading_shear_values(file_name, interest_nodes1, interest_nodes2, interest_nodes3):
    file=open(file_name)
    lst=[]
    for x in file:
        a=re.split(' |, |\|n|\n', x)
        lst.append(a)

    lst2=[]
    for i in range(20,5798):
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
        for j in range(len(interest_nodes1)):
            if lst3[i][0]==interest_nodes1[j][0]:
                lst4.append([lst3[i][0], 10**-6*(1.1*lst3[i][4]+1.1*lst3[i][5])/2])
    lst5=[]
    for i in range(len(lst3)):
        for j in range(len(interest_nodes2)):
            if lst3[i][0]==interest_nodes2[j][0]:
                lst5.append([lst3[i][0], 10**-6*(1.1*lst3[i][4]+1.1*lst3[i][5])/2])
                
    lst6=[]
    for i in range(len(lst3)):
        for j in range(len(interest_nodes3)):
            if lst3[i][0]==interest_nodes3[j][0]:
                lst6.append([lst3[i][0], 10**-6*(1.1*lst3[i][4]+1.1*lst3[i][5])/2])
    return lst4, lst5, lst6

def calculating_shear_center_valida(interest_nodes_z_positive, interest_nodes_z_negative, int_nodes_z_neutral, aileron_height_B737, shear_values_B737_z_positive, shear_values_B737_z_negative, B_737_spar_shear, chord_length_B737):

    radius=aileron_height_B737*0.5
    i=1
    moment_sum=0
    while i<=len(interest_nodes_z_positive)-1:
        shear_value=shear_values_B737_z_positive[i-1][1]*abs(np.arctan(interest_nodes_z_positive[i][2]/interest_nodes_z_positive[i][3])-np.arctan(interest_nodes_z_positive[i-1][2]/interest_nodes_z_positive[i-1][3]) )*radius
        moment_sum+=aileron_height_B737*0.5*shear_value
        i+=1
    
    angle=np.arcsin(aileron_height_B737*0.5/(math.sqrt(chord_length_B737**2+aileron_height_B737**2*0.25)))
    
    i=1
    while i<=len(interest_nodes_z_negative)-1:
        shear_value=shear_values_B737_z_negative[i-1][1]*abs(math.sqrt(interest_nodes_z_negative[i][2]**2+interest_nodes_z_negative[i][3]**2)-math.sqrt(interest_nodes_z_negative[i-1][2]**2+interest_nodes_z_negative[i-1][3]**2))
        if interest_nodes_z_negative[i][2]>0:
            moment_sum+=(interest_nodes_z_negative[i][2]*shear_value*np.cos(angle)-interest_nodes_z_negative[i][3]*shear_value*np.sin(angle))
        elif interest_nodes_z_negative[i][2]<0:
            moment_sum+=-1*(interest_nodes_z_negative[i][2]*shear_value*np.cos(angle)+interest_nodes_z_negative[i][3]*shear_value*np.sin(angle))
          
        i+=1
        
    i=1   
    sumy=0
    sums=[]
    while i<=len(interest_nodes_z_positive)-1:
        shear_value=shear_values_B737_z_positive[i-1][1]*abs(np.arctan(interest_nodes_z_positive[i][2]/interest_nodes_z_positive[i][3])-np.arctan(interest_nodes_z_positive[i-1][2]/interest_nodes_z_positive[i-1][3]) )*radius
        sumy+=shear_value*np.sin(np.pi*0.5-abs(np.arctan(interest_nodes_z_positive[i][2]/interest_nodes_z_positive[i][3])))
        sums.append(sumy)
        i+=1
    sums2=[]
    i=1
    while i<=len(interest_nodes_z_negative)-1:
        shear_value=shear_values_B737_z_negative[i-1][1]*abs(math.sqrt(interest_nodes_z_negative[i][2]**2+interest_nodes_z_negative[i][3]**2)-math.sqrt(interest_nodes_z_negative[i-1][2]**2+interest_nodes_z_negative[i-1][3]**2))
        sumy+=-1*shear_value*np.sin(np.arcsin(aileron_height_B737*0.5/(math.sqrt(chord_length_B737**2+aileron_height_B737**2*0.25))))
        sums2.append(sumy)
        i+=1
    
    sum3=[]
    i=1
    while i<=len(int_nodes_z_neutral)-1:
        shear_value=B_737_spar_shear[i-1][1]*(abs(int_nodes_z_neutral[i][2]-int_nodes_z_neutral[i-1][2]))
        sumy+=-1*shear_value
        sum3.append(sumy)
        i+=1
    
    V_y=sumy
    
    z_pos_S_C= moment_sum/V_y
    
    return z_pos_S_C, V_y, moment_sum,sums,sums2,sum3


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
                lst4.append(lst3[i])
                
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
                lst4.append(lst3[i])
                
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
                lst4.append(lst3[i])
                
    return len(lst4), lst4

def average_y_z_deflections_due_to_torque_validated_data(deflection_y_z_unbent, deflection_y_z_bent, deflection_y_z_bending):
    lst=[]
    for i in range(len(deflection_y_z_unbent)):
        D=(deflection_y_z_unbent[i][3]+(deflection_y_z_bent[i][3]-deflection_y_z_bending[i][3]))/2
        C=(deflection_y_z_unbent[i][4]+(deflection_y_z_bent[i][4]-deflection_y_z_bending[i][4]))/2
        E=math.sqrt(C**2+D**2)
        lst.append([deflection_y_z_unbent[i][0],D,C,E])
    
    return lst

def computing_outputs_twist(deflections_y_and_z_hinge_lines, shear_center_y, shear_center_z):
    #these deflection in y and z are calculated by take the average of the unbent case and the dif. between the bent and pure bending case
    
    lst_final_values=[]
    for i in range(len(deflections_y_and_z_hinge_lines)):
        twist=np.arctan((deflections_y_and_z_hinge_lines[i][3]-shear_center_y)/(shear_center_z-deflections_y_and_z_hinge_lines[i][4]))
        lst_final_values.append([deflections_y_and_z_hinge_lines[i][0],deflections_y_and_z_hinge_lines[i][1], deflections_y_and_z_hinge_lines[i][2], deflections_y_and_z_hinge_lines[i][3], deflections_y_and_z_hinge_lines[i][4], twist])
    return lst_final_values    






lst1, lst2, lst3 = reading_inputs_shear_center_valida("B737(1).INP")
print(len(lst3))
lst4, lst5, lst6 =reading_shear_values("B737(2).RPT", lst1, lst2, lst3)
print(len(lst6[0]))
Z=calculating_shear_center_valida(lst1, lst2,lst3, 205, lst4, lst5, lst6, 605)
print(Z)





