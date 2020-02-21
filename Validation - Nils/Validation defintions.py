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

def reading_outputs_deflections(file_name, interest_nodes):
    
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

def computing_outputs_twist(shear_center_y, shear_center_z):
    
    deflections_y_and_z_hinge_lines=reading_outputs_deflections('B737(2).RPT', reading_inputs('B737(1).INP'))[1]
    print(deflections_y_and_z_hinge_lines)
    lst_final_values=[]
    for i in range(len(deflections_y_and_z_hinge_lines)):
        twist=np.arctan((deflections_y_and_z_hinge_lines[i][3]-shear_center_y)/(shear_center_z-deflections_y_and_z_hinge_lines[i][4]))
        lst_final_values.append([deflections_y_and_z_hinge_lines[i][0],deflections_y_and_z_hinge_lines[i][1], deflections_y_and_z_hinge_lines[i][2], deflections_y_and_z_hinge_lines[i][3], deflections_y_and_z_hinge_lines[i][4], twist])
    return lst_final_values    



