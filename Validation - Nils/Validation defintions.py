# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:25:00 2020

@author: Nils
"""

import numpy as np
import math

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

'''
def reading_input(file_name):
    file=open(file_name)
    lst=[]
    for x in file:
        lst.append(x)
    file.close()
    
    for i in range(len(lst)):


f = open("B737(1).INP")
lst=[]
for x in file:
  for i in range(len(x)):
      for k in range(1,10):
          if x[i]==k:
              lst.append(k)
'''

np.loadtxt('B737(1).INP', usecols=range(1,8))





