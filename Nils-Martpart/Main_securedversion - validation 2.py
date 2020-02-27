# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:42:44 2020

@author: Mart
"""

"Main programme that takes all the functions and makes the main file."

"As inputs for this programme are assumed: moment, torque and shear distributions as lists over x-direction, and section properties Moment of Inertia and centroid"

"As outputs for this programme deflections and max stress is given"

"--------------------------------------------------------------------------------------------------------------------------"
"This section imports stuff, mainly the definitions of functions can then be used"

import math 
import numpy as np
import sys, os
sys.path.clear()
sys.path.append(os.path.realpath('..\\svvproject'))
sys.path.append(os.path.realpath('..\\Nils-Martpart'))
from Mart_function_definitions import normal_stress_x_bending_function, redundant_shear_flow, shear_centre, normal_stress_x_bending_function, von_mises_stress_function
from definitions_max_stress_deflections import relation_shear_1_and_2_torque_and_torque, alternative_q_base_top1, alternative_q_base_spar, alternative_q_base_bottom1, alternative_q_base_sparA, alternative_q_base_top2, alternative_q_base_bottom2, alternative_q_base_sparB, deflection_z_bending_stress, deflection_y_bending_stress, rate_twist_at_x, twist, deflection_due_to_torque_and_bending

"--------------------------------------------------------------------------------------------------------------------------"
"This section gives all the dummy inputs to check whether our model can handle the inputs"

number_of_sections = 266
aileron_length = 2.661
moment_about_y_x_direction =
moment_about_z_x_direction = 
shear_force_y_x_direction = 
shear_force_z_x_direction = 
torque_x_direction = 
spanwise_locations = np.arange(0, aileron_length , round(aileron_length/number_of_sections,4))
moment_of_inertia_y = 8.710953452731196*10**(-5)
moment_of_inertia_z = 1.0510952632762589*10**(-5)
centroid_y = 0.0
centroid_z = 0.1297988808236754 
aileron_height = 0.205
skin_thickness = 0.0011
spar_thickness = 0.0028
chord_length = 0.605
integral_step = 0.01
aileron_radius = aileron_height/2
aileron_angle_radians = 1/5*math.pi
area_circular_section = 0.5*math.pi*aileron_radius**2
area_triangular_section = 0.5*aileron_height*(chord_length-aileron_height/2)
x_location_hinge1 = 0.172
x_location_hinge3 = 2.591
deflection_hinge_1 = 0.01154 *np.cos(28*np.pi/180)
deflection_hinge_3 = 0.0184 *np.cos(28*np.pi/180)
boom_list=[[4.2e-05, 0.0, -0.1025], [4.2e-05, 0.07160610879327349, -0.07334040621298657], [4.2e-05, -0.07160610879327349, -0.07334040621298657], [4.2e-05, 0.10247065573131683, -0.002452491385090335], [4.2e-05, -0.10247065573131683, -0.002452491385090335], [4.2e-05, 0.08714549999999999, 0.0752745], [4.2e-05, -0.08714549999999999, 0.0752745], [4.2e-05, 0.07130105, 0.15295094999999997], [4.2e-05, -0.07130105, 0.15295094999999997], [4.2e-05, 0.055456599999999995, 0.23062739999999998], [4.2e-05, -0.055456599999999995, 0.23062739999999998], [4.2e-05, 0.023766675, 0.38598532499999993], [4.2e-05, -0.023766675, 0.38598532499999993], [4.2e-05, 0.007922224999999991, 0.463661775], [4.2e-05, -0.007922224999999991, 0.463661775]]
elasticity_modulus=73.1*10**9
shear_modulus=28*10**9
"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the shear flow due to torsion. It is the same in every cross-section, so it gives outputs only as fixed values per thingey"

shear_flow_torque_circular_part_list = []
shear_flow_torque_triangular_part_list = []

for torque in torque_x_direction:
    shear_flow_torque_circular_part, shear_flow_torque_triangular_part, ratio_between_redundant_shear_flows = relation_shear_1_and_2_torque_and_torque(torque, aileron_height, skin_thickness, spar_thickness, chord_length)
    shear_flow_torque_circular_part_list.append(shear_flow_torque_circular_part)
    shear_flow_torque_triangular_part_list.append(shear_flow_torque_triangular_part)
#print(shear_flow_torque_circular_part_list)
#print(shear_flow_torque_triangular_part_list)
#print(ratio_between_redundant_shear_flows)

"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the shear flow due to shear, and in that process it will section the cross-section as well"

MOI_y_prime = moment_of_inertia_y
MOI_z_prime = moment_of_inertia_z

maximum_stress_s1 = []
location_maximum_stress_s1 = []
maximum_stress_s2 = []
location_maximum_stress_s2 = []
maximum_stress_s3 = []
location_maximum_stress_s3 = []
maximum_stress_s4 = []
location_maximum_stress_s4 = []
maximum_stress_s5 = []
location_maximum_stress_s5 = []

rate_of_twist_x_list = []
switch=0
for j in range(0, len(shear_force_y_x_direction)):
    shear_force_y = shear_force_y_x_direction[j]
    shear_force_z = shear_force_z_x_direction[j]
    
    z_pos = 0
    y_pos = aileron_height/2
    q_base_list_top_1, s1_list  = alternative_q_base_top1(shear_force_y, shear_force_z, moment_of_inertia_y, moment_of_inertia_z, integral_step, aileron_height, skin_thickness, z_pos, y_pos, boom_list)

    z_pos = 0
    y_pos = -aileron_height/2
    q_base_list_spar_area_1, s5_list = alternative_q_base_spar(shear_force_y, shear_force_z, moment_of_inertia_y, moment_of_inertia_z, chord_length, integral_step, aileron_height, skin_thickness, y_pos)
    
    z_pos = -aileron_height/2
    y_pos = 0
    q_base_list_bottom_1, s3_list = alternative_q_base_bottom1(shear_force_y, shear_force_z, moment_of_inertia_y, moment_of_inertia_z, integral_step, aileron_height, skin_thickness, z_pos, y_pos, boom_list)
    
    z_pos = 0
    y_pos = aileron_height/2
    q_base_list_spar_area_2_A, sA_list = alternative_q_base_sparA(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness, y_pos)

    z_pos = chord_length-aileron_height/2
    y_pos = 0
    q_base_list_top_2, s2_list = alternative_q_base_top2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness, z_pos, y_pos, boom_list)
    
    z_pos = 0
    y_pos = -aileron_height/2
    q_base_list_bottom_2, s4_list = alternative_q_base_bottom2(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness, z_pos, y_pos, boom_list)
    
    z_pos = 0
    y_pos = 0
    q_base_list_spar_area_2_B, sB_list = alternative_q_base_sparB(shear_force_y, shear_force_z, MOI_y_prime, MOI_z_prime, chord_length, integral_step, aileron_height, skin_thickness, y_pos)

    "The upcoming part adds the base flows of the previous parts to the current part"
    
    
    for i in range(0, len(q_base_list_spar_area_1)):
        q_base_list_spar_area_1[i] = q_base_list_spar_area_1[i]+q_base_list_top_1[-1]
    
    for i in range(0, len(q_base_list_bottom_1)):
        q_base_list_bottom_1[i] = q_base_list_bottom_1[i]+q_base_list_spar_area_1[-1]
        
    for i in range(0, len(q_base_list_top_2)):
        q_base_list_top_2[i] = q_base_list_top_2[i]+q_base_list_spar_area_2_A[-1]
        
    for i in range(0, len(q_base_list_bottom_2)):
        q_base_list_bottom_2[i] = q_base_list_bottom_2[i]+q_base_list_top_2[-1]
        
    for i in range(0, len(q_base_list_spar_area_2_B)):
        q_base_list_spar_area_2_B[i] = q_base_list_spar_area_2_B[i]+q_base_list_bottom_2[-1]

    "The upcoming part adds the shear flows of the spars together to get to one value"
    
    q_base_list_spar_A = []
    for i in range(0, len(q_base_list_spar_area_2_A)):
        index=int(-(len(q_base_list_spar_area_1)/2+i))
        q_base_list_spar_A.append(q_base_list_spar_area_2_A[i]-q_base_list_spar_area_1[index])
        
    q_base_list_spar_B = []
    for i in range(0, len(q_base_list_spar_area_2_B)):
        q_base_list_spar_B.append(q_base_list_spar_area_2_B[i]+q_base_list_spar_area_1[-i])
        
        "The upcoming part calculates the redundant shear flows for the sections"
    s1_list=s1_list[:-1]
    s2_list=s2_list[:-1]
    s3_list=s3_list[:-1]
    s4_list=s4_list[:-1]
    sA_list=sA_list[:-1]
    sB_list=sB_list[:-1]
    redundant_shear_flow_circular_section, redundant_shear_flow_triangular_section = redundant_shear_flow(q_base_list_top_1, s1_list, q_base_list_bottom_2, s4_list, aileron_radius, aileron_angle_radians, aileron_height, area_circular_section, area_triangular_section, ratio_between_redundant_shear_flows)
    
    "The upcoming part adds the redundant, base and torque shear flows"    
    q_total_list_top_1 = []    
    for i in range(0, len(q_base_list_top_1)):
        q_total_list_top_1.append(q_base_list_top_1[i]+redundant_shear_flow_circular_section+shear_flow_torque_circular_part_list[j])
        
    q_total_list_top_2 = []
    for i in range(0, len(q_base_list_top_2)):
        q_total_list_top_2.append(q_base_list_top_2[i]+redundant_shear_flow_triangular_section+shear_flow_torque_triangular_part_list[j])
        
    q_total_list_bottom_1 = []
    for i in range(0, len(q_base_list_bottom_1)):
        q_total_list_bottom_1.append(q_base_list_bottom_1[i]+redundant_shear_flow_circular_section+shear_flow_torque_circular_part_list[j])
        
    q_total_list_bottom_2 = []
    for i in range(0, len(q_base_list_bottom_2)):
        q_total_list_bottom_2.append(q_base_list_bottom_2[i]+redundant_shear_flow_triangular_section+shear_flow_torque_triangular_part_list[j])
        
    q_total_list_spar_A = []
    for i in range(0, len(q_base_list_spar_A)):
        q_total_list_spar_A.append(q_base_list_spar_A[i]-redundant_shear_flow_circular_section+redundant_shear_flow_triangular_section-shear_flow_torque_circular_part_list[j]+shear_flow_torque_triangular_part_list[j])
        
    q_total_list_spar_B = []
    for i in range(0, len(q_base_list_spar_B)):
        q_total_list_spar_B.append(q_base_list_spar_B[i]-redundant_shear_flow_circular_section+redundant_shear_flow_triangular_section-shear_flow_torque_circular_part_list[j]+shear_flow_torque_triangular_part_list[j])

    q_total_list_spar = q_total_list_spar_B+q_total_list_spar_A
    s5_list = [sB_list + sA_list]


    "This part will define geometry for all of the 5 parts, creating a z-list and y-list next to the already existing s-lists."
    "s1 is the top circular part"
    "s2 is the top triangular part"
    "s3 is the bottom circular part"
    "s4 is the bottom triangular part"
    "s5 is the spar, bottom to top"
    
    s1_y_list = []
    s1_z_list = []
    for s1 in s1_list:
        s1_y_list.append(np.sin(s1))
        s1_z_list.append(-np.cos(s1))
        
    s2_y_list = []
    s2_z_list = []
    for s2 in s2_list:
        s2_y_list.append(aileron_height/2-s2*np.sin(aileron_angle_radians))
        s2_z_list.append(s2*np.cos(aileron_angle_radians))
        
    s3_y_list = []
    s3_z_list = []
    for s3 in s3_list:
        s3_y_list.append(-s3*np.sin(aileron_angle_radians))
        s3_z_list.append((chord_length-0.5*aileron_height)-s3*np.cos(aileron_angle_radians))
        
    s4_y_list = []
    s4_z_list = []
    for s4 in s4_list:
        s4_y_list.append(-aileron_height/2+np.cos(s4))
        s4_z_list.append(-np.sin(s4))
        
    s5_y_list = []
    s5_z_list = []
    for s5 in s5_list[0]:
        append=-1*aileron_height/2+s5
        s5_y_list.append(append)
        s5_z_list.append(0)
        
    "This part will calculate, for every section, the stress throughout. It will return the maximum stress."
    maximum_stress_stored_s1 = 0
    for i in range(len(s1_list)):
        bending_stress = normal_stress_x_bending_function(moment_about_z_x_direction[j], moment_of_inertia_y, s1_y_list[i]-centroid_y, moment_about_y_x_direction[j], moment_of_inertia_z, s1_z_list[i]-centroid_z)
        shear_stress = q_total_list_top_1[i]*skin_thickness
        von_mises_stress = von_mises_stress_function(bending_stress, shear_stress)
        if abs(von_mises_stress) > abs(maximum_stress_stored_s1):
            maximum_stress_stored_s1 = abs(von_mises_stress)
            location_maximum_stress_s1_value = s1_list[i]
            
    maximum_stress_stored_s2 = 0
    for i in range(len(s2_list)):
        bending_stress = normal_stress_x_bending_function(moment_about_z_x_direction[j], moment_of_inertia_y, s2_y_list[i]-centroid_y, moment_about_y_x_direction[j], moment_of_inertia_z, s2_z_list[i]-centroid_z)
        shear_stress = q_total_list_top_2[i]*skin_thickness
        von_mises_stress = von_mises_stress_function(bending_stress, shear_stress)
        if abs(von_mises_stress) > abs(maximum_stress_stored_s2):
            maximum_stress_stored_s2 = abs(von_mises_stress)
            location_maximum_stress_s2_value = s2_list[i]
            
    maximum_stress_stored_s3 = 0
    for i in range(len(s3_list)):
        bending_stress = normal_stress_x_bending_function(moment_about_z_x_direction[j], moment_of_inertia_y, s3_y_list[i]-centroid_y, moment_about_y_x_direction[j], moment_of_inertia_z, s3_z_list[i]-centroid_z)
        shear_stress = q_total_list_bottom_1[i]*skin_thickness
        von_mises_stress = von_mises_stress_function(bending_stress, shear_stress)
        if abs(von_mises_stress) > abs(maximum_stress_stored_s3):
            maximum_stress_stored_s3 = abs(von_mises_stress)
            location_maximum_stress_s3_value = s3_list[i] 
            
    maximum_stress_stored_s4 = 0
    for i in range(len(s4_list)):
        bending_stress = normal_stress_x_bending_function(moment_about_z_x_direction[j], moment_of_inertia_y, s4_y_list[i]-centroid_y, moment_about_y_x_direction[j], moment_of_inertia_z, s4_z_list[i]-centroid_z)
        shear_stress = q_total_list_bottom_2[i]*skin_thickness
        von_mises_stress = von_mises_stress_function(bending_stress, shear_stress)
        if abs(von_mises_stress) > abs(maximum_stress_stored_s4):
            maximum_stress_stored_s4 = abs(von_mises_stress)
            location_maximum_stress_s4_value = s4_list[i]
            
    maximum_stress_stored_s5 = 0
    for i in range(len(s5_list[0])):
        bending_stress = normal_stress_x_bending_function(moment_about_z_x_direction[j], moment_of_inertia_y, s5_y_list[i]-centroid_y, moment_about_y_x_direction[j], moment_of_inertia_z, s5_z_list[i]-centroid_z)
        shear_stress = q_total_list_spar[i]*spar_thickness
        von_mises_stress = von_mises_stress_function(bending_stress, shear_stress)
        if abs(von_mises_stress) > abs(maximum_stress_stored_s5):
            maximum_stress_stored_s5 = abs(von_mises_stress)
            location_maximum_stress_s5_value = s5_list[0][i]     
            
    "This last part will store everything"
    maximum_stress_s1.append(maximum_stress_stored_s1)
    location_maximum_stress_s1.append(location_maximum_stress_s1_value)
    maximum_stress_s2.append(maximum_stress_stored_s2)
    location_maximum_stress_s2.append(location_maximum_stress_s2_value)
    maximum_stress_s3.append(maximum_stress_stored_s3)
    location_maximum_stress_s3.append(location_maximum_stress_s3_value)
    maximum_stress_s4.append(maximum_stress_stored_s4)
    location_maximum_stress_s4.append(location_maximum_stress_s4_value)
    maximum_stress_s5.append(maximum_stress_stored_s5)
    location_maximum_stress_s5.append(location_maximum_stress_s5_value) 
    
    
    rate_of_twist_x = rate_twist_at_x(shear_flow_torque_circular_part_list[j], shear_flow_torque_triangular_part_list[j], redundant_shear_flow_circular_section, redundant_shear_flow_triangular_section, aileron_height, spar_thickness, skin_thickness, shear_modulus)
    rate_of_twist_x_list.append(rate_of_twist_x)
    



list_max_stress=[max(maximum_stress_s1), max(maximum_stress_s2), max(maximum_stress_s3), max(maximum_stress_s4), max(maximum_stress_s5)]
maximum_stress=max(list_max_stress)
index_max_stress=list_max_stress.index(maximum_stress)

member=index_max_stress+1
print(member, maximum_stress)
"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the computation of shear centre"
q_total_top1 = q_total_list_top_1
q_total_bottom1 = q_total_list_bottom_1
q_total_top2 = q_total_list_top_2
q_total_bottom2 = q_total_list_bottom_2

shear_centre_location_wrt_spar = -shear_centre(q_total_top1, q_total_bottom1, q_total_top2, q_total_bottom2, aileron_height, aileron_angle_radians,s1_list, s2_list, s3_list, s4_list )

print(shear_centre_location_wrt_spar)

shear_center_z=-0.10856995078063854+aileron_height/2

"--------------------------------------------------------------------------------------------------------------------------"
"This part will cover the bending equations and deflection due to bending and torque"

integral_values2_z = deflection_z_bending_stress(moment_about_y_x_direction, spanwise_locations, moment_of_inertia_y, elasticity_modulus)
integral_values2_y = deflection_y_bending_stress(moment_about_z_x_direction, spanwise_locations, moment_of_inertia_z, elasticity_modulus)   



twist_total, twist_list = twist(spanwise_locations, rate_of_twist_x_list)[0], twist(spanwise_locations, rate_of_twist_x_list)[1]
lst_Deflections_y_and_z_respect_x = deflection_due_to_torque_and_bending(twist_list, spanwise_locations, 0, shear_center_z, integral_values2_y, integral_values2_z, x_location_hinge1, x_location_hinge3, deflection_hinge_1, deflection_hinge_3)


#print(maximum_stress, member)
#print(lst_Deflections_y_and_z_respect_x)

#print(twist_list)
#print(list(spanwise_locations))

