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
from Mart_function_definitions import normal_stress_x_bending_function
from definitions_max_stress_deflections import relation_shear_1_and_2_torque_and_torque

"--------------------------------------------------------------------------------------------------------------------------"
"This section gives all the dummy inputs to check whether our model can handle the inputs"

number_of_sections = 100
moment_about_y_x_direction = number_of_sections*[100]
moment_about_z_x_direction = number_of_sections*[300]
shear_force_y_x_direction = number_of_sections*[30]
shear_force_z_x_direction = number_of_sections*[10]
torque_x_direction = number_of_sections*[15]
moment_of_inertia_y = 1000
moment_of_inertia_z = 500
centroid_y = 0.0
centroid_z = 0.4
aileron_height = 1
skin_thickness = 0.001
spar_thickness = 0.005
chord_length = 5


print("reached part 1")
"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the shear flow due to torsion"

shear_flow_torque_circular_part_list = []
shear_flow_torque_triangular_part_list = []

for torque in torque_x_direction:
    shear_flow_torque_circular_part, shear_flow_torque_triangular_part, ratio_between_redundant_shear_flows = relation_shear_1_and_2_torque_and_torque(torque, aileron_height, skin_thickness, spar_thickness, chord_length)
    shear_flow_torque_circular_part_list.append(shear_flow_torque_circular_part)
    shear_flow_torque_triangular_part_list.append(shear_flow_torque_triangular_part)
    
print(shear_flow_torque_circular_part_list)
print(shear_flow_torque_triangular_part_list)
print(ratio_between_redundant_shear_flows)

"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the shear flow due to shear, and in that process it will section the cross-section as well"



"--------------------------------------------------------------------------------------------------------------------------"
"This part will handle the computation of shear centre and following that the deflection due to torque"



"--------------------------------------------------------------------------------------------------------------------------"
"This part will cover the bending equations and deflection due to bending"


    

