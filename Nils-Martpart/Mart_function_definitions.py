# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:24:22 2020

@author: Mart
"""

import math
import numpy as np

def normal_stress_x_bending_function(bending_moment_z, moment_of_inertia_y, y_distance_to_centroid, bending_moment_y, moment_of_inertia_z, z_distance_to_centroid):
    normal_stress_x_bending = ((bending_moment_z*moment_of_inertia_y)*y_distance_to_centroid + bending_moment_y*moment_of_inertia_z*z_distance_to_centroid)/moment_of_inertia_z/moment_of_inertia_y
    return normal_stress_x_bending

def von_mises_stress_function(normal_stress_x_bending, shear_stress_total):
    von_mises_stress = math.sqrt(normal_stress_x_bending**2+3*shear_stress_total**2)
    return von_mises_stress

def redundant_shear_flow(q_b_s1, s1, q_b_s4, s4, aileron_radius, aileron_angle_radians, aileron_height, area_circular_section, area_triangular_section, ratio_between_redundant_shear_flows):
    moment_sum = 0
    for i in range(len(s1)):
        moment_sum = moment_sum+s1[i]*q_b_s1[i]*aileron_radius*2
    for i in range(len(s4)):
        moment_sum = moment_sum+s4[i]*q_b_s4[i]*np.sin(aileron_angle_radians)*aileron_height/2
    redundant_shear_flow_circular_section = -moment_sum/(2*area_circular_section+2*area_triangular_section*ratio_between_redundant_shear_flows)
    redundant_shear_flow_triangular_section = redundant_shear_flow_circular_section*ratio_between_redundant_shear_flows
    return redundant_shear_flow_circular_section, redundant_shear_flow_circular_section

def shear_centre(q_total_top1, q_total_bottom1, q_total_top2, q_total_bottom2, aileron_height, aileron_angle_radians):
    eta = (q_total_top1+q_total_bottom1)*aileron_height/2+(q_total_top2+q_total_bottom2)*np.sin(aileron_angle_radians)*aileron_height/2
    return eta
