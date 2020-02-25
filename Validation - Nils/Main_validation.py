# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:43:19 2020

@author: Nils
"""

import numpy as np
import math
import Validation_defintions

#VALIDATION DATA
    #reading and computing twists, max_stress and deflections
    #deflections are of the hinge line
    #twists are taken about the hinge line

interest_nodes = Validation_defintions.reading_inputs('B737(1).INP')

outputs_deflections_bent= Validation_defintions.reading_outputs_deflections_bent('B737(2).RPT', interest_nodes[0])[1]
outputs_deflections_unbent = Validation_defintions.reading_outputs_deflections_unbent('B737(2).RPT', interest_nodes[0])[1]
outputs_deflections_pure_bending= Validation_defintions.reading_outputs_deflections_pure_bending('B737(2).RPT', interest_nodes[0])[1]

average_y_z_deflections_due_torque=Validation_defintions.average_y_z_deflections_due_to_torque_validated_data(outputs_deflections_unbent, outputs_deflections_bent, outputs_deflections_pure_bending)

shear_center_l_e_z=-0.10856995078063854
hinge_line_to_l_e_z=0.1025
shear_center_hinge_line_z=shear_center_l_e_z+hinge_line_to_l_e_z

twists_validation_data= Validation_defintions.computing_outputs_twist(average_y_z_deflections_due_torque, 0,shear_center_hinge_line_z)
deflections_validation_data=outputs_deflections_bent

max_stress=(423.467*10**(-3)+414.020*10**(-3))/2
max_stress_node=2391
max_stress_location=Validation_defintions.reading_max_stress('B737(1).INP', max_stress_node)
max_stress_member=max_stress_location

min_stress=(178.435*10**(-6)+155.774*10**(-6))/2
min_stress_node= 6519

#NUMERICAL DATA
twists_values=[] #[twist] twist values at each x, x is not specified but is done from 0 to span
deflection_values=[] # [x location, y deflection, z deflection]
x_values=[]

max_stress=
max_stress_node=
max_stress_member=

min_stress=
min_stress_node=


#Adjusting so that the numerical set of twist and deflections have the same number of points along the span
new_x_set=creating_format_x(100, span)
len_twists=len(twists_values)

i=0
while i<=len(new_x_set):
    

#Engaging in the comparison







