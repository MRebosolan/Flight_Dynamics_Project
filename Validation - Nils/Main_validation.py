# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:43:19 2020

@author: Nils
"""
import sys, os
sys.path.clear()
sys.path.append(os.path.realpath('..\\svvproject'))
sys.path.append(os.path.realpath('..\\Validation - Nils'))
import Validation_definitions
import matplotlib.pyplot as plt
#VALIDATION DATA
    #reading and computing twists, max_stress and deflections
    #deflections are of the hinge line
    #twists are taken about the hinge line

interest_nodes = Validation_definitions.reading_inputs('B737(1).INP')

outputs_deflections_bent= Validation_definitions.reading_outputs_deflections_bent('B737(2).RPT', interest_nodes[0])[1]
outputs_deflections_unbent = Validation_definitions.reading_outputs_deflections_unbent('B737(2).RPT', interest_nodes[0])[1]
outputs_deflections_pure_bending= Validation_definitions.reading_outputs_deflections_pure_bending('B737(2).RPT', interest_nodes[0])[1]

average_y_z_deflections_due_torque=Validation_definitions.average_y_z_deflections_due_to_torque_validated_data(outputs_deflections_unbent, outputs_deflections_bent, outputs_deflections_pure_bending)

shear_center_l_e_z=-0.10856995078063854
hinge_line_to_l_e_z=0.1025
shear_center_hinge_line_z=shear_center_l_e_z+hinge_line_to_l_e_z

    #OUTPUTS
twists_validation_data= Validation_definitions.computing_outputs_twist(average_y_z_deflections_due_torque, 0,shear_center_hinge_line_z) #outputs a list of sublists [deflection y and z due to torque, twist]
deflections_validation_data=outputs_deflections_bent

max_stress=(423.467*10**(-3)+414.020*10**(-3))/2
max_stress_node=2391
max_stress_location=Validation_definitions.reading_max_stress('B737(1).INP', max_stress_node)
max_stress_member=max_stress_location
print(max_stress_location)

min_stress=(178.435*10**(-6)+155.774*10**(-6))/2
min_stress_node= 6519

#NUMERICAL DATA
span=2.771
twists_values=
deflection_values=
x_values=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.5700000000000001, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.6900000000000001, 0.7000000000000001, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.8200000000000001, 0.8300000000000001, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.9400000000000001, 0.9500000000000001, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.1300000000000001, 1.1400000000000001, 1.1500000000000001, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.3800000000000001, 1.3900000000000001, 1.4000000000000001, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.6300000000000001, 1.6400000000000001, 1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001, 1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27, 2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54, 2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.7600000000000002, 2.77]


max_stress_actual=

max_stress_member_actual=




#Adjusting so that the numerical set of twist and deflections have the same number of points along the span
new_x_set=Validation_definitions.creating_format_x(108, span)
len_twists=len(twists_values)

    #converting the deflections to 108 nodes
new_deflections_num_model= Validation_definitions.converting_deflection_list(new_x_set, deflection_values)
    #converting the twist to 108 nodes
new_twists_num_model= Validation_definitions.converting_twist_list(new_x_set, twists_values, x_values)


#print(new_deflections_num_model)

#ENGAGING IN THE COMPARISON

    #taking the twist values of twists_list
Twists=[]
nodes=[]
for i in range(len(twists_validation_data)):
    Twists.append(twists_validation_data[i][2])
    nodes.append(twists_validation_data[i][0])
    
    #taking the twist values of twists_list actual
    
Twists_actual=[]
for i in range(len(new_twists_num_model)):
    Twists_actual.append(new_twists_num_model[i][1])
    
    #taking the deflection in y and z
new_deflections_num_model_y=[]
new_deflections_num_model_z=[]
for i in range(len(new_deflections_num_model)):
    new_deflections_num_model_y.append(new_deflections_num_model[i][1])
    new_deflections_num_model_z.append(new_deflections_num_model[i][2])
    
    #validation deflections
deflections_validation_data_y=[]
deflections_validation_data_z=[]
for i in range(len(outputs_deflections_bent)):
    deflections_validation_data_y.append(outputs_deflections_bent[i][3])
    deflections_validation_data_z.append(outputs_deflections_bent[i][4])
    


    #plotting
plt.plot(nodes,Twists, color='blue')
plt.plot(nodes,Twists_actual, color= 'red' )
plt.show()

plt.plot(nodes, new_deflections_num_model_y, color='red')
plt.plot(nodes,deflections_validation_data_y, color= 'blue' )
plt.show()

print(max(deflections_validation_data_y))

plt.plot(nodes, new_deflections_num_model_z, color='red')
plt.plot(nodes,deflections_validation_data_z, color= 'blue' )
plt.show()

#-----------------------------------------------------------------------------
lst_error_deflection= Validation_definitions.deflection_validation(new_deflections_num_model, deflections_validation_data)

lst_error_y=[]
for i in range(len(lst_error_deflection)):
    lst_error_y.append(lst_error_deflection[i][1])

lst_error_z=[]
for i in range(len(lst_error_deflection)):
    lst_error_z.append(lst_error_deflection[i][2])
    

plt.ylim(-0.01, 0.01)
plt.plot(nodes, lst_error_y)
plt.show()

plt.ylim(-0.01, 0.01)
plt.plot(nodes, lst_error_z)
plt.show()



lst_variables_twists= Validation_definitions.twist_validation(Twists_actual, Twists)
lst_error_twists=lst_variables_twists[0]


plt.ylim(-0.01, 0.4)
plt.plot(nodes, lst_error_twists)
plt.show()

maximum_minimum_difference_actual=lst_variables_twists[1]
maximum_minimum_difference_validation=lst_variables_twists[2]

lst_max_stress_variables=Validation_definitions.max_stress_validation(max_stress_actual, max_stress, max_stress_member_actual, max_stress_member)
error_max_stress=lst_max_stress_variables[0]

if lst_max_stress_variables[1]==1:
    print('The stress is on the same member')
else: 
    print('The maximum stress is not on the same member')
    

    

    







