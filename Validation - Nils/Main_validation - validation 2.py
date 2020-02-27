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
span=2.661
twists_values=
deflection_values=
x_values=


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
    

    

    







