# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:25:00 2020

@author: Nils
"""

import numpy as np
import math


def deflection_validation(deflection_set_actual, deflection_set_validation_data):
    error_list=[]
    for i in range(len(deflection_set_actual)):
        y_error=abs(deflection_set_actual[i][1]-deflection_set_validation_data[i][1])*100/(deflection_set_validation_data[i][1])
        z_error=abs(deflection_set_actual[i][2]-deflection_set_validation_data[i][2])*100/(deflection_set_validation_data[i][2])
        error_list.append([deflection_set_actual[i][0], y_error, z_error])
    return error_list