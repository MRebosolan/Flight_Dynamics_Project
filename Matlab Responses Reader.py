import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

#CHANGE VARIABLE NAMES!!!!!!!

#To retrieve a certain series of measurements:
# reference_dataset()[i][j], where i is the index of the parameter and j is: [0] for parameter name, [1] for units,
#[2] for array of measurements.


mat = sio.loadmat("matlab.mat")

def reference_dataset():

    data = mat["flightdata"]
    parameters = []

    temp0 = data[0]
    temp1 = temp0[0]

    for i in range(-1, 47):

        temp2 = temp1[i]
        temp3 = temp2[0]
        temp4 = temp3[0]
        value0 = temp4[0]

        if i == -1:
            value1 = value0[0]
            value = value0.transpose()

        if i != -1:
            value = value0

        if i < 43:
            unit_temp0 = temp4[1]
            unit_temp1 = unit_temp0[0]
            unit_temp2 = unit_temp1[0]
            unit = unit_temp2[0]

        if i >= 43:
            unit = "no units"

        name_temp0 = temp4[2]
        name_temp1 = name_temp0[0]
        name_temp2 = name_temp1[0]

        name = name_temp2[0]


        parameter = [name, unit, value]

        parameters.append(parameter)

        ## Mass f(t)

    mass = 6735.644 #[kg]
    lbs_h_to_kg_s = 0.45359237/3600
    dt = 0.1 #[s]

    mass_val = np.zeros([len(value), 1])
    FF1 = parameters[4][2] * lbs_h_to_kg_s
    FF2 = parameters[5][2] * lbs_h_to_kg_s

    for i in np.arange(len(value)):
        mass = mass - (dt*(FF1[i] + FF2[i]))
        mass_val[i] = mass
        mass_par = ["Aircraft Mass", "kg"], mass_val

    parameters.append(mass_par)

    return parameters

for i in range(0,49):
   print (i, ":", reference_dataset()[i][0], "[", reference_dataset()[i][1], "]")
t1 = 53*60 + 57
t2 = t1 + 360
t = reference_dataset()[0][2]
V = reference_dataset()[41][2]

t_phugoid = []
indexes = []
for i, time_stamp in enumerate(t):
    if time_stamp>t1 and time_stamp<t2:
        t_phugoid.append(time_stamp)
        indexes.append(i)

index1, index2 = indexes[0], indexes[-1]



plt.plot(t_phugoid, V[index1:(index2+1)])
plt.show()

