import numpy as np
def readerx():
    fin = open("aerodynamicloada320.dat")
    data_points = []

    for line in fin:
        row = []
        line = line.split(",")
        for value in line:
            row.append(float(value))
        data_points.append(np.array(row))

    data_points = np.array(data_points)

    return data_points
