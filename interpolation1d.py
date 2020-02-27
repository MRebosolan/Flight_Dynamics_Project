import DataPointReader
import loadfunction
from functools import reduce
import operator
import numpy as np

import math as m
from matplotlib import pyplot
import sympy
value, xlist, z = loadfunction.coordinates()
class interpolation:

    def __init__(self):
        self.data = interpolation.readerx()


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
    def zdirection(self,l,type,magni):
        if type == 'span':
            valuesz = np.array(value[l::41])
        else:
            valuesz = np.array(magni)
        return valuesz
    def zdircord(self,l,type):
        if type == 'span':
            zcord1 = np.array(z[l::41])
        else:
            zcord1 = np.array(xlist[:41])

        return zcord1


    def basis(self,x,j,m,l,type):
        b = (x - self.zdircord(l,type)[m]) / (self.zdircord(l,type)[j] - self.zdircord(l,type)[m])

        return b
    def interpolate(self,x,l,type,magni):
        if type == 'span':
            y = 0
            for j in range(81):
                t = 1
                for m in range(81):
                    if m != j:
                        t = t * self.basis(x, j, m, l, type)

                y += t * self.zdirection(l, type, magni)[j]
            return y
        else:
            y = 0
            for j in range(41):
                t = 1
                for m in range(41):
                    if m != j:
                        t = t * self.basis(x, j, m, l, type)

                y += t * self.zdirection(l, type, magni)[j]
            return y


        #return np.sum([[np.dot(self.basis(x,j),self.zdirection()[j])]for j in range(1,41)])
'''
i = interpolation()

ill = np.linspace(0,-0.5456,400)
y = []
list = []
#for l in range(41):
    #print(l)
    #list.append(l)
l = 10
for element in ill:
    y.append(i.interpolate(element,l,'span',list))
print(len(y))
pyplot.plot(ill,y)
pyplot.show()
'''




