from math import *
import numpy as np


#inputs = interpolated load function w(x,z)

# Integral to find resultant force is \int_{xprime1}^{xprime2} \int_{zprime1}^{zprime2} f(xprime,zprime) dzprime dxprime
#2001 rows, 1001 columns

#Integration is performed row-wise first and then column wise using trapeziodal rule

#loading_function is a dummy

def numerical_integrator(loading_function, nr_of_rows, nr_of_columns):

    zprime1 = 0
    zprime2 = -2.7608485973438723
    xprime1 = 0
    xprime2 = -0.5464859422096491
    total = 0

    zprimes = np.linspace(zprime1, zprime2, nr_of_rows) #array of z coordinates
    xprimes = np.linspace(xprime1, xprime2, nr_of_columns) #array of x coordinates
    dummy_loading_function = np.random.rand(nr_of_rows, nr_of_columns)
    row_sums = []                                           #contains all row integrals
    total_sum = 0

    for row in dummy_loading_function:
        row_sum = 0
        for c, element in enumerate(row[:-1]):         #i is index of element in the row
            row_sum += (row[c] + row[c+1]) * (xprimes[c+1] - xprimes[c]) * 0.5
        row_sums.append(row_sum)

    for r, row in enumerate(row_sums[:-1]):
        total_sum += (row_sums[r] + row_sums[r+1]) * (zprimes[r+1] - zprimes[r]) * 0.5

    return total_sum

print(numerical_integrator(1,2001,1001))
