import numpy as np
from math import *
import matplotlib.pyplot as plt
import xlrd
title = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018'
file = xlrd.open_workbook(title)
sheet = file.sheet_by_name('Sheet1')

def cgcalculator(FU,setting):
    fuelstart = 4050
    fuel = fuelstart - int(FU)
    def cgfuel(fuel):
        fumoments = [298.16, 591.18, 879.08, 1165.42, 1448.40, 1732.53, 2014.80, 2298.84, 2581.92, 2866.3, 3150.18, 3434.52,
               3718.52, 4003.23, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.9, 5994.04, 6278.47, 6562.82,
               6846.96, 7131.00, 7415.33, 7699.6, 7984.34, 8269.06, 8554.05, 8839.04, 9124.8, 9410.62, 9696.97, 9983.4,
               10270.08, 10556.84, 10843.87, 11131.0, 11418.2, 11705.5, 11993.31, 12281.18, 12569.04, 12856.86,
               13144.73, 13432.48, 13720.56, 14008.46, 14320.34]

        mass = np.arange(100, 5100, 100)
        masscrt = (fuel / 100) * 100
        i = np.where(mass == masscrt)[0][0]
        m1 = mass[i]
        m2 = mass[i + 1]
        M1 = fumoments[i]
        M2 = fumoments[i + 1]
        kmom = (M1 + (fuel - m1) * (M2 - M1) / (m2 - m1)) * 100 * i2cm * p2kg
        return kmom
    i2cm = 2.54
    cm2i = 1. / i2cm
    p2kg = 0.453592
    kg2p = 1. / p2kg
    W1 = 95.
    Xcg1 = 131 * i2cm
    W2 = 92.
    Xcg2 = 131 * i2cm
    W3 = 66.
    Xcg3 = 214 * i2cm
    W4 = 61.
    Xcg4 = 214 * i2cm
    W5 = 75.
    Xcg5 = 251 * i2cm
    W6 = 78.
    Xcg6 = 251 * i2cm
    W7 = 86.
    Xcg8 = 288 * i2cm
    W8 = 68.
    if setting != 1:
        Xcg7 = 134 * i2cm
    else:
        Xcg7=288*i2cm
    W10=74.
    Xcg10=170*p2kg
    bew = 9165 * p2kg
    bewcg = 291.65*i2cm
    totcgpayload=W1*Xcg1+W2*Xcg2+W3*Xcg3+W4*Xcg4+W5*Xcg5+W6*Xcg6+W7*Xcg7+W8*Xcg8+W10*Xcg10
    totmasspayload =W1+W2+W3+W4+W5+W6+W7+W8+W10
    cgfuel = cgfuel(fuel)
    CG = (cgfuel + totcgpayload + bewcg) / (totmasspayload + bew + fuel * p2kg)
    Xcg = CG - 261.45 * i2cm
    Weight = fuel * p2kg + totmasspayload + bew
    return CG,Xcg,Weight, (CG*cm2i)




