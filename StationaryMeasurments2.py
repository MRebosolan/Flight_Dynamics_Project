#!/usr/bin/env python3
# -*- coding: utf-8 -*-

<<<<<<< HEAD
@author: Richelle
"""

import math as m 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

print(pd.__version__)
=======
"Created on Fri Mar  6 15:58:53 2020"
import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
>>>>>>> b08434233443e372d3f7b0fbce5e9a67ad96d779
"Import constants"
S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
chord  = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / chord	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * m.pi / 180   # stabiliser angle of incidence [rad]
p0 = 101325 #pa
rho0 = 1.225 #kg/m3
T0 = 288.15 #K
lambda_1 = 1.4 #ratio of heats
lambda0 = -0.0065 #K/m
R = 287 #J/kg*K
mu0 = 1.7894*10**(-5) #kg/m*s
standard_weight = 60500 #Standard weight [N]
engine_inlet_diameter = 0.686 #[m2]

"Calculate total weight at each measure point"
OEW = 9165 * 0.45359237 #lbs to kg
weight_fuel = 4050 * 0.45359237 #lbs to kg
weight_payload = 695 #kg
g0 = 9.80665
"inputs"
# angle of attacks

<<<<<<< HEAD
#title= 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
#file = pd.read_excel(title)
#file1 = file.to_numpy()    
#def valueimporter(row,col1,col2):
#    return file1[row][col1:(col2+1)]
#def valueimport2(row1,row2,col):
#    return file1[row1:(row2+1):col]
#print(valueimport2(26,31,5))
"Transform angle of attack to radians"
def AOAtoRad(AOA):
    AOArad=[]
    for i in AOA:
        i=i*m.pi/180
        AOArad.append(i)
    return(AOArad)

alpha_1 = 1.7 * m.pi / 180 #degree to radians
alpha_2 = 2.4 * m.pi / 180 #degree to radians
alpha_3 = 3.6 * m.pi / 180 #degree to radians
alpha_4 = 5.4 * m.pi / 180 #degree to radians
alpha_5 = 8.7 * m.pi / 180 #degree to radians
alpha_6 = 10.6 * m.pi / 180 #degree to radians
angle_of_attack = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6]
=======
title="Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx"
file = pd.read_excel(title)
file1 = file.to_numpy()
def valueimport2(row1,row2,col,file2):
    for i in range(row1,row2+1):
        file2.append(float(file1[i][col]))
    return file2
IAS=[]
IAS=valueimport2(26,31,4,IAS)
AOA=[]
AOA=valueimport2(26,31,5,AOA)
"Transform angle of attack to radians"
>>>>>>> b08434233443e372d3f7b0fbce5e9a67ad96d779

def degtorad(angledeg):
    aoa=[]
    for i in angledeg:
        anglerad=i*m.pi/180
        aoa.append(anglerad)
    return(aoa)

angle_of_attack0=np.array(degtorad(AOA))

alpha_1 = 1.7 * m.pi / 180 #degree to radians
alpha_2 = 2.4 * m.pi / 180 #degree to radians
alpha_3 = 3.6 * m.pi / 180 #degree to radians
alpha_4 = 5.4 * m.pi / 180 #degree to radians
alpha_5 = 8.7 * m.pi / 180 #degree to radians
alpha_6 = 10.6 * m.pi / 180 #degree to radians
angle_of_attack = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6]


def lbstokg(lbs):
    kg=lbs*0.45359237
    return(kg)
fuel_used1 = 360 * 0.45359237 #lbs to kg
fuel_used2 = 412 * 0.45359237 #lbs to kg
fuel_used3 = 447 * 0.45359237 #lbs to kg
fuel_used4 = 478 * 0.45359237 #lbs to kg
fuel_used5 = 532 * 0.45359237 #lbs to kg
fuel_used6 = 570 * 0.45359237 #lbs to kg

"Calculate weight"
def weight_total(OEW, weight_fuel, weight_payload, fuel_used, g0):
    weight_total = (OEW + weight_fuel + weight_payload - fuel_used) * g0
    return(weight_total)
weight_total1 = (OEW + weight_fuel + weight_payload - fuel_used1) * g0
weight_total2 = (OEW + weight_fuel + weight_payload - fuel_used2) * g0
weight_total3 = (OEW + weight_fuel + weight_payload - fuel_used3) * g0
weight_total4 = (OEW + weight_fuel + weight_payload - fuel_used4) * g0
weight_total5 = (OEW + weight_fuel + weight_payload - fuel_used5) * g0
weight_total6 = (OEW + weight_fuel + weight_payload - fuel_used6) * g0
print(weight_total1)
"Calculate lift"
def lift(weight):
    lift=weight
    return lift
lift1 = weight_total1
lift2 = weight_total2
lift3 = weight_total3
lift4 = weight_total4
lift5 = weight_total5
lift6 = weight_total6

"Transform IAS to CAS"
def IAStoCAS(IAS):
    V_c=(IAS-2)*0.514444444
    return(V_c)
V_c1 = (249 - 2) * 0.514444444 #kts to m/s
V_c2 = (221 - 2) * 0.514444444 #kts to m/s
V_c3 = (192 - 2) * 0.514444444 #kts to m/s
V_c4 = (163 - 2) * 0.514444444 #kts to m/s
V_c5 = (130 - 2) * 0.514444444 #kts to m/s
V_c6 = (118 - 2) * 0.514444444 #kts to m/s

"Transform height to m"
def fttom(hpft):
    hpm=hpft*0.3048
    return(hpm)
hp1 = 5010 * 0.3048 #ft to m
hp2 = 5020 * 0.3048 #ft to m
hp3 = 5020 * 0.3048 #ft to m
hp4 = 5030 * 0.3048 #ft to m
hp5 = 5020 * 0.3048 #ft to m
hp6 = 5110 * 0.3048 #ft to m

"Calculate pressure at specific heights"
def pressure(p0,lambda0,hp,T0,g0,R):
    p=p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    return p
p1 = p0 * (1 + (lambda0*hp1)/T0)**(-g0/(lambda0*R))
p2 = p0 * (1 + (lambda0*hp2)/T0)**(-g0/(lambda0*R))
p3 = p0 * (1 + (lambda0*hp3)/T0)**(-g0/(lambda0*R))
p4 = p0 * (1 + (lambda0*hp4)/T0)**(-g0/(lambda0*R))
p5 = p0 * (1 + (lambda0*hp5)/T0)**(-g0/(lambda0*R))
p6 = p0 * (1 + (lambda0*hp6)/T0)**(-g0/(lambda0*R))

"Calculate Mach number at each measure point"
def mach(p0,p,rho0,V,lambda_1):
    M=m.sqrt((2/(1-lambda_1))*((1 + p0/p1 * ((1 + (1-lambda_1)/(2*lambda_1)*rho0/p0*V_c1**2)**(lambda_1/(1-lambda_1))-1))**((1-lambda_1)/lambda_1) -1))
    return M
M1 = m.sqrt((2/0.4)*((1 + p0/p1 * ((1 + 0.4/2.8 *rho0/p0*V_c1**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
print(M1)
M2 = m.sqrt((2/0.4)*((1 + p0/p2 * ((1 + 0.4/2.8 *rho0/p0*V_c2**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M3 = m.sqrt((2/0.4)*((1 + p0/p3 * ((1 + 0.4/2.8 *rho0/p0*V_c3**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M4 = m.sqrt((2/0.4)*((1 + p0/p4 * ((1 + 0.4/2.8 *rho0/p0*V_c4**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M5 = m.sqrt((2/0.4)*((1 + p0/p5 * ((1 + 0.4/2.8 *rho0/p0*V_c5**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M6 = m.sqrt((2/0.4)*((1 + p0/p6 * ((1 + 0.4/2.8 *rho0/p0*V_c6**2)**(1.4/0.4)-1))**(0.4/1.4) -1))

"Measured total temperature"
def TemptoK(T):
    TAT=T+273.15
    return TAT
TAT1 = 12.5 + 273.15 #Celsius to Kelvin
TAT2 = 10.5 + 273.15 #Celsius to Kelvin
TAT3 = 8.8 + 273.15 #Celsius to Kelvin
TAT4 = 7.2 + 273.15 #Celsius to Kelvin
TAT5 = 6 + 273.15 #Celsius to Kelvin
TAT6 = 5.2 + 273.15 #Celsius to Kelvin

"Converting total temperature to normal temperature"
def normTemp(TAT,M):
    T=TAT/(1+0.2*M**2)
    return T
T1 = TAT1 / (1+0.2*M1**2)
T2 = TAT2 / (1+0.2*M2**2)
T3 = TAT3 / (1+0.2*M3**2)
T4 = TAT4 / (1+0.2*M4**2)
T5 = TAT5 / (1+0.2*M5**2)
T6 = TAT6 / (1+0.2*M6**2)

"Calculate speed of sound"
def SpeedofSound(lambda_1,R,T):
    a=m.sqrt(lambda_1*R*T1)
    return a
a1 = m.sqrt(lambda_1*R*T1)
a2 = m.sqrt(lambda_1*R*T2)
a3 = m.sqrt(lambda_1*R*T3)
a4 = m.sqrt(lambda_1*R*T4)
a5 = m.sqrt(lambda_1*R*T5)
a6 = m.sqrt(lambda_1*R*T6)

"Calculate true airspeed"
def TAS(M,a):
    TAS=M*a
    return TAS
TAS1 = M1*a1
TAS2 = M2*a2
TAS3 = M3*a3
TAS4 = M4*a4
TAS5 = M5*a5
TAS6 = M6*a6

"Calculate density"
def density(p,R,T):
    rho=p/(R*T)
    return rho
rho1 = p1/(R*T1)
rho2 = p2/(R*T2)
rho3 = p3/(R*T3)
rho4 = p4/(R*T4)
rho5 = p5/(R*T5)
rho6 = p6/(R*T6)

"Caclulate equivalent airspeed"
def EAS(TAS,rho,rho0):
    EAS=TAS*m.sqrt(rho/rho0)
    return EAS
EAS1 = TAS1 * m.sqrt(rho1/rho0)
EAS2 = TAS2 * m.sqrt(rho2/rho0)
EAS3 = TAS3 * m.sqrt(rho3/rho0)
EAS4 = TAS4 * m.sqrt(rho4/rho0)
EAS5 = TAS5 * m.sqrt(rho5/rho0)
EAS6 = TAS6 * m.sqrt(rho6/rho0)

"Calculate lift coefficient for each measure point"
lift_coefficient_1 = lift1 / (0.5*rho1*TAS1**2*S)
lift_coefficient_2 = lift2 / (0.5*rho2*TAS2**2*S)
lift_coefficient_3 = lift3 / (0.5*rho3*TAS3**2*S)
lift_coefficient_4 = lift4 / (0.5*rho4*TAS4**2*S)
lift_coefficient_5 = lift5 / (0.5*rho5*TAS5**2*S)
lift_coefficient_6 = lift6 / (0.5*rho6*TAS6**2*S)
lift_coefficients = [lift_coefficient_1, lift_coefficient_2, lift_coefficient_3, 
                     lift_coefficient_4, lift_coefficient_5, lift_coefficient_6]

"Plot CL-alpha figure"

x = np.array(angle_of_attack)
y = np.array(lift_coefficients)
plt.plot(x, y, 'o', color='red', label='Measuring point')
mo, bo = np.polyfit(x, y, 1)
x_test = np.arange(-0.05,max(x),0.001)
plt.title(r'$C_{L} - \alpha$')
plt.xlabel(r'$\alpha$ [rad]')
plt.ylabel('$C_{L} [-]$')
plt.plot(x_test, mo*x_test + bo, color = 'black', label = 'Linear regression line')
plt.legend()


"Calculate intersection point for zero-lift angle of attack"
zero_lift_angle_of_attack = -bo/mo
print("Alpha_CL=0", zero_lift_angle_of_attack)

"Print CL_alpha and alpha_CL0"
print('CL_alpha:', mo)

"Calculate thrust"
Tp1 = 3678.46 + 3784.69
Tp2 = 3007.55 + 3069.61
Tp3 = 2410.92 + 2537.7
Tp4 = 1873.55 + 2026.55
Tp5 = 1903.04 + 2087.09
Tp6 = 2220.99 +2418.08
Tp=[Tp1,Tp2,Tp3,Tp4,Tp5,Tp6]

"Calculate drag"
def drag(angle,Thrust):
    Drag=[]
    for i in range(0,6):
        D=Tp[i]
        Drag.append(D)
    return(Drag)
D1 = Tp1
D2 = Tp2
D3 = Tp3
D4 = Tp4
D5 = Tp5
D6 = Tp6

"Caclulate drag coefficient"
def dragcoeffiecient(drag,rho,TAS,S):
    CD=drag/(0.5*rho*TAS**2*S)
    return CD
CD1 = D1 / (0.5*rho1*TAS1**2*S)
CD2 = D2 / (0.5*rho2*TAS2**2*S)
CD3 = D3 / (0.5*rho3*TAS3**2*S)
CD4 = D4 / (0.5*rho4*TAS4**2*S)
CD5 = D5 / (0.5*rho5*TAS5**2*S)
CD6 = D6 / (0.5*rho6*TAS6**2*S)
drag_coefficients = [CD1, CD2, CD3,CD4, CD5, CD6]

"Plot CD-alpha figure"
plt.figure()
r = np.array(angle_of_attack)
s = np.array(drag_coefficients)
plt.plot(r, s, 'o', color='red', label='Measuring point')
f, g, h = np.polyfit(r, s, 2)
r_test = np.arange(min(r),max(r),0.001)
plt.title(r'$C_{D} - \alpha$')
plt.xlabel(r'$\alpha$ [rad]')
plt.ylabel('$C_{D} [-]$')
plt.plot(r_test, h + g*r_test + f*r_test**2, color = 'black', label = '2nd order polynomial')

plt.legend()

"Plot CL-CD figure"
plt.figure()
b = np.array(drag_coefficients)
c = np.array(lift_coefficients)
plt.plot(b, c, 'o', color='red', label='Measuring point')
i, j, k = np.polyfit(b, c, 2)
b_test = np.arange(min(b),max(b),0.001)
plt.title('$C_{L} - C_{D}$')
plt.xlabel('$C_{D} [-]$')
plt.ylabel('$C_{L} [-]$')
plt.plot(b_test, k + j*b_test + i*b_test**2, color = 'black', label = '2nd order polynomial')
plt.legend()

"Lift coefficients squared"
lift_coefficients_squared = np.array(lift_coefficients)**2

"Plot CD - CL^2 figure"
plt.figure()
d = np.array(lift_coefficients_squared)
e = np.array(drag_coefficients)
plt.plot(d, e, 'o', color='red', label='Measuring point')
g, v = np.polyfit(d, e, 1)
d_test = np.arange(-0.05,max(d),0.001)
plt.title('$C_{L}^2 - C_{D}$')
plt.xlabel('$C_{L}^2 [-]$')
plt.ylabel('$C_{D} [-]$')
plt.plot(d_test, g*d_test + v, color = 'black', label = 'Linear regression line')
plt.show()
plt.legend()

"Print dCL^2_dCD and CD0"
print('dCL^2/dCD:', g)
print('CD0:', v)

"Calculate dynamic viscosity"
mu1 = mu0*(T1/T0)**(3/2)*(T0+S)/(T1+S)
mu2 = mu0*(T2/T0)**(3/2)*(T0+S)/(T2+S)
mu3 = mu0*(T3/T0)**(3/2)*(T0+S)/(T3+S)
mu4 = mu0*(T4/T0)**(3/2)*(T0+S)/(T4+S)
mu5 = mu0*(T5/T0)**(3/2)*(T0+S)/(T5+S)
mu6 = mu0*(T6/T0)**(3/2)*(T0+S)/(T6+S)

"Calculate Reynolds number"
Re1 = (rho1*TAS1*chord) / mu1 / 10**6
Re2 = (rho2*TAS2*chord) / mu2 / 10**6
Re3 = (rho3*TAS3*chord) / mu3 / 10**6
Re4 = (rho4*TAS4*chord) / mu4 / 10**6
Re5 = (rho5*TAS5*chord) / mu5 / 10**6
Re6 = (rho6*TAS6*chord) / mu6 / 10**6
reynolds_number = [Re1, Re2, Re3, Re4, Re5, Re6]

dCD_dCL2 = g #slope of CD-CL^2 plot 

"Calculate Oswald efficiency factor"
oswald = 1 / (m.pi * A * dCD_dCL2)
print('Oswald efficiency factor:', oswald)

"Calculate reduced equivalent airspeed"
def reducedV(rho,TAT, fuel_used,IAS,hp):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    weight_total = (OEW + weight_fuel + weight_payload - fuel_used) * g0
    V_c = (IAS - 2) * 0.514444444  # kts to m/s
    M = m.sqrt((2/0.4)*((1 + p0/p * ((1 + 0.4/2.8 *rho0/p0*V_c**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
    T = (TAT+273.15) / (1 + 0.2 * M ** 2)
    a = m.sqrt(lambda_1 * R * T)
    TAS = M * a
    EAS = TAS * m.sqrt(rho/rho0)
    Ve_bar = EAS * m.sqrt(standard_weight/weight_total)
    return Ve_bar
V_E1=reducedV(rho1,12.5,fuel_used1,249,hp1)
Ve_bar1 = EAS1 * m.sqrt(standard_weight / weight_total1)
print(V_E1, Ve_bar1)
Ve_bar2 = EAS2 * m.sqrt(standard_weight / weight_total2)
Ve_bar3 = EAS3 * m.sqrt(standard_weight / weight_total3)
Ve_bar4 = EAS4 * m.sqrt(standard_weight / weight_total4)
Ve_bar5 = EAS5 * m.sqrt(standard_weight / weight_total5)
Ve_bar6 = EAS6 * m.sqrt(standard_weight / weight_total6)

"Calculate standard thrust"
standard_thrust1 = Tp1 / (0.5 * rho1 * Ve_bar1 * 2 * engine_inlet_diameter)
standard_thrust2 = Tp2 / (0.5 * rho2 * Ve_bar2 * 2 * engine_inlet_diameter)
standard_thrust3 = Tp3 / (0.5 * rho3 * Ve_bar3 * 2 * engine_inlet_diameter)
standard_thrust4 = Tp4 / (0.5 * rho4 * Ve_bar4 * 2 * engine_inlet_diameter)
standard_thrust5 = Tp5 / (0.5 * rho5 * Ve_bar5 * 2 * engine_inlet_diameter)
standard_thrust6 = Tp6 / (0.5 * rho6 * Ve_bar6 * 2 * engine_inlet_diameter)

"Calculate temperature difference"
delta_t1 = T1 - (T0 + lambda0 * hp1)
delta_t2 = T2 - (T0 + lambda0 * hp2)
delta_t3 = T3 - (T0 + lambda0 * hp3)
delta_t4 = T4 - (T0 + lambda0 * hp4)
delta_t5 = T5 - (T0 + lambda0 * hp5)
delta_t6 = T6 - (T0 + lambda0 * hp6)