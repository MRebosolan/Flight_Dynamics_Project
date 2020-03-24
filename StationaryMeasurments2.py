#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"Created on Fri Mar  6 15:58:53 2020"
import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
import subprocess


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
gamma = 1.4 #ratio of heats
lambda0 = -0.0065 #K/m
lambda_r=2.6 #bypass ratio
R = 287 #J/kg*K
mu0 = 1.7894*10**(-5) #kg/m*s
standard_weight = 60500 #Standard weight [N]
engine_inlet_diameter = 0.686 #[m2]
"inputs"
rowbegin=26
rowend=31
lenrow=rowend-rowbegin+1

title="20200310_V4.xlsx" #write import file
file = pd.read_excel(title)
file1 = file.to_numpy()
def valueimport2(row1,row2,col,file2):
    for i in range(row1,row2+1):
        file2.append(float(file1[i][col]))
    return file2

angle_of_attack=[]
angle_of_attack=valueimport2(rowbegin,rowend,5,angle_of_attack)
"Transform angle of attack to radians"
def degtorad(angledeg):
    return angledeg*m.pi/180
for i in range(lenrow):
    angle_of_attack[i]=degtorad(angle_of_attack[i])


def lbstokg(lbs):
    return round(lbs*0.45359237,7)
OEW = lbstokg(9165)
weight_fuel=[]
weight_fuel = lbstokg(sum((valueimport2(16,16,3,weight_fuel))))
weight_payload=[]
weight_payload=sum(valueimport2(6,14,7,weight_payload))
g0 = 9.80665

fuel_used=[]
fuel_used=valueimport2(rowbegin,rowend,8,fuel_used)
for i in range(lenrow):
    fuel_used[i]=lbstokg(fuel_used[i])

"Calculate weight"
def total_weight(fuel_used):
    return (OEW + weight_fuel + weight_payload - fuel_used) * g0
weight_total=[]
for i in fuel_used:
    weight_total.append(total_weight(i))

"Transform IAS to CAS"
IAS=[]
IAS=valueimport2(rowbegin,rowend,4,IAS)

"Transform height to m"
hp=[]
valueimport2(rowbegin,rowend,3,hp)
def fttom(hp):
    return hp*0.3048
for i in range(lenrow):
    hp[i]=fttom(hp[i])
hp=np.array(hp)

"Measured total temperature"
T_org=[]
T_org=valueimport2(rowbegin,rowend,9,T_org)

"density calculation"
def density(hp,IAS,T_org):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444  # kts to m/s
    M = m.sqrt(
        (2 / 0.4) * ((1 + p0 / p * ((1 + 0.4 / 2.8 * rho0 / p0 * V_c ** 2) ** (1.4 / 0.4) - 1)) ** (0.4 / 1.4) - 1))
    T = (T_org + 273.15) / (1 + 0.2 * M ** 2)
    rho = p / (R * T)
    return rho

"Calculate lift coefficient for each measure point"
def lift_coefficient(fuel_used,hp,IAS,T_org) :
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444
    M = m.sqrt((2/(gamma-1))*((1+(p0/p)*((1+(gamma-1)/(2*gamma)*(rho0/p0)*V_c**2)**(gamma/(gamma-1))-1))**((gamma-1)/gamma)-1))
    TAT = T_org + 273.15
    T = TAT / (1 + 0.2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    TAS = M*a
    rho = p / (R * T)
    lift=(OEW + weight_fuel + weight_payload - fuel_used) * g0
    return lift/(0.5*rho*S*TAS**2)
lift_coefficients=[]
for i in range(lenrow):
    lift_coefficients.append(lift_coefficient(fuel_used[i],hp[i],IAS[i],T_org[i]))
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

"Calculate temperature difference"
delta_t=[]
def deltat(T_org,hp,IAS):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444
    M = m.sqrt((2 / (gamma - 1)) * ((1 + (p0 / p) * (
                (1 + (gamma - 1) / (2 * gamma) * (rho0 / p0) * V_c ** 2) ** (gamma / (gamma - 1)) - 1)) ** (
                                                (gamma - 1) / gamma) - 1))
    TAT = T_org + 273.15
    T = TAT / (1 + 0.2 * M ** 2)
    return T - (T0 + lambda0 * hp)
for i in range(lenrow):
    delta_t.append(deltat(T_org[i],hp[i],IAS[i]))


"getting all inputs for the thrust file"
def mach(hp,IAS) :
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444
    M = m.sqrt((2/(gamma-1))*((1+(p0/p)*((1+(gamma-1)/(2*gamma)*(rho0/p0)*V_c**2)**(gamma/(gamma-1))-1))**((gamma-1)/gamma)-1))
    return M
M=[]
for i in range(lenrow):
    M.append(mach(hp[i],IAS[i]))
FFl=[]
valueimport2(rowbegin,rowend,6,FFl)
FFr=[]
valueimport2(rowbegin,rowend,7,FFr)
matlab=[]
def lbshrtokgs(FuelFlow):
    return lbstokg(FuelFlow)/3600

for i in range(lenrow):
    lol=hp[i], M[i], delta_t[i], lbshrtokgs(FFl[i]), lbshrtokgs(FFr[i])
    matlab.append(lol)

f = open("matlab.dat", 'w+')
for row in matlab:
    for item in row:
        f.write(str(item) + " ")
    f.write("\n")#Nieuwe regel
f.close()
subprocess.Popen([r"C:\Users\lizzy\OneDrive\Documenten\Universiteit\2019-2020\SVV\Flight_Dynamics_Project\thrust(1).exe"])
data = np.genfromtxt("thrust.dat")
thrust=[]
for i in range(lenrow):
    thrust.append(sum(data[i]))

"Caclulate drag coefficient"
def dragcoeffiecient(thrust,hp,IAS,T_org):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444
    M = m.sqrt((2 / (gamma - 1)) * ((1 + (p0 / p) * (
                (1 + (gamma - 1) / (2 * gamma) * (rho0 / p0) * V_c ** 2) ** (gamma / (gamma - 1)) - 1)) ** (
                                                (gamma - 1) / gamma) - 1))
    TAT = T_org + 273.15
    T = TAT / (1 + 0.2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    TAS = M * a
    rho = p / (R * T)
    return thrust/(0.5*rho*S*TAS**2)
drag_coefficients=[]
for i in range(lenrow):
    drag_coefficients.append(dragcoeffiecient(thrust[i],hp[i],IAS[i],T_org[i]))

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
plt.legend()

"Print dCL^2_dCD and CD0"
print('dCL^2/dCD:', g)
print('CD0:', v)

"Calculate dynamic viscosity"
def reynoldsnumber(T_org,hp,IAS):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444  # kts to m/s
    M = m.sqrt((2 / 0.4) * ((1 + p0 / p * ((1 + 0.4 / 2.8 * rho0 / p0 * V_c ** 2) ** (1.4 / 0.4) - 1)) ** (0.4 / 1.4) - 1))
    T = (T_org + 273.15) / (1 + 0.2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    rho = p / (R * T)
    TAS = M * a
    rho = p / (R * T)
    mu=mu0*(T/T0)**(3/2)*(T0+S)/(T+S)
    Re=(rho * TAS * chord) / mu / 10 ** 6
    return mu, Re
Dynamic_viscolity=[]
Reynolds_number=[]
for i in range(lenrow):
    mu, Re =reynoldsnumber(T_org[i],hp[i],IAS[i])
    Dynamic_viscolity.append(mu)
    Reynolds_number.append(Re)


dCD_dCL2 = g #slope of CD-CL^2 plot

"Calculate Oswald efficiency factor"

oswald = 1 / (m.pi * A * dCD_dCL2)
print('Oswald efficiency factor:', oswald)

"Calculate reduced equivalent airspeed"
def reducedV(T_org, fuel_used,IAS,hp):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    weight_total = total_weight(fuel_used)
    V_c = (IAS - 2) * 0.514444444  # kts to m/s
    M = m.sqrt((2/0.4)*((1 + p0/p * ((1 + 0.4/2.8 *rho0/p0*V_c**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
    T = (T_org+273.15) / (1 + 0.2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    rho = p / (R * T)
    TAS = M * a
    EAS = TAS * m.sqrt(rho/rho0)
    Ve_bar = EAS * m.sqrt(standard_weight/weight_total)
    return Ve_bar
Ve_bar=[]
for i in range(lenrow):
    Ve_bar.append(reducedV(T_org[i],fuel_used[i],IAS[i],hp[i]))

"Calculate standard thrust"
def stanthrust(thrust,hp,IAS,T_org,fuel_used):
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444  # kts to m/s
    M = m.sqrt(
        (2 / 0.4) * ((1 + p0 / p * ((1 + 0.4 / 2.8 * rho0 / p0 * V_c ** 2) ** (1.4 / 0.4) - 1)) ** (0.4 / 1.4) - 1))
    T = (T_org + 273.15) / (1 + 0.2 * M ** 2)
    rho = p / (R * T)
    Ve_bar=reducedV(T_org, fuel_used,IAS,hp)
    return thrust / (0.5 * rho * Ve_bar * 2 * engine_inlet_diameter)
standard_thrust=[]
for i in range(lenrow):
    standard_thrust.append(stanthrust(thrust[i],hp[i],IAS[i],T_org[i],fuel_used[i]))

