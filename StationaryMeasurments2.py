#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"Created on Fri Mar  6 15:58:53 2020"
import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import csv

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


title="20200310_V4.xlsx" #write import file
file = pd.read_excel(title)
file1 = file.to_numpy()
def valueimport2(row1,row2,col,file2):
    for i in range(row1,row2+1):
        file2.append(float(file1[i][col]))
    return file2

angle_of_attack=[]
angle_of_attack=valueimport2(26,31,5,angle_of_attack)
"Transform angle of attack to radians"
def degtorad(angledeg):
    return angledeg*m.pi/180
for i in range(6):
    angle_of_attack[i]=degtorad(angle_of_attack[i])
alpha_1 = angle_of_attack[0]
alpha_2 = angle_of_attack[1]
alpha_3 = angle_of_attack[2]
alpha_4 = angle_of_attack[3]
alpha_5 = angle_of_attack[4]
alpha_6 = angle_of_attack[5]

def lbstokg(lbs):
    return lbs*0.45359237
OEW = lbstokg(9165)
weight_fuel=[]
weight_fuel = lbstokg(sum((valueimport2(16,16,3,weight_fuel))))
weight_payload=[]
weight_payload=sum(valueimport2(6,14,7,weight_payload))
g0 = 9.80665

fuel_used=[]
fuel_used=valueimport2(26,31,8,fuel_used)
for i in range(6):
    fuel_used[i]=lbstokg(fuel_used[i])

"Calculate weight"
def total_weight(fuel_used):
    return (OEW + weight_fuel + weight_payload - fuel_used) * g0
weight_total=[]
for i in fuel_used:
    weight_total.append(total_weight(i))
print(weight_total[5])
weight_total1 = weight_total[0]
weight_total2 = weight_total[1]
weight_total3 = weight_total[2]
weight_total4 = weight_total[3]
weight_total5 = weight_total[4]
weight_total6 = weight_total[5]

"Calculate lift"
def lift(weight):
    return weight
lift=weight_total
print(lift)
lift1 = weight_total1
lift2 = weight_total2
lift3 = weight_total3
lift4 = weight_total4
lift5 = weight_total5
lift6 = weight_total6

"Transform IAS to CAS"
IAS=[]
IAS=valueimport2(26,31,4,IAS)
def IAStoCAS(IAS):
    V_c=(IAS-2)*0.514444444
    return(V_c)
V_c1 = (249 - 2) * 0.514444444 #kts to m/s
V_c2 = (219 - 2) * 0.514444444 #kts to m/s
V_c3 = (193 - 2) * 0.514444444 #kts to m/s
V_c4 = (160 - 2) * 0.514444444 #kts to m/s
V_c5 = (131 - 2) * 0.514444444 #kts to m/s
V_c6 = (113 - 2) * 0.514444444 #kts to m/s

"Transform height to m"
hp=[]
valueimport2(26,31,3,hp)
def fttom(hp):
    return hp*0.3048
for i in range(6):
    hp[i]=fttom(hp[i])
hp=np.array(hp)
hp1 = 6990 * 0.3048 #ft to m
hp2 = 6990 * 0.3048 #ft to m
hp3 = 7000 * 0.3048 #ft to m
hp4 = 6960 * 0.3048 #ft to m
hp5 = 6950 * 0.3048 #ft to m
hp6 = 7000 * 0.3048 #ft to m
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
def mach(p,V_c):
    M=m.sqrt((2/(1-lambda_1))*((1 + p0/p1 * ((1 + (1-lambda_1)/(2*lambda_1)*rho0/p0*V_c1**2)**(lambda_1/(1-lambda_1))-1))**((1-lambda_1)/lambda_1) -1))
    return M
M1 = m.sqrt((2/0.4)*((1 + p0/p1 * ((1 + 0.4/2.8 *rho0/p0*V_c1**2)**(1.4/0.4)-1))**(0.4/1.4) -1))

M2 = m.sqrt((2/0.4)*((1 + p0/p2 * ((1 + 0.4/2.8 *rho0/p0*V_c2**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M3 = m.sqrt((2/0.4)*((1 + p0/p3 * ((1 + 0.4/2.8 *rho0/p0*V_c3**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M4 = m.sqrt((2/0.4)*((1 + p0/p4 * ((1 + 0.4/2.8 *rho0/p0*V_c4**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M5 = m.sqrt((2/0.4)*((1 + p0/p5 * ((1 + 0.4/2.8 *rho0/p0*V_c5**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
M6 = m.sqrt((2/0.4)*((1 + p0/p6 * ((1 + 0.4/2.8 *rho0/p0*V_c6**2)**(1.4/0.4)-1))**(0.4/1.4) -1))
print(M1, M2, M3, M4, M5, M6)
"Measured total temperature"
T_org=[]
T_org=valueimport2(26,31,9,T_org)

def TemptoK(T):
    TAT=T+273.15
    return TAT
TAT1 = 13.2 + 273.15 #Celsius to Kelvin
TAT2 = 11.2 + 273.15 #Celsius to Kelvin
TAT3 = 9.2 + 273.15 #Celsius to Kelvin
TAT4 = 7.8 + 273.15 #Celsius to Kelvin
TAT5 = 6.5 + 273.15 #Celsius to Kelvin
TAT6 = 5.8 + 273.15 #Celsius to Kelvin

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
    a=m.sqrt(gamma*R*T1)
    return a
a1 = m.sqrt(gamma*R*T1)
a2 = m.sqrt(gamma*R*T2)
a3 = m.sqrt(gamma*R*T3)
a4 = m.sqrt(gamma*R*T4)
a5 = m.sqrt(gamma*R*T5)
a6 = m.sqrt(gamma*R*T6)

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
print(rho1)

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
for i in range(6):
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
for i in range(6):
    delta_t.append(deltat(T_org[i],hp[i],IAS[i]))
print(delta_t[0])

"getting all inputs for the thrust file"
def mach(hp,IAS) :
    p = p0 * (1 + (lambda0 * hp) / T0) ** (-g0 / (lambda0 * R))
    V_c = (IAS - 2) * 0.514444444
    M = m.sqrt((2/(gamma-1))*((1+(p0/p)*((1+(gamma-1)/(2*gamma)*(rho0/p0)*V_c**2)**(gamma/(gamma-1))-1))**((gamma-1)/gamma)-1))
    return M
M=[]
for i in range(6):
    M.append(mach(hp[i],IAS[i]))
FFl=[]
valueimport2(26,31,6,FFl)
FFr=[]
valueimport2(26,31,7,FFr)
matlab=[]
def lbshrtokgs(FuelFlow):
    return lbstokg(FuelFlow)/3600

for i in range(6):
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
for i in range(6):
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
for i in range(6):
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
for i in range(6):
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
for i in range(6):
    standard_thrust.append(stanthrust(thrust[i],hp[i],IAS[i],T_org[i],fuel_used[i]))

