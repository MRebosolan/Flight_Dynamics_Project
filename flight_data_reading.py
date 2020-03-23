from math import *
import numpy as np
from matplotlib import pyplot as plt
import control as control
import pandas as pd
#Asymmetric flight
from scipy.io import loadmat
from scipy.signal import find_peaks
S=30  #m^2
c = 2.0569  #m
b = 15.911 #m
p0 = 101325 #pa
rho0 = 1.225 #kg/m3
T0 = 288.15 #K
lambda_1 = 1.4 #ratio of heats
lambda0 = -0.0065 #K/m
R = 287 #J/kg*K
mu0 = 1.7894*10**(-5) #kg/m*s
OEW = 9165 * 0.45359237 #lbs to kg     CHECK!!!weights
weight_fuel = 4050 * 0.45359237 #lbs to kg
weight_payload = 695 #kg 
g0 = 9.80665


weight_payloadFL = 733 #kg
weight_fuelEIG = 2750*0.45359237 #lbs to kg

Flightdata = loadmat('FTISxprt-20200310_flight4.mat')
#Weight
data = [[row.flat[0] for row in line] for line in Flightdata['flightdata'][0]]

columns = ['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'column_SE', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2', 'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate', 'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas', 'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time' ]
dftrain = pd.DataFrame(data, columns=columns)
DATA = dftrain.to_numpy()    #converting data
#recovering data
time = np.array(DATA[0][48][0])
timesFL=time[0,:]   #time array

deflight = np.array(DATA[0][17][0])
de_FL = deflight[:,0]

daflight = np.array(DATA[0][16][0])
da_FL = daflight[:,0]

drflight = np.array(DATA[0][18][0])
dr_FL = drflight[:,0]

alphaflight = np.array(DATA[0][0][0])
alpha_FL = alphaflight[:,0]

thetaflight = np.array(DATA[0][22][0])
theta_FL = thetaflight[:,0]


phiflight = np.array(DATA[0][21][0])
phi_FL = phiflight[:,0]

qflight = np.array(DATA[0][27][0])
q_FL = qflight[:,0]

altflight = np.array(DATA[0][37][0])
alt_FL = altflight[:,0]*0.3048

Machflight = np.array(DATA[0][40][0])
Mach_FL = Machflight[:,0]

tatflight = np.array(DATA[0][36][0])
tat_FL = tatflight[:,0]

V_tasflight = np.array(DATA[0][42][0])
V_tas_FL = V_tasflight[:,0]*0.514444

pflight = np.array(DATA[0][26][0])
p_FL = pflight[:,0]

rflight = np.array(DATA[0][28][0])
r_FL = rflight[:,0]

fuel_right_usedflight = np.array(DATA[0][15][0])
fuel_right_used_FL = fuel_right_usedflight[:,0]

fuel_left_usedflight = np.array(DATA[0][14][0])
fuel_left_used_FL = fuel_left_usedflight[:,0]

#For plotting
t_DutchRoll = 3163  #CHECKED 3164
t_spiral = 3315-30    #tbd
t_aperiodic = 3091-4   #tbd
t_phugoid = 2812   #CHECKED 2815 for the data
t_shortper = 2979-2   #


for i in range(20000,len(timesFL)):
    if abs(t_shortper - timesFL[i]) < 0.12:
        shortper_index = i
    if abs(t_phugoid - timesFL[i]) < 0.12:
        phugoid_index = i
    if abs(t_aperiodic - timesFL[i]) < 0.12:
        aperiodic_index = i
    if abs(t_DutchRoll - timesFL[i]) < 0.12:
        DutchRoll_index = i
    if abs(t_spiral - timesFL[i]) < 0.12:
        spiral_index = i

#Lengths of eigenmotions (approx)   currently RANDOM numbers
#t_lenDutchRoll =
#t_lenspiral = 
#t_lenaperiodic = 
#t_lenshortper = 
#t_lenphugoid = 
#adjust time interval depending on current index
######################################
current_index = DutchRoll_index  #INPUT
#################################        
if current_index == phugoid_index or current_index==spiral_index:
    time_eigenmotion = 150  #150        #INPUT
if current_index == shortper_index:
    time_eigenmotion = 8
if current_index==DutchRoll_index or current_index==aperiodic_index:
    time_eigenmotion = 19

de_raw = de_FL[current_index: current_index + time_eigenmotion*10]        
de = pi/180*(de_raw - de_raw[0]) #- de_raw[0]*pi/180

da_raw = da_FL[current_index: current_index + time_eigenmotion*10] 
dr_raw = dr_FL[current_index: current_index + time_eigenmotion*10] 

da=pi/180*(da_raw-da_raw[0])
dr=pi/180*(dr_raw - dr_raw[0])

times = timesFL[current_index: current_index + time_eigenmotion*10]

alpha_raw = alpha_FL[current_index: current_index + time_eigenmotion*10]
alphaFL = pi/180*(alpha_raw - alpha_raw[0])
alpha_0 = alpha_raw[0]*pi/180

theta_raw = theta_FL[current_index: current_index + time_eigenmotion*10]
thetaFL = pi/180*(theta_raw - theta_raw[0])
theta_0 = theta_raw[0]*pi/180

q_raw = q_FL[current_index: current_index + time_eigenmotion*10]
qFL = pi/180*(q_raw-q_raw[0])

phi_raw = phi_FL[current_index: current_index + time_eigenmotion*10]
#phi = pi/180*(phi_raw.reshape(time_eigenmotion*10))
phi = pi/180*(phi_raw - phi_raw[0])

p_raw = p_FL[current_index: current_index + time_eigenmotion*10]
p = pi/180*(p_raw - p_raw[0])
#p = pi/180*(p_raw.reshape(time_eigenmotion*10))

r_raw = r_FL[current_index: current_index + time_eigenmotion*10]
r = pi/180*(r_raw - r_raw[0])

#r = pi/180*(r_raw.reshape(time_eigenmotion*10))

V_tas_raw = V_tas_FL[current_index: current_index + time_eigenmotion*10]
V_tas = (V_tas_raw.reshape(time_eigenmotion*10))

V0 = V_tas_FL[current_index]
alt = alt_FL[current_index]
Mach = Mach_FL[current_index]
tat = tat_FL[current_index]+273.15
pressureFL = p0 * (1 + (lambda0*alt)/T0)**(-g0/(lambda0*R))
TFL = tat/ (1+0.2*Mach**2)
dynamic_viscosity = mu0*(TFL/T0)**(3/2)*((T0+S)/(TFL+S))
rhoFL = pressureFL/R/TFL

ReFL=rhoFL*V0*c/dynamic_viscosity

fuel_used_FL = 0.45359237*(fuel_right_used_FL[current_index] + fuel_left_used_FL[current_index]) 
weight_totalFL = (OEW + weight_fuelEIG + weight_payloadFL - fuel_used_FL) * g0

#aileron, rudder deflections /these lines mby not needed
aileron = pi/180*da.reshape(time_eigenmotion*10) - da[0]*pi/180
rudder = pi/180*dr.reshape(time_eigenmotion*10) - dr[0]*pi/180
elevator = pi/180*de.reshape(time_eigenmotion*10) - de[0]*pi/180

#aircraft dimensions
S=30  #m^2
c = 2.0569  #m
b = 15.911 #m
rho=rhoFL#Weigh   #INPUT
W=weight_totalFL        #INPUT
m=W/9.80665
#m=13600*0.45359237
mju_b = m/(rho*S*b)#7000*2*V/b
mju_c = m/(rho*S*c)

#aircraft Inertia
Kxx2 = 0.019
Kyy2= 1.3925
Kzz2 = 0.042
Kxz = 0.002

V=V0
#Initial vector

u_init = 0*V0 
alpha_init = alphaFL[0]
theta_init = thetaFL[0]
q_init = 0*q_raw[0]*pi/180

INIT_s = [u_init, alpha_init, theta_init, q_init]


#Beta for spiral non-zero????
beta_init =0
phi_init = 0*phi_raw[0]*pi/180
p_init = 0*p_raw[0]*pi/180
r_init = 0*r_raw[0]*pi/180
INIT_a = [beta_init, phi_init, p_init, r_init]

t_shaped = (times - times[0]).reshape(time_eigenmotion*10)
V_tas2=(V_tas-V_tas[0])/V_tas[0]
'''
plt.figure()
plt.title('elevator/rudder/aileron')
plt.plot(t_shaped, de)
'''
for i in range(2,len(V_tas2)):
    if V_tas2[i]>0 and V_tas2[i]<V_tas2[i-1] and V_tas2[i-1]>V_tas2[i-2]:
        #print(V_tas2[i-1])
        es=1
'''
plt.figure()
plt.xlabel("time(s)")
plt.ylabel("u")
plt.plot(t_shaped, thetaFL)
plt.show()
'''
'''
#plt.figure()
#plt.plot(t_shaped, alphaFL)
#plt.show()
#plt.figure()
#plt.plot(t_shaped, thetaFL)
#plt.figure()
#plt.plot(t_shaped, qFL)
#plt.show()
#Array = Flightdata['flightdata'][0][0]['Ahrs1_bYawRate']
'''
#which function to analyse?
'''
Function=r

if Function[10]==thetaFL[10] or Function[10]==phi[10] or Function[10]==p[10] or Function[10]==r[10]:
    corr = 1
    if Function[10]==thetaFL[10]:
        height1=0.06
    else:
        height1=0.006
else:
    corr=0
    height1=0.01
peaks = find_peaks(Function, height=height1, distance=30)
fit_x = np.zeros(len(peaks[0])-corr)
damping = np.zeros(len(peaks[0])-corr)
for i in range(corr,len(peaks[0])):
    fit_x[i-corr] = t_shaped[peaks[0][i]]
    damping[i-corr] = Function[peaks[0][i]]

#Get frequency, period and Damping ratio for flight data and then deduce parameters and eigs
damp_funct=np.polyfit(fit_x, np.log(damping), 1)
#function in form Ae^(bx)
damp_A = exp(damp_funct[1])
damp_b = damp_funct[0]   #real part of the eigenvalue!

period = (fit_x[-1]-fit_x[0])/(len(fit_x)-1)
T_12 = log(0.5)/damp_b

eig_real = damp_b
eig_imag = 2*pi/period

wn = 2*pi/period
w0=sqrt(eig_real**2 + eig_imag**2)
damp_ratio = -eig_real/w0
Cn_bFL = w0**2*2*mju_b*Kzz2*(b/V)**2
Cn_rFL = -damp_ratio*4*sqrt(2*mju_b*Kzz2*Cn_bFL)

ytab=[]
xtab=[]
s=0
duration=time_eigenmotion*10
for i in range(duration):
    xtab.append(s)
    y=damp_A*exp(damp_b*s)
    ytab.append(y)
    s=s+0.1
#plt.plot(xtab,ytab)
#plt.show()'''

