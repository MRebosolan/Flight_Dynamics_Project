from math import *
import numpy as np
from matplotlib import pyplot as plt
import control as control
import pandas as pd
#Asymmetric flight
from scipy.io import loadmat
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

#For plotting
t_DutchRoll = 3164  #CHECKED 3164
t_spiral = 3280    #tbd
t_aperiodic = 3040   #tbd
t_phugoid = 2812   #CHECKED 2812
t_shortper = 2982   #CHECKED2982  tbd..


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

current_index = phugoid_index  #INPUT        
if current_index == phugoid_index:
    time_eigenmotion = 150         #INPUT
if current_index == shortper_index or current_index == spiral_index or current_index==DutchRoll_index:
    time_eigenmotion = 30 #actually 8 and 15

de_raw = de_FL[current_index: current_index + time_eigenmotion*10]        
de = pi/180*de_raw.reshape(time_eigenmotion*10) #- de_raw[0]*pi/180

da = pi/180*da_FL[current_index: current_index + time_eigenmotion*10] 
dr = pi/180*dr_FL[current_index: current_index + time_eigenmotion*10] 
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
tat = tat_FL[current_index]
#aileron, rudder deflections /these lines mby not needed
aileron = pi/180*da.reshape(time_eigenmotion*10) - da[0]*pi/180
rudder = pi/180*dr.reshape(time_eigenmotion*10) - dr[0]*pi/180
elevator = pi/180*de.reshape(time_eigenmotion*10) - de[0]*pi/180
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
plt.plot(t_shaped, dr)
plt.plot(t_shaped, da)


V_tas2=(V_tas-V_tas[0])/V_tas[0]
#plt.figure()
#plt.plot(t_shaped, beta)
#plt.show()
plt.figure()
plt.plot(t_shaped, phi)
plt.show()
plt.figure()
plt.plot(t_shaped, p)
plt.show()
plt.figure()
plt.plot(t_shaped, r)
plt.show()'''
#Array = Flightdata['flightdata'][0][0]['Ahrs1_bYawRate']

'''
mju_c = 7000   ##INPUT ?
#velocity
V=82.82
#aircraft dimensions
S=30  #m^2
C = 2.0569  #m
b = 15.911 #m
rho=1.01 #SI
#Weight
mju_c = 6800/(rho*S*b)#7000*2*V/b   ##INPUT ?


#aircraft Inertia
Kxx2 = 0.019
Kyy2= 1.3925
Kzz2 = 0.042
Kxz = 0.002

#Aerodynamic coeffs
CD_0 = 0.04
CL_a = 5.084
e=0.8
CL=0.6    #INPUT ?

#clean cruise (flaps up, gear up)
xcg=0.25*C

#Assymetric
#Lateral force derivatives
CX_u = -0.0279
mju_c =
Dc =
CX_a = -0.4797
CZ_0 =
CX_q = -0.2817
CZ_u = -0.3762
CZ_a = -5.7434
CZ_adot = -0.0035
CX_0 =
CZ_q = -5.6629
CM_u = 0.0699
CM_a = -0.5626
CM_adot = 0.1780
CM_q = -8.7941




#Differential operators
Dc = 1

symM = np.matrix([[CX_u -2*mju_c*Dc, CX_a, CZ_0, CX_q],
                   [CZ_u, CZ_a + (CZ_adot-2*mju_c)*Dc, -CX_0, CZ_q + 2*mju_c],
                   [0, 0, -Dc, 1],
                   [CM_u, CM_a + (CM_adot*Dc), 0, CM_q - 2*(mju_c*Kyy2*Dc)]])
#Angles = [u,alpha,theta,q*c/V]
RHS = np.matrix([[-CX_d],
                 [-CZ_d],
                 [0],
                 [-CM_d]])
#AilRud = [de]

#GETTING TO STATE SPACE MODEL
#xdot = [udot, adot, thetadot, qdot]
#x[u,a,theta,q]

#x[beta,phi,p,r]

#make dimensionless mass
#state space A=-C1^(-1)C2, B=-C1^(-1)C3

C1 = np.matrix([[-2*mju_c*c/V,0,0,0],
                [0,(CZ_adot-2*mju_c)*c/V,0,0],
                [0,0,-u/V,0],
                [0,CM_adot*u/V,0,-2*mju_c*Kyy2*u/V]])
C2 = np.matrix([[CY_b, CL, b/2/V*CY_p, b/2/V*(CY_r - 4*mju_b)],
                [0, 0, b/2/V, 0],
                [Cl_b, 0, Cl_p*b/2/V, Cl_r*b/2/V],
                [Cn_b, 0, Cn_p*b/2/V, Cn_r*b/2/V]])

C3 = -RHS
A=-np.linalg.inv(C1)*C2
B=-np.linalg.inv(C1)*C3
C=np.matrix([[1,0,0,0],
             [0,1,0,0],
             [0,0,1,0],
             [0,0,0,1]])
D=np.matrix([[0,0],
             [0,0],
             [0,0],
             [0,0]])
D=RHS
eigs=np.linalg.eigvals(A)
t=np.zeros(200)
for i in range(1,200):
    t[i]=t[i-1]+0.1
sys=control.StateSpace(A,B,C,D)
T0, yout = control.initial_response(sys, t, X0=[0.5,0.5,0.5,0.5])
T1,yout1=control.step_response(sys,t,[0,0,0,0])
plt.plot(T1,yout[3])
plt.show()
'''
