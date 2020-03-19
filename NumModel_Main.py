from math import *
import numpy as np
from matplotlib import pyplot as plt
import control as control
from control.matlab import lsim
from flight_data_reading import phi,p,r, alphaFL, thetaFL, qFL, V_tas2, alpha_0, theta_0,de,da,dr,times, time_eigenmotion, INIT_s, INIT_a, t_shaped, V0
from scipy.interpolate import griddata 
from StationaryMeasurements import weight_totalPH, rhoFL
#Asymmetric flight
pFL=p
rFL=r
phiFL=phi
#velocity  
V=V0         #INPUTS
V=104.7832256327374         #TAS
theta0=theta_0
#alpha0=0
#aircraft dimensions
S=30  #m^2
c = 2.0569  #m
b = 15.911 #m
rho=1.061626963061385#Weigh   #INPUT
W=63610.51533364605         #INPUT
m=W/9.80665
#m=13600*0.45359237
mju_b = m/(rho*S*b)#7000*2*V/b
mju_c = m/(rho*S*c)

#aircraft Inertia
Kxx2 = 0.019
Kyy2= 1.3925
Kzz2 = 0.042
Kxz = 0.002

#Aerodynamic coeffs
alpha0 = alpha_0             #INPUT!
CD_0 = 0.021271
CL_a = 4.8092325
e=0.76406
CL=2 * W / (rho * V ** 2 * S)    #INPUT 
#CD = CD_0 + (CL_a * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

#clean cruise (flaps up, gear up)
xcg=0.25*c

#Stability derivatives
#Assymetric
#Lateral force derivatives
CY_b=-0.75
CY_bdot = 0
CY_p = -0.0304
CY_r = 0.8495
CY_da = -0.04
CY_dr = 0.2300
#Roll moment derivatives
Cl_b = -0.1026
Cl_p = -0.71085
Cl_r = 0.2376
Cl_da = -0.23088
Cl_dr = 0.0344
#Yam moment derivatives
Cn_b = 0.1348
Cn_bdot = 0
Cn_p = -0.0602
Cn_r = -0.2061
Cn_da = -0.0120
Cn_dr = -0.0939

#Symmetric
CX_u = -0.095

CX_a = 0.4797
CZ_0 =-W*cos(theta0)/0.5/(rho*V**2*S)
CX_q = -0.2817
CZ_u = -0.37616
CZ_a = -5.7434
CZ_adot = -0.0035
CX_0 =W*sin(theta0)/0.5/(rho*V**2*S)
CZ_q = -5.6629
Cm_u = 0.0699
Cm_a = -0.5626              #INPUT
Cm_adot = 0.1780             #INPUT
Cm_q = -8.7941
CX_de=-0.0373
CZ_de=-0.6961 
Cm_de=-1.1642          #INPUT?
CX_dt=0
CZ_dt=0
Cm_dt=0

        
      #ASSYMETRIC
#Differential operators
Db=1

AssyM = np.matrix([[CY_b + (CY_bdot - 2*mju_b)*Db, CL, CY_p, CY_r - 4*mju_b],
                   [0, -0.5*Db, 1, 0],
                   [Cl_b, 0, Cl_p - 4*mju_b*Kxx2*Db, Cl_r + 4*mju_b*Kxz*Db],
                   [Cn_b + Cn_bdot*Db, 0, Cn_p + 4*mju_b*Kxz*Db, Cn_r - 4*mju_b*Kzz2*Db]])
#Angles = [beta,phi, p*b/2/V, r*b/2/V]

#AilRud = [da, dr]

#GETTING TO STATE SPACE MODEL
#xdot = [betadot, phidot, pdot, rdot]
#x[beta,phi,p,r]
#make dimensionless mass
#state space A=-C1^(-1)C2, B=-C1^(-1)C3

C1 = np.matrix([[(CY_bdot -2*mju_b)*b/V, 0, 0, 0],
                [0, -0.5*b/V, 0, 0],
                [0, 0, -4*0.5*b/V*Kxx2*mju_b*b/(V), 4*0.5*b/V*Kxz*mju_b*b/(V)],
                [Cn_bdot*b/V, 0, 4*0.5*b/V*Kxz*mju_b*b/(V), -4*0.5*b/V*Kzz2*mju_b*b/(V)]])
C2 = np.matrix([[CY_b, CL, b/2/V*CY_p, b/2/V*(CY_r - 4*mju_b)],
                [0, 0, b/2/V, 0],
                [Cl_b, 0, Cl_p*b/2/V, Cl_r*b/2/V],
                [Cn_b, 0, Cn_p*b/2/V, Cn_r*b/2/V]])
C3 = np.matrix([[-CY_da, -CY_dr],
                 [0,0],
                 [-Cl_da, -Cl_dr],
                 [-Cn_da, -Cn_dr]])
A=-np.linalg.inv(C1)*C2   
B=np.linalg.inv(C1)*C3
C=np.matrix([[1,0,0,0],
             [0,1,0,0],
             [0,0,1,0],
             [0,0,0,1]])
D=np.matrix([[0,0],
             [0,0],
             [0,0],
             [0,0]])
eigs=np.linalg.eigvals(A)

#time vector
T=time_eigenmotion    #choose time range for plotting
dt = 0.008       #INPUT
t=np.arange(0,T,dt)

sys_Ass=control.StateSpace(A,B,C,C3)
#T0, yout = control.initial_response(sys_Ass, t, X0=[0.5,0.5,0.5,0.5])
#T1,yout1=control.step_response(sys_Ass,t,[0,0,0,0])
'''plt.figure()
plt.plot(T1,yout[0])
plt.figure()
plt.plot(T1,yout[1])
plt.figure()
plt.plot(T1,yout[2])
plt.figure()
plt.plot(T1,yout[3])
plt.show()'''
                 #Symmetric

C1s = np.matrix([[-2*mju_c*c/V,0,0,0],
                [0,(CZ_adot-2*mju_c)*c/V,0,0],
                [0,0,-c/V,0],
                [0,Cm_adot*c/V,0,-2*mju_c*Kyy2*c/V*c/V]])
C2s = np.matrix([[CX_u, CX_a, CZ_0,0],
                [CZ_u, CZ_a, -CX_0, c/V*(CZ_q+2*mju_c)],
                [0, 0, 0, c/V],
                [Cm_u, Cm_a, 0, Cm_q*c/V]])
C3s = np.matrix([[CX_de],
                 [CZ_de],
                 [0],
                 [Cm_de]])

As=-np.linalg.inv(C1s)*C2s   
Bs=-np.linalg.inv(C1s)*C3s

Cs=np.matrix([[1,0,0,0],
             [0,1,0,0],
             [0,0,1,0],
             [0,0,0,1]])
Ds=np.matrix([[0,0],
             [0,0],
             [0,0],
             [0,0]])
sys_Sym=control.StateSpace(As,Bs,Cs,C3s)
eigsSYM = np.linalg.eigvals(As)

#inital vector INPUT
Xs = np.matrix([[INIT_s[0]],
                [INIT_s[1]],
                [INIT_s[2]],
                [INIT_s[3]]])
Xa = np.matrix([[INIT_a[0]],
                [INIT_a[1]],
                [INIT_a[2]],
                [INIT_a[3]]])

#de_list = griddata(times, de, t, method = "linear")

#interpolation for elevator deflection
de_list=[]
interpolants=[]
for i in range(1,len(t_shaped)):
    interpolant = (de[i]-de[i-1])/(t_shaped[i]-t_shaped[i-1])
    interpolants.append(interpolant)
de_list.append(de[0])
count=0
t[0]=t_shaped[0]
for j in range(1,len(t)):
    t[j]=dt+t[j-1]
    if count<time_eigenmotion*10-1:
        de_list.append(de[count]+interpolants[count]*(t[j]-t_shaped[count]))
    else:
        de_list.append(de[count]+interpolants[count-1]*(t[j]-t_shaped[count]))
    if count<time_eigenmotion*10-1 and t[j]>=t_shaped[count+1]:
       count=count+1
       
da_list=[]
interpolantsA=[]
for i in range(1,len(t_shaped)):
    interpolantA = (da[i]-da[i-1])/(t_shaped[i]-t_shaped[i-1])
    interpolantsA.append(interpolantA)
da_list.append(da[0])
countA=0
t[0]=t_shaped[0]
for j in range(1,len(t)):
    t[j]=dt+t[j-1]
    if countA<time_eigenmotion*10-1:
        da_list.append(da[countA]+interpolantsA[countA]*(t[j]-t_shaped[countA]))
    else:
        da_list.append(da[countA]+interpolantsA[countA-1]*(t[j]-t_shaped[countA]))
    if countA<time_eigenmotion*10-1 and t[j]>=t_shaped[countA+1]:
       countA=countA+1
dr_list=[]
interpolantsR=[]
for i in range(1,len(t_shaped)):
    interpolantR = (dr[i]-dr[i-1])/(t_shaped[i]-t_shaped[i-1])
    interpolantsR.append(interpolantR)
dr_list.append(dr[0])
countR=0
t[0]=t_shaped[0]
for j in range(1,len(t)):
    t[j]=dt+t[j-1]
    if countR<time_eigenmotion*10-1:
        dr_list.append(dr[countR]+interpolantsR[countR]*(t[j]-t_shaped[countR]))
    else:
        dr_list.append(dr[countR]+interpolantsR[countR-1]*(t[j]-t_shaped[countR]))
    if countR<time_eigenmotion*10-1 and t[j]>=t_shaped[countR+1]:
       countR=countR+1

u=[]
alpha=[]
theta=[]
q=[]

for i in range(len(t)):
    u.append(Xs[0,0])
    alpha.append(Xs[1,0])
    theta.append(Xs[2,0])
    q.append(Xs[3,0])
    
    Us = de_list[i]
    
    DXs = As*Xs + Bs*Us
    Xs = Xs + DXs*dt
    
beta=[]
phi=[]
p=[]
r=[]

for i in range(len(t)):
    beta.append(Xa[0,0])
    phi.append(Xa[1,0])
    p.append(Xa[2,0])
    r.append(Xa[3,0])
    
    Ua = np.matrix([[-da_list[i]],
                    [-dr_list[i]]])
    #print(Xa[1,0])
    DXa = A*Xa + B*Ua
    Xa = Xa + DXa*dt
#check   VERIFICATION!!!!
UA=np.zeros((len(da_list),2))
for rr in range(len(da_list)):
    UA[rr,0] = da_list[rr]
    UA[rr,1] = dr_list[rr]
yout,Ts,xout = control.matlab.lsim(sys_Sym, de_list, t, INIT_s)
youta,Ta,xouta = control.matlab.lsim(sys_Ass, UA, t, INIT_a)

plt.figure()
plt.title('simulation')
plt.plot(Ta,youta[:,1])
plt.show()

#a2=pi/180*(np.array(alpha))
#theta2=pi/180*(np.array(theta))

plt.figure()
plt.title('elevator')
plt.plot(t,de_list)
plt.figure()
plt.title('velocity')
plt.plot(t,u)
plt.plot(t_shaped,V_tas2)
plt.figure()
plt.title('roll, AoA')
plt.plot(t,alpha)
plt.plot(t_shaped,alphaFL)
plt.figure()
plt.title('RollRate, pitch')
plt.plot(t,theta)
plt.plot(t_shaped, thetaFL)
plt.figure()
plt.title('yawRate, pitchRate')
plt.plot(t,q)
plt.plot(t_shaped, qFL)
plt.show()




#VERIFICAITON

#simplified analytical eigenvalues for 
#SHORT PERIOD
Asp=2*mju_c*Kyy2*(2*mju_c - CZ_adot)
Bsp=-2*mju_c*Kyy2*CZ_a - (2*mju_c + CZ_q)*Cm_adot - (2*mju_c - CZ_adot)*Cm_q
Csp = CZ_a*Cm_q - Cm_a*(2*mju_c+CZ_q)

sp_real=(-Bsp/2/Asp)
sp_imag=(sqrt(4*Asp*Csp - Bsp**2)/2/Asp)

SPREAL = (sp_real*V/c)
SPIMAG = (sp_imag*V/c)
#PHUGOID
Aph = -4*mju_c**2
Bph = 2*mju_c*CX_u
Cph = -CZ_u*CZ_0

ph_real=(-Bph/2/Aph)
ph_imag=(sqrt(4*Aph*Cph - Bph**2)/2/Aph)
PHREAL=ph_real*V/c
PHIMAG = ph_imag*V/c           #absolute value

#PHUGOID2
Aph2=2*mju_c*(CZ_a*Cm_q - 2*mju_c*Cm_a)
Bph2 = 2*mju_c*(CX_u*Cm_a - Cm_u*CX_a) + Cm_q*(CZ_u*CX_a - CX_u*CZ_a)
Cph2 = CZ_0*(Cm_u*CZ_a - CZ_u*Cm_a)

ph2_real=(-Bph2/2/Aph2)
ph2_imag=(sqrt(4*Aph2*Cph2 - Bph2**2)/2/Aph2)
PH2REAL=ph2_real*V/c
PH2IMAG = ph2_imag*V/c       #absolute value

#A-PERIODIC ROLL
lambdaAPER = Cl_p/(4*mju_b*Kxx2)*V/b

#DUTCH ROLL
Adr = 8*mju_b**2*Kzz2
Bdr = -2*mju_b*(Cn_r+2*Kzz2*CY_b)
Cdr = 4*mju_b*Cn_b + CY_b*Cn_r

dr_real=(-Bdr/2/Adr)
dr_imag=(sqrt(4*Adr*Cdr - Bdr**2)/2/Adr)
drREAL=dr_real*V/b
drIMAG = dr_imag*V/b       #absolute value

#SPIRAL
lambdaSPIRAL = V/b*2*CL*(Cl_b*Cn_r - Cn_b*Cl_r)/(Cl_p*(CY_b*Cn_r + 4*mju_b*Cn_b) - Cn_p*(CY_b*Cl_r + 4*mju_b*Cl_b))


lambdaDUTCH2 = (Cn_b+Cn_r)/(4*mju_b*Kzz2)*V/b
