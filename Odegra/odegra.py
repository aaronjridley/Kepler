#!/opt/local/bin/python

## Import needed modules
import matplotlib.pyplot as plt
from pylab import *
from scipy.integrate import odeint
import re
import sys

# Set some constants:

G  = 6.67384e-11 # Gravitational Constants in m3 kg-1 s-2

# Sun-Earth Distance
AU = 149597870.0 * 1000.0 # m
Re = 6372.0 * 1000.0

dt = 60.0

StartTime = datetime.datetime(2021, 9, 22, 0, 0, 0)
StartThrust = datetime.datetime(2025, 1, 1, 12, 0, 0)
EndThrust = datetime.datetime(2025, 3, 31, 23, 59, 59)

start_thrust_time = (StartThrust - StartTime).total_seconds()
end_thrust_time = (EndThrust - StartTime).total_seconds()

mass = 1.0e15 # kg
thrust = 1.0e7 # Newtons

def propagate(StartPos, StartVel, TotalTime):

    t0 = 0.
    # create time array starting at t0, 
    t = arange(t0, TotalTime+np.abs(dt), np.abs(dt))   
    if (dt < 0):
        t = -t
    Time = []
    for ti in t:
        Time.append(StartTime + datetime.timedelta(seconds=ti))

    # Solve the ODE with odeint returns an 2-dimensional array with
    # the first index specifying the time and the second index
    # specifying the component of the state vector

    X0 = np.concatenate((np.array(StartPos), np.array(StartVel)))
    
    X = odeint(f_func, X0, t)

    x = X[:,0]/AU
    y = X[:,1]/AU  # x, y, vx, and vy are arrays that specify their 
    z = X[:,2]/AU

    vx = X[:,3] # values at any given time index
    vy = X[:,4]
    vz = X[:,5]

    position = [x, y, z]
    velocity = [vx, vy, vz]

    return position, velocity, Time
    
def set_initial_conditions(Apogee, Perigee, Rotation, Inc = 0.0):
    
    SemiMajor = (Perigee + Apogee) / 2.0
    Vsma = sqrt(mu/SemiMajor)
    Dsma = 2*pi*SemiMajor
    Period = Dsma/Vsma
    v0 = sqrt(mu * (2/Apogee - 1/SemiMajor))

    x = cos(Rotation) * cos(Inc)
    y = sin(Rotation) * cos(Inc)
    z = sin(Inc)

    StartPos = Apogee * np.array([x, y, z])
    StartVel = v0 * np.array([-sin(Rotation), cos(Rotation), 0.0])
    
    return StartPos, StartVel, Period


def thrust_acceleration(state, time):

    if ((time >= start_thrust_time) and
        (time < end_thrust_time)):
        vx = state[3]
        vy = state[4]
        vz = state[5]
        v = sqrt(vx*vx + vy*vy + vz*vz)
        vx_hat = vx / v
        vy_hat = vy / v
        vz_hat = vz / v

        ax = vx_hat * thrust / mass
        ay = vy_hat * thrust / mass
        az = vz_hat * thrust / mass
        a = np.array([ax, ay, az])

    else:
        a = [0.0, 0.0, 0.0]

    return a
    
## Define function to return f(X,t)
def f_func(state, time):
    f=zeros(6)
    f[0] = state[3]
    f[1] = state[4]
    f[2] = state[5]
    r = sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    #[fdragx,fdragy,fdragz] = calc_drag(state,time)

    #j2pert = J2*(3.0/2.0)*(R0/r)**2
    #j2sub  = 5.0*(state[2]/R0)**2

    #print(j2pert*(j2sub-1))

    acc_thrust = thrust_acceleration(state, time)

    f[3] = -(mu*state[0])/r**3 + acc_thrust[0]
    f[4] = -(mu*state[1])/r**3 + acc_thrust[1]
    f[5] = -(mu*state[2])/r**3 + acc_thrust[2]

    return f


def calc_miss_info(EarthPos, EarthVel,
                   CometPos, CometVel,
                   Time):

    ceDist = np.sqrt((EarthPos[0] - CometPos[0])**2 +
                     (EarthPos[1] - CometPos[1])**2 +
                     (EarthPos[2] - CometPos[2])**2)

    minDist = np.min(ceDist)
    i = np.arange(len(Time))
    iPt = i[(ceDist == minDist)]
    iPt = iPt[0]
    tClosest = Time[iPt]

    velDiff = np.sqrt( (CometVel[0][iPt] - EarthVel[0][iPt]) ** 2 +
                       (CometVel[1][iPt] - EarthVel[1][iPt]) ** 2 +
                       (CometVel[2][iPt] - EarthVel[2][iPt]) ** 2)

    xy = [ EarthPos[0][iPt], EarthPos[1][iPt]]

    # deltaTime = (tCollide - StartTime).total_seconds()/86400.0 / 365.25

    return iPt, minDist, velDiff, tClosest, xy


bodies = ['sun', 'earth']
M = {}

M['sun'] = 1.989e30 # kg
M['earth']  = 5.97219e24  # Mass of the Earth in kg

R = {'sun' : 695500.0 * 1000.0,
     'earth': 1.0 * AU}

body = 'sun'

Mp = M[body]
Rp = R[body]

pi = 3.141592
dtor = pi/180.0

mu = G * Mp

# ---------------------------------------------------------------------------------
# Calculate Odegra's position
#
# STK
# collision -> September 5, 2028  18:05 

# This seems to work!

CometPerigee = 0.50 * AU
CometApogee =  5.00 * AU
Rotation = np.pi/2.0 + 3.17941 * np.pi/12.0
Inc = np.pi/7600.0

# inc:
# 3000 -> 4.03 Re
# 3100 -> 1.06
# 3200 -> 1.00
# 3300 -> 0.976
# 3500 -> 0.8854
# 5000 -> 0.6102
# 6200 -> 0.65
# 7500 -> 0.44
# 7600 -> 0.394  <<---
# 7900 -> 0.401
# 8100 -> 0.72
# 8700 -> 0.48
# 10000 -> 0.54
# 15000 -> 1.02
# 20000 -> 1.30

# rot:
# 3.178 -> 7.809
# 3.179 -> 2.47
# 3.17925 -> 1.265
# 3.17931 -> 1.1097
# 3.17934 -> 1.07
# 3.179355 -> 1.0638  <---- 0.391
# 3.179356 ->         <---- 0.388
# 3.179360 ->         <---- 0.370
# 3.179365 ->         <---- 0.343
# 3.1794 ->         <---- 0.268
# 3.17941 ->         <---- 0.264
# 3.17937 -> 1.0645
# 3.17943 -> 1.13
# 3.1795  -> 1.346
# 3.180   -> 3.96
# 3.181 -> 9.7


sma = (CometPerigee + CometApogee) / 2.0
ecc = 1.0 - (CometPerigee / sma)

print('Odegra Inclination : ', Inc)
print('Odegra Rotation : ', Rotation)

StartPos, StartVel, Period = set_initial_conditions(CometApogee,
                                                    CometPerigee,
                                                    Rotation,
                                                    Inc)
TotalTime = Period*1.75
print('Odegra Period : ', Period/86400.0/365.25, ' years')
CometPos, CometVel, Time = propagate(StartPos, StartVel, TotalTime)

thrust = 0.0
CometPosNoThrust, CometVelNoThrust, Time = propagate(StartPos, StartVel, TotalTime)

# ---------------------------------------------------------------------------------
# Calculate Earth's position
#
# https://ssd.jpl.nasa.gov/horizons/app.html#/
# StartTime = datetime.datetime(2021, 09, 22, 00, 00, 0)
# See bottom for more times
StartPos = [1.489117942779481E+08, -2.212111523556938E+06, 2.340664510791190E+04]
StartVel = [7.239067223484917E-02, 2.964591101669518E+01, -1.740084474842973E-03]

StartPos = np.array(StartPos) * 1000.0
StartVel = np.array(StartVel) * 1000.0

thrust = 0.0
EarthPos, EarthVel, Time = propagate(StartPos, StartVel, TotalTime)


iPt, minDist, velDiff, tClosest, xy = calc_miss_info(EarthPos, EarthVel,
                                                     CometPos, CometVel, Time)

print('#---------------------------------------')
print('With Thrust')
print('Minimum distance in Earth Radii : ', minDist*AU/Re)
print("Velocity Diff (km/s) : ", velDiff/1000.0)
print('Time of closest approach : ', tClosest)

iPtn, minDistn, velDiffn, tClosestn, xyn = calc_miss_info(EarthPos, EarthVel,
                                                     CometPosNoThrust, CometVelNoThrust, Time)

print('#---------------------------------------')
print('With NO Thrust')
print('Minimum distance in Earth Radii : ', minDistn*AU/Re)
print("Velocity Diff (km/s) : ", velDiffn/1000.0)
print('Time of closest approach : ', tClosestn)

xThrust = []
yThrust = []
zThrust = []

for i, t in enumerate(Time):
    if ((t >= StartThrust) and
        (t < EndThrust)):
        xThrust.append(CometPos[0][i])
        yThrust.append(CometPos[1][i])
        zThrust.append(CometPos[2][i])
        
# Plot the results
plt.figure(1,figsize=(10,10))

ax = plt.subplot(221,aspect=1.0)
    
ax.plot(CometPos[0], CometPos[1], 'b', label = 'Odegra')
ax.plot(EarthPos[0], EarthPos[1], 'r', label = 'Earth')

if (len(xThrust) > 0):
    ax.plot(xThrust, yThrust, 'y.', label = 'Thrust')

mindist = minDist * 1.05

ax.plot(CometPos[0][iPt], CometPos[1][iPt], 'k.')
ax.plot(EarthPos[0][iPt], EarthPos[1][iPt], 'k.')
ax.set_xlabel('Sun Centered - AU')
ax.set_ylabel('Sun Centered - AU')
ax.legend()
ax.grid()

ax.set_xlim([-4, 1])
ax.set_ylim([-1, 4])

# ------------------------------------------------------------------------------
# Zoomed a bit

ax2 = plt.subplot(222,aspect=1.0)
    
ax2.plot(CometPos[0], CometPos[1], 'b', label = 'Odegra')
ax2.plot(EarthPos[0], EarthPos[1], 'r', label = 'Earth')

if (len(xThrust) > 0):
    ax2.plot(xThrust, yThrust, 'y.', label = 'Thrust')

mindist = minDist * 1.05

ax2.plot(CometPos[0][iPt], CometPos[1][iPt], 'k.')
ax2.plot(EarthPos[0][iPt], EarthPos[1][iPt], 'k.')
ax2.set_xlabel('Sun Centered - AU')
ax2.set_ylabel('Sun Centered - AU')
ax2.legend()
ax2.grid()

ax2.set_xlim([0, 1])
ax2.set_ylim([0, 1])

# ------------------------------------------------------------------------------
# Zoomed a bit

ax3 = plt.subplot(223,aspect=1.0)
    
xCenter = EarthPos[0][iPt]
yCenter = EarthPos[1][iPt]

nPts = int(4 * 3600.0 / dt)
i = range(iPt-nPts, iPt+nPts)
xC = (CometPos[0] - xCenter) * AU / Re
yC = (CometPos[1] - yCenter) * AU / Re
xE = (EarthPos[0] - xCenter) * AU / Re
yE = (EarthPos[1] - yCenter) * AU / Re

xCenter = EarthPos[0][iPtn]
yCenter = EarthPos[1][iPtn]
xCn = (CometPosNoThrust[0] - xCenter) * AU / Re
yCn = (CometPosNoThrust[1] - yCenter) * AU / Re

ax3.plot(xC[i], yC[i], 'b', label = 'Odegra')
ax3.plot(xCn[i], yCn[i], 'c', label = 'Odegra (No Thrust)')
ax3.plot(xE[i], yE[i], 'r', label = 'Earth')

ax3.plot(xC[iPt], yC[iPt], 'k.')
#ax3.plot(xE[iPt], yE[iPt], 'k.')
ax3.plot(xCn[iPtn], yCn[iPtn], 'ko')
#ax3.plot(xE[iPtn], yE[iPtn], 'ko')
ax3.set_xlabel('Earth Centered - Re')
ax3.set_ylabel('Earth Centered - Re')
ax3.legend()
ax3.grid()

d = 66.0 
ax3.set_xlim([-d, d])
ax3.set_ylim([-d, d])


# ------------------------------------------------------------------------------
# Zoomed within GeoSync

ax4 = plt.subplot(224,aspect=1.0)
    
nPts = int(1 * 3600.0 / dt)

ax4.plot(xC[i], yC[i], 'b', label = 'Odegra')
ax4.plot(xCn[i], yCn[i], 'c', label = 'Odegra (No Thrust)')
ax4.plot(xE[i], yE[i], 'r', label = 'Earth')

ax4.plot(xC[iPt], yC[iPt], 'k.', label = 'Closest')
#ax4.plot(xE[iPt], yE[iPt], 'k.')
ax4.plot(xCn[iPtn], yCn[iPtn], 'ko', label = 'Closest (No Thrust)')
#ax4.plot(xE[iPtn], yE[iPtn], 'ko')
ax4.set_xlabel('Earth Centered - Re')
ax4.set_ylabel('Earth Centered - Re')
ax4.legend()
ax4.grid()

d = 6.6 
ax4.set_xlim([-d, d])
ax4.set_ylim([-d, d])

plt.savefig('asteroid_zoomed.png')

# This is Earth's position and velocity at different times:
# 459477.500000000 = A.D. 2021-Sep-20 00:00:00.0000 TDB 
# X = 1.488119671168073E+08 Y =-7.332252086316249E+06 Z = 2.367041463139933E+04
# VX= 1.082603873923091E+00 VY= 2.960962255535691E+01 VZ=-1.294788604189634E-03
# LT= 4.969854731411649E+02 RG= 1.489924965832828E+08 RR=-3.758634602559419E-01
# 459478.500000000 = A.D. 2021-Sep-21 00:00:00.0000 TDB 
# X = 1.488837013297233E+08 Y =-4.772965812936675E+06 Z = 2.354818476424902E+04
# VX= 5.778067271562314E-01 VY= 2.963178097561621E+01 VZ=-1.530969139100335E-03
# LT= 4.968777108071541E+02 RG= 1.489601902482899E+08 RR=-3.719484336050012E-01
# 459479.500000000 = A.D. 2021-Sep-22 00:00:00.0000 TDB 
# X = 1.489117942779481E+08 Y =-2.212111523556938E+06 Z = 2.340664510791190E+04
# VX= 7.239067223484917E-02 VY= 2.964591101669518E+01 VZ=-1.740084474842973E-03
# LT= 4.967710891129521E+02 RG= 1.489282258685089E+08 RR=-3.679643459448752E-01
# 459480.500000000 = A.D. 2021-Sep-23 00:00:00.0000 TDB 
# X = 1.488961930903900E+08 Y = 3.496146286696873E+05 Z = 2.324853455396302E+04
# VX=-4.336269815401860E-01 VY= 2.965195074658891E+01 VZ=-1.913301766196085E-03
# LT= 4.966656144494224E+02 RG= 1.488966053598726E+08 RR=-3.640022240498886E-01
# 459481.500000000 = A.D. 2021-Sep-24 00:00:00.0000 TDB 
# X = 1.488368473153884E+08 Y = 2.911510224069408E+06 Z = 2.307724247805681E+04
# VX=-9.402070597552906E-01 VY= 2.964981853716674E+01 VZ=-2.044378895822874E-03
# LT= 4.965612694231133E+02 RG= 1.488653235079554E+08 RR=-3.601359315935472E-01
# 459482.500000000 = A.D. 2021-Sep-25 00:00:00.0000 TDB 
# X = 1.487337108797298E+08 Y = 5.472865252589719E+06 Z = 2.289659341042372E+04
# VX=-1.447287502652985E+00 VY= 2.963942163688433E+01 VZ=-2.129425861836509E-03
# LT= 4.964580178087652E+02 RG= 1.488343694526975E+08 RR=-3.564225732789025E-01
# 459483.500000000 = A.D. 2021-Sep-26 00:00:00.0000 TDB 
# X = 1.485867439963445E+08 Y = 8.032961465116822E+06 Z = 2.271066108334856E+04
# VX=-1.954785379075334E+00 VY= 2.962066301719807E+01 VZ=-2.166487313230547E-03
# LT= 4.963558090410125E+02 RG= 1.488037280349838E+08 RR=-3.529051321221880E-01
# 459484.500000000 = A.D. 2021-Sep-27 00:00:00.0000 TDB 
# X = 1.483959148141791E+08 Y = 1.059107220716252E+07 Z = 2.252362109164009E+04
# VX=-2.462600600108289E+00 VY= 2.959344560686620E+01 VZ=-2.155089033744062E-03
# LT= 4.962545817482808E+02 RG= 1.487733808560790E+08 RR=-3.496162944948154E-01
# 459485.500000000 = A.D. 2021-Sep-28 00:00:00.0000 TDB 
# X = 1.481612007252686E+08 Y = 1.314646249689538E+07 Z = 2.233963978061453E+04
# VX=-2.970619914244524E+00 VY= 2.955767373490762E+01 VZ=-2.095870387735843E-03
# LT= 4.961542661432853E+02 RG= 1.487433069942817E+08 RR=-3.465823908594555E-01
# 459486.500000000 = A.D. 2021-Sep-29 00:00:00.0000 TDB 
# X = 1.478825893600919E+08 Y = 1.569838911861306E+07 Z = 2.216278792418819E+04
# VX=-3.478719941191281E+00 VY= 2.951325213843759E+01 VZ=-1.990389477988685E-03
# LT= 4.960547853744198E+02 RG= 1.487134834100597E+08 RR=-3.438265265336531E-01
# 459487.500000000 = A.D. 2021-Sep-30 00:00:00.0000 TDB 
# X = 1.475600793991559E+08 Y = 1.824610055001538E+07 Z = 2.199696177314874E+04
# VX=-3.986768244583784E+00 VY= 2.946008325183934E+01 VZ=-1.841153558874709E-03
# LT= 4.959560561851303E+02 RG= 1.486838851437263E+08 RR=-3.413701353220114E-01
# 459488.500000000 = A.D. 2021-Oct-01 00:00:00.0000 TDB 
# X = 1.471936814053875E+08 Y = 2.078883661660334E+07 Z = 2.184579084678181E+04
# VX=-4.494621666732802E+00 VY= 2.939806379980486E+01 VZ=-1.651896402247388E-03
# LT= 4.958579894343059E+02 RG= 1.486544854714486E+08 RR=-3.392322835012217E-01
# 459489.500000000 = A.D. 2021-Oct-02 00:00:00.0000 TDB 
# X = 1.467834189383161E+08 Y = 2.332582787074235E+07 Z = 2.171251150561683E+04
# VX=-5.002121418423893E+00 VY= 2.932708210571035E+01 VZ=-1.428085998011142E-03
# LT= 4.957604912140604E+02 RG= 1.486252562403506E+08 RR=-3.374261338069189E-01
# 459490.500000000 = A.D. 2021-Oct-03 00:00:00.0000 TDB 
# X = 1.463293302347410E+08 Y = 2.585629483611088E+07 Z = 2.159978929987550E+04
# VX=-5.509084949671387E+00 VY= 2.924701802867118E+01 VZ=-1.177581752415335E-03
# LT= 4.956634654418056E+02 RG= 1.485961686455970E+08 RR=-3.359522358878174E-01
# 459491.500000000 = A.D. 2021-Oct-04 00:00:00.0000 TDB 
# X = 1.458314706966284E+08 Y = 2.837944744303337E+07 Z = 2.150948425776511E+04
# VX=-6.015295702230696E+00 VY= 2.915774783935536E+01 VZ=-9.112568405083010E-04
# LT= 4.955668188110942E+02 RG= 1.485671947146186E+08 RR=-3.347890243429079E-01
# 459492.500000000 = A.D. 2021-Oct-05 00:00:00.0000 TDB 
# X = 1.452899162616430E+08 Y = 3.089448517529249E+07 Z = 2.144237469217181E+04
# VX=-6.520493561341524E+00 VY= 2.905915606141821E+01 VZ=-6.432745905566861E-04
# LT= 4.954704686926348E+02 RG= 1.485383096757770E+08 RR=-3.338822911576386E-01
# 459493.500000000 = A.D. 2021-Oct-06 00:00:00.0000 TDB 
# X = 1.447047674055584E+08 Y = 3.340059856162151E+07 Z = 2.139788606277853E+04
# VX=-7.024370668114898E+00 VY= 2.895115455911975E+01 VZ=-3.906314232740016E-04
# LT= 4.953743538066417E+02 RG= 1.485094951578548E+08 RR=-3.331372225154728E-01
# 459494.500000000 = A.D. 2021-Oct-07 00:00:00.0000 TDB 
# X = 1.440761530833600E+08 Y = 3.589697252886048E+07 Z = 2.137390228569135E+04
# VX=-7.526577799183216E+00 VY= 2.883370565277365E+01 VZ=-1.716800649553818E-04
# LT= 4.952784462468092E+02 RG= 1.484807427947518E+08 RR=-3.324177816363418E-01
# 459495.500000000 = A.D. 2021-Oct-08 00:00:00.0000 TDB 
# X = 1.434042335475515E+08 Y = 3.838279167850153E+07 Z = 2.136674768803269E+04
# VX=-8.026744029378269E+00 VY= 2.870684208053826E+01 VZ=-3.712026659385970E-06
# LT= 4.951827621439605E+02 RG= 1.484520574223673E+08 RR=-3.315571169869956E-01
# 459496.500000000 = A.D. 2021-Oct-09 00:00:00.0000 TDB 
# X = 1.426892009763783E+08 Y = 4.085724683620325E+07 Z = 2.137139872745052E+04
# VX=-8.524506717748958E+00 VY= 2.857067520506082E+01 VZ= 9.980515310381577E-05
# LT= 4.950873677326910E+02 RG= 1.484234588973333E+08 RR=-3.303785087244785E-01
# 459497.500000000 = A.D. 2021-Oct-10 00:00:00.0000 TDB 
# X = 1.419312773687655E+08 Y = 4.331954159818677E+07 Z = 2.138191637404449E+04
# VX=-9.019543009090222E+00 VY= 2.842538645929045E+01 VZ= 1.314178216116346E-04
# LT= 4.949923785131520E+02 RG= 1.483949818457242E+08 RR=-3.287209983668796E-01
# 459498.500000000 = A.D. 2021-Oct-11 00:00:00.0000 TDB 
# X = 1.411307100326878E+08 Y = 4.576889747036675E+07 Z = 2.139201455830224E+04
# VX=-9.511592529337820E+00 VY= 2.827120433963588E+01 VZ= 9.056609748014921E-05
# LT= 4.948979513855252E+02 RG= 1.483666733050311E+08 RR=-3.264611344725704E-01
# 459499.500000000 = A.D. 2021-Oct-12 00:00:00.0000 TDB 
# X = 1.402877657320068E+08 Y = 4.820455668655673E+07 Z = 2.139564082444459E+04
# VX=-1.000046533683720E+01 VY= 2.810837533841187E+01 VZ=-1.697805780409567E-05
# LT= 4.948042718979955E+02 RG= 1.483385889012004E+08 RR=-3.235244904008383E-01
# 459500.500000000 = A.D. 2021-Oct-13 00:00:00.0000 TDB 
# X = 1.394027247692434E+08 Y = 5.062578259486955E+07 Z = 2.138746139762551E+04
# VX=-1.048603613664999E+01 VY= 2.793713802465398E+01 VZ=-1.806525594929553E-04
# LT= 4.947115398353291E+02 RG= 1.483107885281982E+08 RR=-3.198859493181354E-01
# 459501.500000000 = A.D. 2021-Oct-14 00:00:00.0000 TDB 
# X = 1.384758759642374E+08 Y = 5.303185815120239E+07 Z = 2.136319533862919E+04
# VX=-1.096823043318758E+01 VY= 2.775770559597343E+01 VZ=-3.869491410792847E-04
# LT= 4.946199559571536E+02 RG= 1.482833323722468E+08 RR=-3.155622332101438E-01
# 459502.500000000 = A.D. 2021-Oct-15 00:00:00.0000 TDB 
# X = 1.375075129505773E+08 Y = 5.542208330084413E+07 Z = 2.131979584314674E+04
# VX=-1.144700874743462E+01 VY= 2.757025733409864E+01 VZ=-6.210478241435169E-04
# LT= 4.945297114232066E+02 RG= 1.482562777415938E+08 RR=-3.106014728630102E-01
 
