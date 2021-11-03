#!/opt/local/bin/python

## Import needed modules
import matplotlib.pyplot as plt
from pylab import *
from scipy.integrate import odeint
import re
import sys

# Set some constants:

G  = 6.67384e-11 # Gravitational Constants in m3 kg-1 s-2

# Sun
AU = 1.496e11 # m
M = 1.989e30 # kg

pi = 3.141592
dtor = pi/180.0

mu = G * M


rEarth = 1.0 * AU
rMars = 1.52 * AU

v1 = sqrt(mu/rEarth)
v2 = sqrt(mu/rMars)

print(v1, v2)

x = [rEarth]
y = [0.0]

vx = [0.0]
vy = [v1]

r = [np.sqrt(x[0]*x[0] + y[0]*y[0])]

ax = [-mu*x[0] / r[0]**3]
ay = [-mu*y[0] / r[0]**3]

dt = 60.0

t = 0.0

aThrust = 0.005
tThrust = 1000 * 24.0*3600.0

rHalf = (rEarth + rMars)*0.51

DoThrustOut = 1
DoThrustIn = 0

Mass0 = 1000.0
Ve = 30000.0

FuelMass = 3000.0

# athrust = 0.005
#   Days = 83
#   Ve = 30,000
#   Fuel = 2,800
#   Ve = 20,000
#   Fuel = 6,200
#   Ve = 10,000
#   Fuel = 48,500

# athrust = 0.01
#   Days = 54
#   Ve = 30,000
#   Fuel = 5,750
#   Ve = 20,000
#   Fuel = 16,500
#
# athrust = 0.02
#   Days = 42
#   Ve = 30,000
#   Fuel = 13,000


i = 1
timeV = 0
timeR = 0
while (t < 100*24.0*3600.0):

    x.append(x[i-1] + dt * vx[i-1] + 0.5 * dt * dt * ax[i-1])
    y.append(y[i-1] + dt * vy[i-1] + 0.5 * dt * dt * ay[i-1])

    vx.append(vx[i-1] + dt * ax[i-1])
    vy.append(vy[i-1] + dt * ay[i-1])
    
    r.append(np.sqrt(x[i]*x[i] + y[i]*y[i]))

    ax.append(-mu*x[i] / r[i]**3)
    ay.append(-mu*y[i] / r[i]**3)

    v = np.sqrt(vx[i]*vx[i] + vy[i]*vy[i])
    vxHat = vx[i]/v
    vyHat = vy[i]/v

    xHat = x[i]/r[i]
    yHat = y[i]/r[i]

    if (r[i] > rHalf):
        # This makes it so we only turn around once!
        if (DoThrustOut):
            DoThrustIn = 1
        DoThrustOut = 0

    if (r[i] < r[i-1]):
        DoThrustIn = 0
        
    if (DoThrustOut):
#        ax[i] = ax[i] + aThrust * vxHat
#        ay[i] = ay[i] + aThrust * vyHat
        ax[i] = ax[i] + aThrust * xHat
        ay[i] = ay[i] + aThrust * yHat
        Force = aThrust * (Mass0 + FuelMass)
        MassFlowRate = Force/Ve
        FuelMass = FuelMass - MassFlowRate * dt
        
    if (DoThrustIn):

        vxmHat = -yHat
        vymHat = xHat

        vxGoal = v2 * vxmHat
        vyGoal = v2 * vymHat

        dvx = vxGoal - vx[-1]
        dvy = vyGoal - vy[-1]
        dv = np.sqrt(dvx*dvx + dvy*dvy)
        dvxHat = dvx / dv
        dvyHat = dvy / dv
               
        #print(vxGoal, vx[-1], vyGoal, vy[-1], dv)

        ax[i] = ax[i] + aThrust * dvxHat #xHat
        ay[i] = ay[i] + aThrust * dvyHat #yHat
        Force = aThrust * (Mass0 + FuelMass)
        MassFlowRate = Force/Ve
        FuelMass = FuelMass - MassFlowRate * dt
        if (dv < 10.0):
            DoThrustIn = 0.0
            timeV = t/86400.0

    if (FuelMass < 0.0):
        DoThrustIn = 0.0
        DoThrustOut = 0.0
        
    if ((r[i]/rMars > 0.98) & (timeR == 0)):
        timeR = t/86400.0
        
    t = t + dt
    i = i + 1

x = np.array(x)/AU
y = np.array(y)/AU
r = np.array(r)/AU

print("TimeR : ",timeR)
print("TimeV : ",timeV)
if (timeV == 0):
    print("Delta V",dv)
print("Mass of Fuel Left : ", FuelMass)
print(MassFlowRate)

plt.figure(1)
ax1 = plt.subplot(1,1,1)
ax1.plot(x, y)

rmax = rMars*1.1/AU # np.max(r)
ax1.set_xlim([-rmax,rmax])
ax1.set_ylim([-rmax,rmax])
ax1.set_aspect(1.0)

theta = np.arange(361) * dtor
xEarth = rEarth * cos(theta)/AU
yEarth = rEarth * sin(theta)/AU

xMars = rMars * cos(theta)/AU
yMars = rMars * sin(theta)/AU

ax1.plot(xEarth, yEarth, linestyle='dotted', color = 'blue')
ax1.plot(xMars, yMars, linestyle='dotted', color = 'red')

plt.savefig('test.png')
