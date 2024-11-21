#!/usr/bin/env python3

"""
Calculation of orbital insertion given initial conditions
"""

## Import needed modules
import matplotlib.pyplot as plt
import numpy as np

# set some constants:

G  = 6.67384e-11 # Gravitational Constants in m3 kg-1 s-2
R0 = 6378137.0   # radius of the Earth in r
M  = 5.97219e24  # Mass of the Earth in kg
mu = G * M
pi = np.pi
dtor = pi / 180.0
rtod = 180.0 / pi

# these things can be altered:

soundSpeed = 343.0 # m/s <-- speed of sound at the surface
scaleHeight = 10.0 * 1000.0 # km -> converted to meters
rho0 = 1.29 # kg/m3
vRotation = 465.1  # m/s <- rotational speed of Earth at equator

orbitAltitude = 400.0 * 1000.0  # km -> converted to meters

dt = 0.01
radius = 0.3 # m
Cd = 0.29 # drag coefficient
mass0 = 250 # kg this is the wet pass
payloadMass = 25 # kg
rocketMass = mass0 * 0.1 # assume rocket body + engine = 10% of wet mass
alt0 = 6000 # meters
v0Vertical = 8.0 * soundSpeed  # assume initial launch is vertical
v0Horizontal = vRotation
Isp = 1700.0 # specific impulse in s
Ve = Isp * 9.8 # exhaust velocity in m/s

area = pi * radius**2
ballCoef = Cd * area

massFlowRate = 5.0 # kg/s <---- this is a critical setting!
maxAcc = 5.0 * 9.8 # in m/s2 <--- basically don't want parts to break

def calc_altitude_reached(currentAltitude, currentSpeed):
    r1 = R0 + currentAltitude
    r2 = mu / (-0.5 * currentSpeed**2 + mu / r1)
    finalAltitude = r2 - R0
    return finalAltitude

def calc_drag(currentAltitude, vHor, vVer, Cb):
    rho = rho0 * np.exp(-currentAltitude / scaleHeight)
    vHorNet = vHor - vRotation
    totalSpeed = np.sqrt(vHorNet**2 + vVer**2) 
    forceMag = 0.5 * rho * Cb * totalSpeed ** 2
    ratio = np.array([-vHorNet / totalSpeed, -vVer / totalSpeed])
    force = forceMag * ratio
    return force

def calc_orbital_speed(currentAltitude):
    r1 = R0 + currentAltitude
    return np.sqrt(mu/r1)

def gravity(currentAltitude):
    r = R0 + currentAltitude
    g = -mu / r**2
    return g

altReached = calc_altitude_reached(alt0, v0Vertical)
orbitalSpeed = calc_orbital_speed(alt0)
print('Potential Altitude Reached on Launch (km) : ', altReached/1000.0)
print('Orbital speed and current horizontal speed (m/s) : ', \
      orbitalSpeed, v0Horizontal)

# First do a super simple calculation ignoring all forces:
deltaVneeded = orbitalSpeed - v0Horizontal
massF = rocketMass + payloadMass
deltaVsimple = Ve * np.log(mass0 / massF)
print('Dry Mass + Payload (kg) : ', massF)
print('DeltaV needed and DeltaV from Rocket Equation (m/s) : ', \
      deltaVneeded, deltaVsimple)

drag = calc_drag(alt0, v0Horizontal, v0Vertical, ballCoef)

massFlowRate = np.abs(drag[1]) / Ve * 1.001

print('mfr : ', massFlowRate)
totalThrust = Ve * massFlowRate
# Want to thrust in the vertical direction just enough to counteract drag
vThrust = -drag[1]
hThrust = np.sqrt(totalThrust**2 - vThrust**2)

hAcc0 = drag[0] + hThrust/mass0
vAcc0 = drag[1]/mass0 + gravity(alt0) + vThrust/mass0

print('initial Drag (N) : ', drag)
print('Initial Accelerations (hor, ver, in m/s) : ', hAcc0, vAcc0)

times = [0.0]
alts = [alt0]
hPoss = [0.0]
vVels = [v0Vertical]
hVels = [v0Horizontal]
vAccs = [vAcc0]
hAccs = [hAcc0]
masses = [mass0]
thrustAcc = [totalThrust / mass0]

hVelGoal = calc_orbital_speed(orbitAltitude)

i = 0
while ((masses[-1] > rocketMass + payloadMass) & \
       (alts[-1] < orbitAltitude)):
    
    # Positions:
    alt = alts[-1]
    hPos = hPoss[-1]
    # Velocities:
    vV = vVels[-1]
    hV = hVels[-1]
    # current mass:
    mass = masses[-1]

    # Drag:
    drag = calc_drag(alt, hV, vV, ballCoef)

    if (hV < hVelGoal):
        accMag = np.sqrt(hAccs[-1]**2 + vAccs[-1]**2)
        if ((accMag > maxAcc) & (Ve * massFlowRate * 0.98 > np.abs(drag[1]))):
            # slowly reduce mass flow rate if accelerating too fast
            massFlowRate = massFlowRate * 0.99
            print('Limiting Mfr : ', massFlowRate)
        if (Ve * massFlowRate < np.abs(drag[1])):
            massFlowRate = np.abs(drag[1]) / Ve * 1.1
    else:
        massFlowRate = 0.0
 
    # thrust:
    totalThrust = Ve * massFlowRate

    if (totalThrust > 0):
        # Want to thrust in the vertical direction just enough to counteract drag
        vThrust = -drag[1]
        # put the rest of the thrust in the horizontal direction
        hThrust = np.sqrt(totalThrust**2 - vThrust**2)
    else:
        vThrust = 0.0
        hThrust = 0.0
        
    # calculate accelerations:
    hAcc = drag[0]/mass + hThrust/mass
    vAcc = drag[1]/mass + gravity(alt) + vThrust/mass

    thrustAcc.append(totalThrust / mass)
    vAccs.append(vAcc)
    hAccs.append(hAcc)
    
    accMag = np.sqrt(hAcc**2 + vAcc**2) / 9.8
    
    # update velocities:
    vVels.append(vV + vAcc * dt)
    hVels.append(hV + hAcc * dt)

    # update positions:

    alt = alt + vV * dt + 0.5 * vAcc * dt * dt
    hPos = hPos + hV * dt + 0.5 * hAcc * dt * dt

    alts.append(alt)
    hPoss.append(hPos)

    mass = mass - massFlowRate * dt
    masses.append(mass)

    times.append(times[-1] + dt)
    
    percentSpeed = int(hV / hVelGoal * 100.0)
    print('time, alt, hVel :',
          int(times[-1]),
          int(alt/1000.0),
          int(hV), percentSpeed, int(mass), accMag)

    i += 1

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

ax.axhline(orbitAltitude/1000.0, label = 'Orbital Alt.')
ax.plot(times, np.array(alts)/1000.0, label = 'Altitude')
ax.plot(times, np.array(hPoss)/1000.0, label = 'Down Range')
ax.set_ylabel('Distance (km)')
ax.legend()

ax2.plot(times, np.array(vVels), label = 'Vertical Velocity')
ax2.plot(times, np.array(hVels), label = 'Horizontal Velocity')
ax2.set_ylabel('Velocity (m/s)')
ax2.legend()

ax3.plot(times, np.array(vAccs)/9.8, label = 'Vertical Accel.')
ax3.plot(times, np.array(hAccs)/9.8, label = 'Horizontal Accel.')
ax3.plot(times, np.array(thrustAcc)/9.8, label = 'Thrust Acc')
ax3.set_ylabel('Acceleration (g)')
ax3.legend()

ax4.plot(times, np.array(masses), label = 'Total Mass (kg)')
ax4.axhline(payloadMass + rocketMass, label = 'Dry Mass', color = 'k')
ax4.set_ylim(0.0, mass0)
ax4.set_ylabel('Mass (kg)')
ax4.set_xlabel('Time (s)')
ax4.legend()


plotfile = 'alt.png'
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()

