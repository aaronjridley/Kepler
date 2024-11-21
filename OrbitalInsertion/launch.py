#!/usr/bin/env python3

"""
Dave Jarvis
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

dt = 0.1
radius = 0.25 # m
Cd = 0.29 # drag coefficient
mass0 = 250 # kg this is the wet pass
payloadMass = 5 # kg
rocketMass = mass0 * 0.1 # assume rocket body + engine = 10% of wet mass
alt0 = 6000 # meters
v0Vertical = 8.0 * soundSpeed  # assume initial launch is vertical
v0Horizontal = vRotation
Isp = 1700.0 # specific impulse in s
Ve = Isp * 9.8 # exhaust velocity in m/s

area = pi * radius**2
ballCoef = Cd * area

massFlowRate = 10.0 # kg/s <---- this is the critical setting!

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
vVel = [v0Vertical]
hVel = [v0Horizontal]
vAcc = []
hAcc = []
mass = [mass0]



