#!/usr/bin/env python3

"""
kepler.py
Script to solve the orbit problem using the scipy.integrate.odeint().
http://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html
"""

## Import needed modules
import matplotlib.pyplot as plt
import argparse
#from mpl_toolkits.basemap import Basemap
from pylab import *
from scipy.integrate import odeint
import re
import sys
import os

# Set some constants:

G  = 6.67384e-11 # Gravitational Constants in m3 kg-1 s-2
R0 = 6378137.0   # radius of the Earth in r
M  = 5.97219e24  # Mass of the Earth in kg
J2 = 1082.645e-6
J3 =   -2.546e-6
J4 =   -1.649e-6
J5 =   -0.210e-6
J6 =    0.646e-6
J7 =   -0.333e-6
J8 =   -0.270e-6
J9 =    0.053e-6

rho0 = 1.0e-12 # kg/m3 (400 km altitude)
Z0 = 400000.0 # m
H = 60000.0 # m
D = 0.0 # 0.5 * 0.1 * 3.0 # 1/2 * A * Cd / M 

pi = np.pi
dtor = pi / 180.0
rtod = 180.0 / pi

mu = G * M

Drag = []

# ----------------------------------------------------------------------
# Figure out arguments
# ----------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Code to make orbits (around Earth)!')
    
    parser.add_argument('-sat',  \
                        default = "gdcA", \
                        help = 'satellite name')

    parser.add_argument('-perigee', default = 400.0, \
                        type = float, \
                        help = 'perigee (as altitude) in km')

    parser.add_argument('-apogee', default = 400.0, \
                        type = float, \
                        help = 'apogee (as altitude) in km')

    parser.add_argument('-alt', default = -1.0, \
                        type = float, \
                        help = 'altitude in km')
    
    parser.add_argument('-delay', default = 0.0, \
                        type = float, \
                        help = 'delay (in +/- fraction of a period) ADDED to start time')
    
    parser.add_argument('-inc', metavar = 'inclination', default = 82.0, \
                        type = float, \
                        help = 'inclination in degrees')
    
    parser.add_argument('-ltan', default = 0.0, \
                        type = float, \
                        help = 'local time of ascending node (hours)')

    # time stuff
    parser.add_argument('-dt', default = 60.0, \
                        type = double, \
                        help = 'delta-t')
    parser.add_argument('-time', default = 24.0, \
                        type = double, \
                        help = 'time to simulate in hours')
    parser.add_argument('-norbits', default = 0, \
                        type = int, \
                        help = 'number of orbits to simulate')
    parser.add_argument('-starttime', default = '20031118:000000', \
                        help = 'start time in YYYYMMDD:HHMMSS')

    # files
    parser.add_argument('-plotfile', default = '', \
                        help = 'plot file to output')
    parser.add_argument('-orbitfile', default = '', \
                        help = 'output file')
    parser.add_argument('-restartin', default = '', \
                        help = 'restart output file')
    parser.add_argument('-logfile', default = False, \
                        action = 'store_true', \
                        help = 'write a logfile-style file')

    # modify the orbit
    parser.add_argument('-modify', default = False, \
                        action = 'store_true', \
                        help = 'modify the orbit at the ascending eq crossing')
    parser.add_argument('-newinc', metavar = 'newinclination', \
                        default = -1.0, \
                        type = float, \
                        help = 'inclination in degrees')
    parser.add_argument('-newalt', metavar = 'newaltitude', \
                        default = -1.0, \
                        type = float, \
                        help = 'new altitude (in km)')
        
    args = parser.parse_args()

    return args

#------------------------------------------------------------------------
# Calculate Drag using an ideal (shitty) atmosphere
#------------------------------------------------------------------------

### Calculate drag given time and position
def calc_drag(state, time):
    x = state[0]
    y = state[1]
    z = state[2]
    r = sqrt(x*x+y*y+z*z)

    dr = r-(R0+Z0)
    rho = rho0*np.exp(-dr/H)

    vx = state[3]
    vy = state[4]
    vz = state[5]

    fdragx = -D*rho*vx**2 * np.sign(vx)
    fdragy = -D*rho*vy**2 * np.sign(vy)
    fdragz = -D*rho*vz**2 * np.sign(vz)

    drag = [fdragx, fdragy, fdragz]

    return drag

#------------------------------------------------------------------------
# ODE solver - this calculates the forces
#------------------------------------------------------------------------

## Define function to return f(X,t)
def f_func(state, time):

    f=zeros(6)
    f[0] = state[3]
    f[1] = state[4]
    f[2] = state[5]
    r = sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    [fdragx, fdragy, fdragz] = calc_drag(state, time)

    j2pert = J2*(3.0/2.0)*(R0/r)**2
    j2sub  = 5.0*(state[2]/R0)**2

    #print(j2pert*(j2sub-1))

    f[3] = -(G*M*state[0])/r**3 * (1.0-j2pert*(j2sub-1)) + fdragx
    f[4] = -(G*M*state[1])/r**3 * (1.0-j2pert*(j2sub-1)) + fdragy
    f[5] = -(G*M*state[2])/r**3 * (1.0-j2pert*(j2sub-3)) + fdragz

    return f

#------------------------------------------------------------------------
# Find the next equatorward cross for ascending or descending node
#------------------------------------------------------------------------

def find_eq_crossing(X, isAscending):
    z = X[:,2]
    vz = X[:,5]
    if (not isAscending):
        vz = -vz
    i = 1
    iCross_ = -1
    isFound = False
    while (not isFound):
        if ((z[i] * z[i-1] <= 0.0) and \
            (vz[i] > 0)):
            iCross_ = i
            isFound = True
        i += 1
        if (i >= len(z)):
            break
    return iCross_

#------------------------------------------------------------------------
# find approximate period.  This is not really correct, since it
# assumes a circular orbit. 
#------------------------------------------------------------------------

def calc_period(X0):
    x0 = np.array(X0)
    r0 = np.sqrt(np.sum(x0[0:3] * x0[0:3]))
    v0 = np.sqrt(np.sum(x0[3:6] * x0[3:6]))
    circ = 2 * np.pi * r0
    period = circ / v0
    return period

#------------------------------------------------------------------------
# store x, y, z, ut, and time in a dictionary
#------------------------------------------------------------------------

def store_vars_in_dict(StartTime, time, X, iStop_):

    # These variables will be needed for outputting:
    t = []
    ut = np.zeros(iStop_)
    x = np.zeros(iStop_)
    y = np.zeros(iStop_)
    z = np.zeros(iStop_)
    vx = np.zeros(iStop_)
    vy = np.zeros(iStop_)
    vz = np.zeros(iStop_)
    for i in range(iStop_):
        t.append(StartTime + \
                 datetime.timedelta(seconds = time[i]))
        ut[i] = \
            float(t[-1].hour) + \
            float(t[-1].minute) / 60.0 + \
            float(t[-i].second) / 3600.0
        x[i] = X[i, 0]
        y[i] = X[i, 1]
        z[i] = X[i, 2]
        vx[i] = X[i, 3]
        vy[i] = X[i, 4]
        vz[i] = X[i, 5]

    data = {'time' : t,
            'ut' : ut,
            'x' : x,
            'y' : y,
            'z' : z,
            'vx' : vx,
            'vy' : vy,
            'vz' : vz}

    return data

#------------------------------------------------------------------------
# change the altitude of the orbit
#------------------------------------------------------------------------

def change_altitude(preData, newR):

    print('Changing Altitude')

    X0 = [preData['x'][-1],
          preData['y'][-1],
          preData['z'][-1],
          preData['vx'][-1],
          preData['vy'][-1],
          preData['vz'][-1]]
    x0 = np.array(X0)
    r0 = np.sqrt(np.sum(x0[0:3] * x0[0:3]))
    v0 = np.sqrt(np.sum(x0[3:6] * x0[3:6]))

    print('  --> Current Alt : ', (r0 - R0)/1000.0, ' km')
    print('  --> Target Alt : ', (newR - R0)/1000.0, ' km')
    
    sma = (newR + r0)/2.0
    vChange = sqrt(mu * (2/r0 - 1/sma))
    
    x0[3:6] = x0[3:6] * vChange / v0

    print('  --> Current Velocity : ', v0)
    print('  --> New Velocity : ', vChange)
    
    period = calc_period(x0)
    print('  --> Approximate Orbital Period : ', period/3600.0, ' hours')
    
    # create time array to go around a couple of times:
    dt = (preData['time'][1] - preData['time'][0]).total_seconds()
    Tstep2 = arange(0.0, 2 * period, dt)
    
    # calculate orbit:
    Xstep2 = odeint(f_func, x0, Tstep2)

    # next - find the equatorward crossing of the descending node:
    iCross_ = find_eq_crossing(Xstep2, False)
    print('  --> descending equatorward crossing found : ', Tstep2[iCross_])
    if (iCross_ > -1):
        preData2 = store_vars_in_dict(preData['time'][-1], Tstep2, Xstep2, iCross_)
        preDataFinal = concat_dict(preData, preData2)
        x0 = np.array(Xstep2[iCross_])
        r0 = np.sqrt(np.sum(x0[0:3] * x0[0:3]))
        v0 = np.sqrt(np.sum(x0[3:6] * x0[3:6]))
        print('  --> Current Alt (descending) : ', (r0 - R0)/1000.0, ' km')
        print('  --> Target Alt (descending)  : ', (newR - R0)/1000.0, ' km')
        # Circularize Orbit
        sma = (newR + r0)/2.0
        vChange = sqrt(mu * (2/r0 - 1/sma))
        x0[3:6] = x0[3:6] * vChange / v0
    else:
        print('Could not find equatorial crossing.... exiting')
        exit()
    
    return preDataFinal, x0
    
#------------------------------------------------------------------------
# Change the inclination of the orbit:
#------------------------------------------------------------------------

def change_inclination(newInclination, Vold):

    Vmag = sqrt(np.sum(Vold * Vold))
    oldInclination = np.arcsin(Vold[2] / Vmag)
    print('Old Inclination : ', oldInclination * rtod)

    Vxy = np.sqrt(Vold[0]**2 + Vold[1]**2)
    phi = np.arccos(Vold[0]/Vxy)
    if (Vold[1] < 0.0):
        phi = 2 * pi - phi

    # Now derive new velocity:
    Vz = Vmag * np.sin(newInclination * dtor)
    Vx = Vmag * np.cos(newInclination * dtor) * np.cos(phi)
    Vy = Vmag * np.cos(newInclination * dtor) * np.sin(phi)
    
    Vnew = [Vx, Vy, Vz]
    return Vnew

#------------------------------------------------------------------------
# Derive X, Y, Z, Vx, Vy, Vz from geographic coordinates (lat, lon, alt)
#------------------------------------------------------------------------
    
def derive_geo(x, y, z, vx, vy, vz, UT, ltan, dt):
    
    xy = np.sqrt(x**2 + y**2)
    lon, lat, alt = convert_xyz_to_lla(x, y, z, UT)

    xGeo = xy * np.cos(lon*dtor)
    yGeo = xy * np.sin(lon*dtor)
    zGeo = z

    nP = len(xGeo)
    vxGeo = xGeo*0.0
    vxGeo[1:nP-2] = (xGeo[2:nP-1] - xGeo[0:nP-3]) / (dt*2)
    vxGeo[0]      = (xGeo[1]      - xGeo[0])      / dt
    vxGeo[nP-1]   = (xGeo[nP-1]   - xGeo[nP-2])   / dt

    vyGeo = yGeo*0.0
    vyGeo[1:nP-2] = (yGeo[2:nP-1] - yGeo[0:nP-3]) / (dt*2)
    vyGeo[0]      = (yGeo[1]      - yGeo[0])      / dt
    vyGeo[nP-1]   = (yGeo[nP-1]   - yGeo[nP-2])   / dt

    vzGeo = vz
    
    return xGeo, yGeo, zGeo, vxGeo, vyGeo, vzGeo

#------------------------------------------------------------------------
# Derive Lon, Lat, Alt from x, y, z (local time) coordinates
#------------------------------------------------------------------------

def convert_xyz_to_lla(x, y, z, UT):

    r = np.sqrt(x**2 + y**2 + z**2)/1000.0
    alt = (r - R0/1000.0)

    lat = np.arcsin(z/1000.0/r)/dtor
    xy = np.sqrt(x**2 + y**2)
    # local time in hours:
    localtimes = np.arccos(x/xy)/dtor/15
    localtimes[y < 0] = 360.0 - localtimes[y < 0]
    lon = ((localtimes - UT) * 15.0 + 360.0) % 360

    return lon, lat, alt
    
#------------------------------------------------------------------------
# Write data to file:
#------------------------------------------------------------------------

def write_lines(fpout, times, lon, lat, alt, useSpaces = False):
    deli = ' '
    if (not useSpaces):
        deli = ',' + deli
    for i, t in enumerate(times):
        timeS = t.strftime('%Y' + deli +
                           '%m' + deli +
                           '%d' + deli +
                           '%H' + deli +
                           '%M' + deli +
                           '%S' + deli)
        fpout.write(timeS)
        xs = "{:.3f}".format(lon[i])
        ys = "{:.3f}".format(lat[i])
        zs = "{:.3f}".format(alt[i])
        fpout.write(xs + deli + ys + deli + zs + "\n")
    return

#------------------------------------------------------------------------
# Write out a data file
#------------------------------------------------------------------------

def open_data_file(orbitfile, sat, isLog = False):

    print('Writing output file : ', orbitfile)
        
    fpout = open(orbitfile, 'w')

    if (isLog):
        fpout.write("\n")
        fpout.write("This orbit file was created with kepler.py\n")
        fpout.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))
        fpout.write("\n")
        fpout.write("#SATELLITE\n")
        fpout.write(sat + "\n")
        fpout.write("\n")
        fpout.write("#DIRECTORY\n")
        fpout.write(os.getcwd() + "\n")
        fpout.write("\n")
        fpout.write("#VARIABLES\n")
        fpout.write("Year\n")
        fpout.write("Month\n")
        fpout.write("Day\n")
        fpout.write("Hour\n")
        fpout.write("Minute\n")
        fpout.write("Second\n")
        fpout.write("Longitude (deg)\n")
        fpout.write("Latitude (deg)\n")
        fpout.write("Altitude (km)\n")
        fpout.write("\n")
        fpout.write("#START\n")    
    else:
        fpout.write(sat + "\n")
        fpout.write("year, mon, day, hr, min, sec, ")
        fpout.write("lon (deg), lat (deg), alt (km)\n")
    
    return fpout

#------------------------------------------------------------------------
# concatinate two dictionaries 
#------------------------------------------------------------------------

def concat_dict(dict1, dict2):
    dictOut = {}
    for key in dict1.keys():
        dictOut[key] = np.concatenate((dict1[key], dict2[key]))
    return dictOut
        
#------------------------------------------------------------------------
# Main code is here:
#------------------------------------------------------------------------

args = get_args()

if (len(args.restartin) > 4):
    print('Restarting from file : ', args.restartin)
    fpin = open(args.restartin, 'r')
    sat = fpin.readline().strip('\n')
    line = fpin.readline()
    fpin.close()
    s = line.split(',')
    year = int(s[0])
    month = int(s[1])
    day = int(s[2])
    hour = int(s[3])
    minute = int(s[4])
    second = int(s[5])

    # can't use nOrbits in restart, since you don't really know the
    # period of the orbit....
    TotalTime = args.time * 3600.0

    x0 = double(s[6])
    y0 = double(s[7])
    z0 = double(s[8])

    vx0 = double(s[9])
    vy0 = double(s[10])
    vz0 = double(s[11])

    StartTime = datetime.datetime(year,month,day,hour,minute,second)
    EndTime = StartTime + datetime.timedelta(seconds = TotalTime)
    
else:

    sat = args.sat
    year = int(args.starttime[0:4])
    month = int(args.starttime[4:6])
    day = int(args.starttime[6:8])
    hour = 0
    minute = 0
    second = 0
    if (len(args.starttime) > 9):
        hour = int(args.starttime[9:11])
    if (len(args.starttime) > 11):
        minute = int(args.starttime[11:13])
    if (len(args.starttime) > 13):
        second = int(args.starttime[13:15])
    
    # Set initial position and speed of satellite.
    # Compute and display some values for a circular orbit to help decide 
    # on initial conditions.

    if (args.alt > 0.0):
        Perigee = args.alt * 1000 + R0
        Apogee = args.alt * 1000 + R0
    else:
        Perigee = args.perigee * 1000 + R0
        Apogee = args.apogee * 1000 + R0

    nOrbits = args.norbits

    SemiMajor = (Perigee + Apogee)/2.0
    if (nOrbits > 0.0):
        Vsma = sqrt(mu/SemiMajor)
        Dsma = 2*pi*SemiMajor
        Period = Dsma/Vsma
        print('Orbital Period : ', Period)
        TotalTime = Period * nOrbits
    else:
        TotalTime = args.time * 3600.0
    
    v = sqrt(mu * (2/Perigee - 1/SemiMajor))

    print('Perigee Radius  : ', Perigee)
    print('Apogee Radius   : ', Apogee)
    print('SemiMajor Axis  : ', SemiMajor)
    print('V0              : ', v)

    # Get initial position and speed.

    # ltan is in hours -> convert to radians:
    #LTAN = args.ltan * 15.0 * dtor

    #x0 = Perigee * np.cos(LTAN)
    #y0 = Perigee * np.sin(LTAN)
    x0 = Perigee 
    y0 = 0.0
    z0 = 0.0
    print('Position Components : ', x0, y0, z0)
    
    #vx0 = -v*np.cos(args.inc * dtor) * np.sin(LTAN)
    #vy0 = v*np.cos(args.inc * dtor) * np.cos(LTAN)
    vx0 = 0.0
    vy0 = v*np.cos(args.inc * dtor)
    vz0 = v*np.sin(args.inc * dtor)
    print('Velocity Components : ', vx0, vy0, vz0)

StartTime = datetime.datetime(year,month,day,hour,minute,second)
print('Start Time set to :', StartTime)
if (np.abs(args.delay) > 0.0):
    print('  --> delay entered!')
    X0 = [ x0, y0, z0, vx0, vy0, vz0]
    period = calc_period(X0)
    print('  --> Approximate Orbital Period : ', period/3600.0, ' hours')
    delay = args.delay * period
    print('  --> delay set to: ', delay, ' seconds')
    StartTime = StartTime + datetime.timedelta(seconds = delay)
    print('  --> new start time : ', StartTime)
    
EndTime = StartTime + datetime.timedelta(seconds = TotalTime)


if (len(args.orbitfile) < 4):
    orbitfile = sat.replace(" ", "") + "_" + \
        StartTime.strftime("%Y%m%d_%H%M%S")
    if (args.logfile):
        orbitfile += '.txt'
    else:
        orbitfile += ".csv"
else:
    orbitfile = args.orbitfile

if (args.modify):
    # Set initial conditions and define needed array.
    X0 = [ x0, y0, z0, vx0, vy0, vz0]
    
    period = calc_period(X0)
    print('Approximate Orbital Period : ', period/3600.0, ' hours')

    # create time array to go around a couple of times:
    Tpre = arange(0.0, 2 * period, args.dt)
    
    # calculate orbit:
    Xpre = odeint(f_func, X0, Tpre) 

    # find first equatorward cross moving northward:
    iCross_ = find_eq_crossing(Xpre, True)
    
    if (iCross_ > -1):
        print('Ascending node, equator crossing found!')
        preData = store_vars_in_dict(StartTime, Tpre, Xpre, iCross_)

        NewStartTime = StartTime + datetime.timedelta(seconds = Tpre[iCross_])
        TotalTime = (EndTime - NewStartTime).total_seconds()
        print('New Start Time : ', NewStartTime)
        StartTime = NewStartTime
        x0, y0, z0, vx0, vy0, vz0 = Xpre[iCross_, :]
        
    else:
        print('Ascending node, equator crossing not found!')
        print('Exiting!')
        exit()

    # Modify the inclination of the orbit:
    if (args.newinc > -1):
        Vold = np.array([vx0, vy0, vz0])
        Vnew = change_inclination(args.newinc, Vold)
        print('Vold : ', Vold)
        print('Vnew : ', Vnew)
        vx0 = Vnew[0]
        vy0 = Vnew[1]
        vz0 = Vnew[2]

    # Modify the altitude of the orbit:
    if (args.newalt > -1):
        X0 = Xpre[iCross_, :]
        newR = args.newalt * 1000 + R0
        preData, x0 = change_altitude(preData, newR)
        x0, y0, z0, vx0, vy0, vz0 = x0
        StartTime = preData['time'][-1] + datetime.timedelta(seconds = args.dt)
        TotalTime = (EndTime - StartTime).total_seconds()
        print('  --> Real New Start Time : ', StartTime)
        
# Set initial conditions and define needed array.
X0 = [ x0, y0, z0, vx0, vy0, vz0]       # set initial state of the system

t0 = 0.

# create time array starting at t0, 
t = arange(t0, TotalTime + args.dt, args.dt)   

Time = []
for ti in t:
    Time.append(StartTime + datetime.timedelta(seconds = ti))

# We need UT so we can convert to local time
UT = []
for ti in Time:
    UT.append(float(ti.hour)+float(ti.minute)/60.0+float(ti.second)/3600.0)


print('Going into the ODE solver...')
## Solve the ODE with odeint -
# returns an 2-dimensional array with the first index specifying the
# time and the second index specifying the component of the state
# vector

X, outs = odeint(f_func, X0, t, full_output = 1) 

# rotate for local time:
if (args.ltan != 0):
    print('Rotating solution to local time : ', args.ltan) 
    LTAN = args.ltan * 15.0 * dtor
    xn = X[:,0] * np.cos(LTAN) - X[:,1] * np.sin(LTAN)
    yn = X[:,0] * np.sin(LTAN) + X[:,1] * np.cos(LTAN)
    vxn = X[:,3] * np.cos(LTAN) - X[:,4] * np.sin(LTAN)
    vyn = X[:,3] * np.sin(LTAN) + X[:,4] * np.cos(LTAN)

    X[:,0] = xn
    X[:,1] = yn
    X[:,3] = vxn
    X[:,4] = vyn
    
# positions
x = X[:,0] 
y = X[:,1]  
z = X[:,2]

# velocities
vx = X[:,3] 
vy = X[:,4]
vz = X[:,5]

if (len(args.plotfile) > 4):

    print('Making figures...')

    # Plot the results
    plt.figure(1,figsize=(6,9))

    # -------------------------------------------------
    # View in Local Time Frame, North Polar
    # -------------------------------------------------

    ax1 = plt.subplot(3, 2, 1, aspect = 1.0)
    theta = arange(0,361,2) * dtor
    for colat in arange(10,91,10):
        rColat = sin(colat*dtor)
        xEarth = rColat*cos(theta)
        yEarth = rColat*sin(theta)

        ax1.plot(xEarth,yEarth,color='#808080',linestyle=':')

    rColat = 1.0
    xEarth = rColat*cos(theta)
    yEarth = rColat*sin(theta)
    ax1.plot(xEarth,yEarth,color='#000000')
    xEarth = [1.0,-1.0]
    yEarth = [0.0,0.0]
    ax1.plot(xEarth,yEarth,color='#808080',linestyle=':')
    xEarth = [0.0,0.0]
    yEarth = [1.0,-1.0]
    ax1.plot(xEarth,yEarth,color='#808080',linestyle=':')

    # Plot the orbit
    #ax1.scatter(x[z>0]/R0,y[z>0]/R0)

    for i in np.arange(1,len(x)-1):
        if (z[i] > 0):
            if (z[i]-z[i-1] > 0):
                c = 'b.'
            else:
                c = 'g.'
            ax1.plot(x[i]/R0, y[i]/R0, c)

    rmax = np.max(np.absolute([x/R0,y/R0]))
    plt.ylim([-rmax,rmax])
    plt.xlim([-rmax,rmax])
    plt.axis('off')
    
    plt.text( 1.15, 0.0, '06', rotation = -90,
              horizontalalignment='center',
              verticalalignment='center')
    plt.text(-1.15, 0.0, '18', rotation = 90,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(0.0, 1.15, '12', rotation = 0,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(0.0, -1.15, '00', rotation = 180,
             horizontalalignment='center',
             verticalalignment='center')

    plt.text(0.0, 1.30, 'North', rotation = 0,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(-1.30, 0.0, 'Local Time', rotation = 90,
             horizontalalignment='center',
             verticalalignment='center')

    # -------------------------------------------------
    # View in Local Time Frame, South Polar
    # -------------------------------------------------

    ax2 = plt.subplot(3,2,2,aspect=1.0)
    theta = arange(0,361,2)*3.1415/180.0
    for colat in arange(10,91,10):
        rColat = sin(colat*dtor)
        xEarth = rColat*cos(theta)
        yEarth = rColat*sin(theta)

    ax2.plot(xEarth,yEarth,color='#808080',linestyle=':')

    rColat = 1.0
    xEarth = rColat*cos(theta)
    yEarth = rColat*sin(theta)
    ax2.plot(xEarth,yEarth,color='#000000')
    xEarth = [1.0,-1.0]
    yEarth = [0.0,0.0]
    ax2.plot(xEarth,yEarth,color='#808080',linestyle=':')
    xEarth = [0.0,0.0]
    yEarth = [1.0,-1.0]
    ax2.plot(xEarth,yEarth,color='#808080',linestyle=':')

    # Plot the orbit
    # ax2.scatter(x[z<0]/R0,y[z<0]/R0)

    for i in np.arange(1,len(x)-1):
        if (z[i] < 0):
            if (z[i]-z[i-1] > 0):
                c = 'b.'
            else:
                c = 'g.'
            ax2.plot(x[i]/R0, y[i]/R0, c)

    rmax = np.max(np.absolute([x/R0,y/R0]))
    plt.ylim([-rmax,rmax])
    plt.xlim([-rmax,rmax])
    plt.axis('off')
    
    plt.text( 1.15, 0.0, '06', rotation = -90,
              horizontalalignment='center',
              verticalalignment='center')
    plt.text(-1.15, 0.0, '18', rotation = 90,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(0.0, 1.15, '12', rotation = 0,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(0.0, -1.15, '00', rotation = 180,
             horizontalalignment='center',
             verticalalignment='center')

    plt.text(0.0, 1.30, 'South', rotation = 0,
             horizontalalignment='center',
             verticalalignment='center')
    plt.text(1.30, 0.0, 'Local Time', rotation = -90,
             horizontalalignment='center',
             verticalalignment='center')

    plt.savefig(args.plotfile)


# This opens the output file and writes the header:
fpout = open_data_file(orbitfile, args.sat, isLog = args.logfile)
    
# If we modified the orbit, then there is a change at the
# equatorward cross, so write out before this time:
if (args.modify):
    preLon, preLat, preAlt = convert_xyz_to_lla(preData['x'], \
                                                preData['y'], \
                                                preData['z'], \
                                                preData['ut'])
    write_lines(fpout, preData['time'], preLon, preLat, preAlt, args.logfile)

lon, lat, alt = convert_xyz_to_lla(x, y, z, UT)
write_lines(fpout, Time, lon, lat, alt, args.logfile)

fpout.close()

# -------------------------------------------------------
# Write last point to a restart file
    
RestartOut = \
    ".restart_" + \
    sat.replace(" ", "") + "_" + \
    Time[-1].strftime("%Y%m%d_%H%M%S") + ".csv"
rout = open(RestartOut, 'w')
i = -1
timeS = Time[i].strftime(' %Y, %m, %d, %H, %M, %S')
xs = ''
for j in range(6):
    xs = xs + ", %f" % X[i,j] 
rout.write(sat)
rout.write("\n")
rout.write(timeS)
rout.write(xs)
rout.close
