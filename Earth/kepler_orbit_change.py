#!/Users/ridley/bin/Enthought/User/bin/python

"""
kepler.py
Script to solve the orbit problem using the scipy.integrate.odeint().
http://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html
"""

## Import needed modules
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from scipy.integrate import odeint
import re
import sys

# Set some constants:

G  = 6.67384e-11 # Gravitational Constants in m3 kg-1 s-2

# Earth
R0 = 6378137.0   # radius of the Earth in r
M  = 5.97219e24  # Mass of the Earth in kg

# Sun
AU = 1.496e11 # m
M = 1.989e30 # kg

pi = 3.141592
dtor = pi/180.0

mu = G * M

#-----------------------------------------------------------------------------
# Figure out arguments
#-----------------------------------------------------------------------------

def get_args(argv):

    # set initial values

    inclination = 0.0

    year = 2003
    month = 11
    day = 18
    hour = 0
    minute = 0
    second = 0

    dt = 60.0
    time = 24.0*365.0*1.4

    Plotfile = 'kepler'
    Orbitfile = 'kepler.csv'
    DoOutputFile = 0

    help = 0

    r1 = 1.0
    r2 = 1.52
    
    for arg in argv:

        m = re.match(r'-help',arg)
        if m:
            help = 1

        m = re.match(r'-plotfile=(.*)',arg)
        if m:
            Plotfile = m.group(1)

        m = re.match(r'-outputfile=(.*)',arg)
        if m:
            Orbitfile = m.group(1)
            DoOutputFile = 1

        m = re.match(r'-ymd=(.*)',arg)
        if m:
            Ymd = m.group(1)
            m = re.match(r'(\d\d\d\d)(\d\d)(\d\d)',Ymd)
            if m:
                year  = int(m.group(1))
                month = int(m.group(2))
                day   = int(m.group(3))
            else:
                m = re.match(r'(\d\d)(\d\d)(\d\d)',Ymd)
                if m:
                    year  = 2000 + int(m.group(1))
                    month = int(m.group(2))
                    day   = int(m.group(3))
                else:
                    print("Can not understand format of -ymd=YYYYMMDD")
                    help = 1

        m = re.match(r'-hour=(.*)',arg)
        if m:
            hour = int(m.group(1))

        m = re.match(r'-minute=(.*)',arg)
        if m:
            minute = int(m.group(1))

        m = re.match(r'-second=(.*)',arg)
        if m:
            second = int(m.group(1))

        m = re.match(r'-inclination=(.*)',arg)
        if m:
            inclination = float(m.group(1))
            print('Inclination set to : ',inclination, ' deg')

        m = re.match(r'-r1=(.*)',arg)
        if m:
            r1 = float(m.group(1))
            print('r1 set to : ',r1, ' AU')

        m = re.match(r'-r2=(.*)',arg)
        if m:
            r2 = float(m.group(1))
            print('r2 set to : ',r2, ' AU')

        m = re.match(r'-dt=(.*)',arg)
        if m:
            dt = float(m.group(1))
            print('time step set to : ',dt, ' s')

        m = re.match(r'-time=(.*)',arg)
        if m:
            time = float(m.group(1))
            print('Total time set to : ',time, ' hours')

    if (help == 1):
        print('./kepler.py options:')
        print('            -plotfile=plotfile')
        print('            -outputfile=textfile')
        print('            -ymd=YYYYMMDD or YYMMDD')
        print('            -hour=HH (hour to start run)')
        print('            -minute=MM ()')
        print('            -second=SS ()')
        print('            -time=number of hours to run (hours)')
        print('            -dt=time step (seconds)')
        print('            -r1=r1 (inner circular orbit in AU)')
        print('            -r2=r2 (outer circular orbit in AU)')
        print('')

    StartTime = datetime.datetime(year,month,day,hour,minute,second)
    TotalTime = time * 3600.0
    r1 = r1 * AU
    r2 = r2 * AU

    args = {'StartTime':StartTime, 
            'TotalTime':TotalTime,
            'r1':r1,
            'r2':r2,
            'Inclination':inclination,
            'Dt':dt,
            'Plotfile':Plotfile,
            'Orbitfile':Orbitfile,
            'DoOutputFile':DoOutputFile,
            'Quit':help}

    return args

## Define function to return f(X,t)
def f_func(state,time):
    f=zeros(7)
    f[0] = state[3]
    f[1] = state[4]
    f[2] = state[5]
    r = sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    f[3] = -(G*M*state[0])/r**3
    f[4] = -(G*M*state[1])/r**3
    f[5] = -(G*M*state[2])/r**3

    return f


#------------------------------------------------------------------------
# Main code is here:
#------------------------------------------------------------------------

args = get_args(sys.argv)

if (args['Quit'] == 0):

    print('Start Time set to :',args['StartTime'])

    # Set initial position and speed of satellite.
    # Compute and display some values for a circular orbit to help decide 
    # on initial conditions.

    EndTime = args['StartTime']+datetime.timedelta(seconds=args['TotalTime'])

    SemiMajor = (args['r1']+args['r2'])/2.0

    v1 = sqrt(mu/args['r1'])
    v2 = sqrt(mu/args['r2'])
    v = sqrt(mu * (2/args['r1'] - 1/SemiMajor)) * 1.0

    dvList = np.arange(0,v, 1000)

    dVStart = dvList * 0.0
    dVEnd = dvList * 0.0
    Days = dvList * 0.0
    
    iList = 0
    
    for dv in dvList:

        #print (dv)
        v = sqrt(mu * (2/args['r1'] - 1/SemiMajor)) + dv

        r1 = args['r1']
        r2 = args['r2']
    
        #print('Inner Radius    : ',r1)
        #print('Outer Radius    : ',r2)
        #print('SemiMajor Axis  : ',SemiMajor)
        #print('V0              : ',v)

        # Get initial position and speed.

        x0 = args['r1']
        y0 = 0.0
        z0 = 0.0

        vx0 = 0.0
        vy0 = v*np.cos(args['Inclination'] * 3.14159/180.0)
        vz0 = v*np.sin(args['Inclination'] * 3.14159/180.0)

        deltaV1 = np.sqrt(vx0**2 + (vy0-v1)**2 + vz0**2)
        print('Delta V1 : ',deltaV1,deltaV1*3.6/1.6,deltaV1/2)
        dVStart[iList] = deltaV1
        
        # Set initial conditions and define needed array.
        X0 = [ x0, y0, z0, vx0, vy0, vz0]       # set initial state of the system

        t0 = 0.

        # create time array starting at t0, 
        t = arange(t0,args['TotalTime']+args['Dt'],args['Dt'])   

        Time = []
        for ti in t:
            Time.append(args['StartTime'] + datetime.timedelta(seconds=ti))

        # We need UT so we can convert to local time
        UT = []
        for ti in Time:
            UT.append(float(ti.hour)+float(ti.minute)/60.0+float(ti.second)/3600.0)

        if (dv == dvList[0]):
            args['TotalTime'] = args['TotalTime']/2
            
        ## Solve the ODE with odeint
        X = odeint(f_func,X0,t) # returns an 2-dimensional array with 
                                # the first index specifying the time
                                # and the second index specifying the 
                                # component of the state vector

        x = X[:,0]  # Giving a ':' as the index specifies all of the 
                    # elements for that index so
        y = X[:,1]  # x, y, vx, and vy are arrays that specify their 
        z = X[:,2]

        vx = X[:,3] # values at any given time index
        vy = X[:,4]
        vz = X[:,5]

        r = np.sqrt(x**2 + y**2 + z**2)
        d = abs(r-r2)
        #print(min(d))
        ind = np.where(d==min(d))[0][0]
        Days[iList] = ind * args['Dt']/(24.0*3600.0)
        print('Time to get there : ',Days[iList],' days')

        # This is where it is on the circle at r2:
        xHat = x[ind]/r[ind]
        yHat = y[ind]/r[ind]
        zHat = z[ind]/r[ind]

        # The direction of travel is 90 deg CCW for a circular orbit:
        v2x = - v2 * yHat
        v2y =   v2 * xHat
        v2z =   v2 * zHat
    
        #print(v2x, vx[ind])
        #print(v2y, vy[ind])
        #print(v2z, vz[ind])

        deltaV2 = np.sqrt((v2x-vx[ind])**2 + (v2y-vy[ind])**2 + (v2z-vz[ind])**2)
        dVEnd[iList] = deltaV2
        print('delta V2 : ',deltaV2,deltaV2*3.6/1.6,deltaV2/2)

        # Plot the results
        fig = plt.figure()

        # Plot the orbit
        ax1 = fig.add_subplot(1,1,1,aspect=1.0)
        ax1.plot(x/AU,y/AU)
        rmax = r2*1.1/AU
        ax1.set_ylim([-rmax,rmax])
        ax1.set_xlim([-rmax,rmax])
    
        theta = arange(0,361,2)*3.1415/180.0
        x1 = r1*cos(theta)/AU
        y1 = r1*sin(theta)/AU
        ax1.plot(x1,y1)

        x2 = r2*cos(theta)/AU
        y2 = r2*sin(theta)/AU
        ax1.plot(x2,y2)
        plotfile = args['Plotfile']+str(iList)+'.png'
        plt.savefig(plotfile)

        #fig.close()

        iList = iList + 1
    

        
    # ----
    # Plot out variables as a function of initial delta-V

    plt.figure(1)
    ax1 = plt.subplot(1,1,1)
    ax1.plot(dVStart, Days)
    plt.savefig('days.png')

    plt.figure(1)
    ax2 = plt.subplot(1,1,1)
    ax2.plot(dVStart, dVEnd)
    #ax1.plot(dVEnd)
    plt.savefig('test.png')
        
    if (args['DoOutputFile'] == 1):

        fpout = open(args['Orbitfile'],'w')
        fpout.write("Satellite: A\n")
        fpout.write("Time (UTCG), x (km), y (km), z (km), ")
        fpout.write("vx (m/sec), vy (m/sec), vz (m/sec)\n")

        i = 0
        for x in xGeo:
            #timeS = Time.strftime('%Y-%m-%dT%H:%M:%SZ')
            timeS = Time[i].strftime('%d %b %Y %H:%M:%S.000')
            fpout.write(timeS+", ")
            xs = "{:.3f}".format(x[i]/AU)
            ys = "{:.3f}".format(y[i]/AU)
            zs = "{:.3f}".format(z[i]/AU)
            fpout.write(xs+", "+ys+", "+zs+", ")
            xs = "{:.3f}".format(vx[i])
            ys = "{:.3f}".format(vy[i])
            zs = "{:.3f}".format(vz[i])
            fpout.write(xs+", "+ys+", "+zs+"\n")
            i=i+1

        fpout.close()

