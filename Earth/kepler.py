#!/opt/local/bin/python

"""
kepler.py
Script to solve the orbit problem using the scipy.integrate.odeint().
http://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html
"""

## Import needed modules
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from pylab import *
from scipy.integrate import odeint
import re
import sys

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

pi = 3.141592
dtor = pi/180.0

mu = G * M

Drag = []

#-----------------------------------------------------------------------------
# Figure out arguments
#-----------------------------------------------------------------------------

def get_args(argv):

    # set initial values

    perigee = 400.0
    apogee = 400.0

    inclination = 82.0
    LTAN = 0.0

    year = 2003
    month = 11
    day = 18
    hour = 0
    minute = 0
    second = 0

    dt = 60.0
    time = 24.0
    nOrbits = -1.0
    
    Plotfile = 'kepler.png'
    Orbitfile = 'kepler.csv'
    DoOutputFile = 0

    help = 0

    for arg in argv:

        m = re.match(r'-help',arg)
        if m:
            help = 1

        m = re.match(r'-h',arg)
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

        m = re.match(r'-LTAN=(.*)',arg)
        if m:
            LTAN = int(m.group(1))

        m = re.match(r'-perigee=(.*)',arg)
        if m:
            perigee = float(m.group(1))
            print('Perigee altitude set to : ',perigee, ' km')

        m = re.match(r'-apogee=(.*)',arg)
        if m:
            apogee = float(m.group(1))
            print('Apogee altitude set to : ',apogee, ' km')

        m = re.match(r'-dt=(.*)',arg)
        if m:
            dt = float(m.group(1))
            print('time step set to : ',dt, ' s')

        m = re.match(r'-time=(.*)',arg)
        if m:
            time = float(m.group(1))
            print('Total time set to : ',time, ' hours')

        m = re.match(r'-orbits=(.*)',arg)
        if m:
            nOrbits = float(m.group(1))
            print('Number of orbits set to : ',nOrbits)
            
    if (help == 1):
        print('./kepler.py options:')
        print('            -plotfile=plotfile')
        print('            -outputfile=textfile')
        print('            -ymd=YYYYMMDD or YYMMDD')
        print('            -hour=HH (hour to start run)')
        print('            -minute=MM ()')
        print('            -second=SS ()')
        print('            -time=number of hours to run (hours)')
        print('            -orbits=number of orbits to run')
        print('            -dt=time step (seconds)')
        print('            -perigee=perigee (lowest part of orbit)')
        print('            -apogee=apogee (highest part of orbit)')
        print('            -inclination=inclination')
        print('')

    StartTime = datetime.datetime(year,month,day,hour,minute,second)
    TotalTime = time * 3600.0
    Perigee = perigee*1000.0 + R0
    Apogee  = apogee*1000.0 + R0

    args = {'StartTime':StartTime, 
            'TotalTime':TotalTime,
            'nOrbits':nOrbits,
            'LTAN':LTAN,
            'Perigee':Perigee,
            'Apogee':Apogee,
            'Inclination':inclination,
            'Dt':dt,
            'Plotfile':Plotfile,
            'Orbitfile':Orbitfile,
            'DoOutputFile':DoOutputFile,
            'Quit':help}

    return args

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

## Define function to return f(X,t)
def f_func(state,time):
    f=zeros(6)
    f[0] = state[3]
    f[1] = state[4]
    f[2] = state[5]
    r = sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    [fdragx,fdragy,fdragz] = calc_drag(state,time)

    #print(r,dr,H,fdragx,fdragy)

    j2pert = J2*(3.0/2.0)*(R0/r)**2
    j2sub  = 5.0*(state[2]/R0)**2

    #print(j2pert*(j2sub-1))

    f[3] = -(G*M*state[0])/r**3 * (1.0-j2pert*(j2sub-1)) + fdragx
    f[4] = -(G*M*state[1])/r**3 * (1.0-j2pert*(j2sub-1)) + fdragy
    f[5] = -(G*M*state[2])/r**3 * (1.0-j2pert*(j2sub-3)) + fdragz

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

    SemiMajor = (args['Perigee']+args['Apogee'])/2.0
    if (args['nOrbits'] > 0.0):
        Vsma = sqrt(mu/SemiMajor)
        Dsma = 2*pi*SemiMajor
        Period = Dsma/Vsma
        print(Period)
        args['TotalTime'] = Period * args['nOrbits']
    
    EndTime = args['StartTime']+datetime.timedelta(seconds=args['TotalTime'])

    LTAN = args['LTAN'] * 15.0
    
    v =sqrt(mu * (2/args['Perigee'] - 1/SemiMajor))

    print('Perigee Radius  : ',args['Perigee'])
    print('Apogee Radius   : ',args['Apogee'])
    print('SemiMajor Axis  : ',SemiMajor)
    print('V0              : ',v)

    # Get initial position and speed.

    x0 = args['Perigee']
    y0 = 0.0
    z0 = 0.0

    vx0 = 0.0
    vy0 = v*np.cos(args['Inclination'] * 3.14159/180.0)
    vz0 = v*np.sin(args['Inclination'] * 3.14159/180.0)

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


    print('Going into the ODE solver... this could take a bit...')
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

    print('Making figures...')
    
    # Plot the results
    plt.figure(1,figsize=(6,9))

    # -------------------------------------------------
    # View in Local Time Frame, North Polar
    # -------------------------------------------------

    ax1 = plt.subplot(3,2,1,aspect=1.0)
    theta = arange(0,361,2)*3.1415/180.0
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

    r = np.sqrt(x**2 + y**2 + z**2)/1000.0
    alt = (r - R0/1000.0)

    lat = np.arcsin(z/1000.0/r)/dtor
    xy = np.sqrt(x**2 + y**2)
    lon = np.arccos(x/xy)/dtor
    i = 0
    for lons in lon:
        if (y[i] < 0.0):
            lon[i] = 360-lons
        lon[i] = (lon[i] + UT[i]*15.0 + LTAN) % 360
        i+=1

    xGeo = xy * np.cos(lon*dtor)
    yGeo = xy * np.sin(lon*dtor)
    zGeo = z

    nP = len(xGeo)
    vxGeo = xGeo*0.0
    vxGeo[1:nP-2] = (xGeo[2:nP-1] - xGeo[0:nP-3])/args['Dt']/2
    vxGeo[0]      = (xGeo[1]      - xGeo[0])     /args['Dt']
    vxGeo[nP-1]   = (xGeo[nP-1]   - xGeo[nP-2])  /args['Dt']

    vyGeo = yGeo*0.0
    vyGeo[1:nP-2] = (yGeo[2:nP-1] - yGeo[0:nP-3])/args['Dt']/2
    vyGeo[0]      = (yGeo[1]      - yGeo[0])     /args['Dt']
    vyGeo[nP-1]   = (yGeo[nP-1]   - yGeo[nP-2])  /args['Dt']

    vzGeo = vz

#    ax3 = plt.subplot(3,2,3,aspect=1.0)
#    map = Basemap( lat_0=90, lon_0=-90)
#    #map = Basemap(projection='ortho', lat_0=90, lon_0=-90)
#    map.drawmapboundary(fill_color='aqua')
#    map.fillcontinents(color='coral',lake_color='aqua')
#    map.drawcoastlines()
#    parallels = np.arange(0.,81,10.)
#    map.drawparallels(parallels)
#    meridians = np.arange(0.,360.,45.)
#    map.drawmeridians(meridians)
#
#    for i in np.arange(1,len(lon)-1):
#        if (lat[i] > 0):
#            xm,ym = map(lon[i],lat[i])
#            if (lat[i]-lat[i-1] > 0):
#                c = 'b.'
#            else:
#                c = 'g.'
#            map.plot(xm, ym, c)
#
#    xm,ym = map(-180.0, 0.0)
#    plt.text(xm-637200.0, ym, 'Longitude', rotation = 90,
#             horizontalalignment='right',
#             verticalalignment='center')
#    
#    ax3 = plt.subplot(3,2,4,aspect=1.0)
#    map = Basemap(projection='ortho', lat_0=-90, lon_0=90)
#    map.drawmapboundary(fill_color='aqua')
#    map.fillcontinents(color='coral',lake_color='aqua')
#    map.drawcoastlines()
#    parallels = np.arange(-80.,0,10.)
#    map.drawparallels(parallels) # ,labels=[False,True,True,False])
#    meridians = np.arange(0.,360.,45.)
#    map.drawmeridians(meridians)
#
#    for i in np.arange(1,len(lon)-1):
#        if (lat[i] < 0):
#            xm,ym = map(lon[i],lat[i])
#            if (lat[i]-lat[i-1] > 0):
#                c = 'b.'
#            else:
#                c = 'g.'
#            map.plot(xm, ym, c)
#
#    xm,ym = map(-180.0, 0.0)
#    plt.text(xm+637200.0, ym, 'Longitude', rotation = -90,
#             horizontalalignment='left',
#             verticalalignment='center')
#    
#    ax3 = plt.subplot(3,1,3,aspect=0.5)
#    map = Basemap(projection='robin', lat_0 = 0, lon_0 = -90)
#    map.drawmapboundary(fill_color='aqua')
#    map.fillcontinents(color='coral',lake_color='aqua')
#    map.drawcoastlines()
#    parallels = np.arange(-80.,81,10.)
#    map.drawparallels(parallels)
#    meridians = np.arange(0.,360.,45.)
#    map.drawmeridians(meridians)
#                    
#    for i in np.arange(1,len(lon)-1):
#        xm,ym = map(lon[i],lat[i])
#        if (lat[i]-lat[i-1] > 0):
#            c = 'b.'
#        else:
#            c = 'g.'
#        map.plot(xm, ym, c)

#    ax4 = subplot(2,2,4)
#
#    #---------------------------------------
#    # if you have basemap, use these lines:
#    map = Basemap(projection='hammer',lon_0=275.0)
#    map.drawmapboundary(fill_color='#99ffff')
#    map.fillcontinents(color='#cc9966',lake_color='#99ffff')
#    xm,ym = map(lon,lat)
#    map.plot(xm, ym, 'b.')

# if not, use this line:
    #plot(Drag,marker='.',linestyle='None')
#    plot(lon,lat,marker='.',linestyle='None')

    #print(Drag)

#eps_plot_min = 1.1*min(eps)
#ax2.axis([min(t),max(t),eps_plot_min,0])
#ax2.xlabel('t')
#ax2.ylabel('E/m')
    plt.savefig(args['Plotfile'])

    if (args['DoOutputFile'] == 1):

        print('Writing output file')
        
        fpout = open(args['Orbitfile'],'w')
        fpout.write("Satellite: A\n")
        fpout.write("year, mon, day, hr, min, sec, ") 
        fpout.write("lon (deg), lat (deg), alt (km), ")
        fpout.write("x (km), y (km), z (km), ")
        fpout.write("vx (km/s), vy (km/s), vz (km/s)\n")

        i = 0
        for x in xGeo:
            timeS = Time[i].strftime(' %Y, %m, %d, %H, %M, %S, ')
            #timeS = Time[i].strftime('%d %b %Y %H:%M:%S.000')
            fpout.write(timeS)
            xs = "{:.3f}".format(lon[i])
            ys = "{:.3f}".format(lat[i])
            zs = "{:.3f}".format(alt[i])
            fpout.write(xs+", "+ys+", "+zs+", ")
            xs = "{:.3f}".format(xGeo[i]/1000.0)
            ys = "{:.3f}".format(yGeo[i]/1000.0)
            zs = "{:.3f}".format(zGeo[i]/1000.0)
            fpout.write(xs+", "+ys+", "+zs+", ")
            xs = "{:.3f}".format(vxGeo[i]/1000.0)
            ys = "{:.3f}".format(vyGeo[i]/1000.0)
            zs = "{:.3f}".format(vzGeo[i]/1000.0)
            fpout.write(xs+", "+ys+", "+zs+"\n")
            i=i+1

        fpout.close()

