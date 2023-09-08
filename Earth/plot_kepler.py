#!/usr/bin/env python3

"""

"""

## Import needed modules
import matplotlib.pyplot as plt
import argparse
import numpy as np
import datetime as dt
import glob

# ----------------------------------------------------------------------
# Figure out arguments
# ----------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Code to make orbits (around Earth)!')
    
    parser.add_argument('filelist', nargs='*', \
                        help = 'list files to use for generating plots')

    parser.add_argument('-sats', default = '', \
                        help = 'satellites to plot (using glob)')
    
    parser.add_argument('-cat', default = False, \
                        action = 'store_true', \
                        help = 'concatenate the files together for one long plot')

    parser.add_argument('-diff', default = False, \
                        action = 'store_true', \
                        help = 'compare different satellites to the first')

    parser.add_argument('-output', default = False, \
                        action = 'store_true', \
                        help = 'write data to a new file')
    
    args = parser.parse_args()

    return args

#------------------------------------------------------------------------
# concatinate two dictionaries 
#------------------------------------------------------------------------

def concat_dict(dict1, dict2):
    dictOut = {}
    for key in dict1.keys():
        if (not isinstance(dict1[key], str)):
            dictOut[key] = np.concatenate((dict1[key], dict2[key]))
        else:
            dictOut[key] = dict1[key]
    return dictOut
        
# ----------------------------------------------------------------------
# read in one csv kepler file
# ----------------------------------------------------------------------

def read_kepler_file(filename):

    fpin = open(filename, 'r')
    sat = fpin.readline().strip()
    header = fpin.readline().strip().replace(' ','')

    keys = header.split(',')
    data = {}
    for key in keys:
        k = key.replace(' ','')
        data[k] = []
    data['times'] = []
    data['sat'] = sat
    line = fpin.readline()
    while (len(line) > 10):
        vals = line.split(',')
        for i, key in enumerate(keys):
            data[key].append(float(vals[i]))
        data['times'].append(dt.datetime(int(data['year'][-1]), \
                                         int(data['mon'][-1]), \
                                         int(data['day'][-1]), \
                                         int(data['hr'][-1]), \
                                         int(data['min'][-1]), \
                                         int(data['sec'][-1])))
        line = fpin.readline()
    fpin.close()

    return data


#------------------------------------------------------------------------
# Read all kepler files:
#------------------------------------------------------------------------

def read_all_kepler_files(files, doCat):

    alldata = []
    for iFile, file in enumerate(files):

        print('Reading file : ', file)
        alldata.append(read_kepler_file(file))

    if (doCat):
        if (len(alldata) >= 2):
            print('Concatenating files...')
            data = concat_dict(alldata[0], alldata[1])
            for iFile in np.arange(2, len(alldata)):
                data = concat_dict(data, alldata[iFile])
        else:
            data = alldata[0]
        alldata = [data]

    return alldata

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def great_circle(lat1, lon1, lat2, lon2):

    l = len(lat1)
    if (len(lat2) < l):
        l = len(lat2)
    
    la1 = lat1[:l] * np.pi/180.0
    la2 = lat2[:l] * np.pi/180.0
    lo1 = lon1[:l] * np.pi/180.0
    lo2 = lon2[:l] * np.pi/180.0

    dist = np.arccos(np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2) +
                     np.sin(la1) * np.sin(la2)) * 180.0 / np.pi
    
    return dist

#------------------------------------------------------------------------
# Find all equatorward crossings from south to north
#------------------------------------------------------------------------

def equatorward_crossing(lat):

    iC2_ = []
    for i in np.arange(1,len(lat)):
        if ((lat[i] >= 0.0) & (lat[i-1] < 0.0)):
            iC2_.append(i)

    return iC2_

#------------------------------------------------------------------------
# Write data to file:
#------------------------------------------------------------------------

def write_lines(fpout, times, lon, lat, alt):
    for i, t in enumerate(times):
        timeS = t.strftime(' %Y, %m, %d, %H, %M, %S, ')
        fpout.write(timeS)
        xs = "{:.3f}".format(lon[i])
        ys = "{:.3f}".format(lat[i])
        zs = "{:.3f}".format(alt[i])
        fpout.write(xs+", "+ys+", "+zs+"\n")
    return

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()

    if (len(args.sats) > 0):
        satList = args.sats.split(',')
        filelists = []
        for sat in satList:
            list = glob.glob(sat+'*.csv')
            filelists.append(sorted(list))
        print(filelists)
    else:
        filelists = [args.filelist]
        
    for iFile, files in enumerate(filelists):

        if (len(files) > 0):
        
            alldata = read_all_kepler_files(files, args.cat)

            if ((args.output) & (args.cat)):
                orbitfile = alldata[0]['sat'] + '.csv'
                print('Writing output file : ', orbitfile)
                fpout = open(orbitfile, 'w')
                fpout.write(sat + "\n")
                fpout.write("year, mon, day, hr, min, sec, ")
                fpout.write("lon (deg), lat (deg), alt (km)\n")
                write_lines(fpout,
                            alldata[0]['times'],
                            alldata[0]['lon(deg)'],
                            alldata[0]['lat(deg)'],
                            alldata[0]['alt(km)'])
                fpout.close()
                
            
            if (args.diff):
                # We want to plot the differences between the satellites
                # make some assumptions here:
                if (iFile == 0):
                    iCrossRef_ = equatorward_crossing(alldata[0]['lat(deg)'])
                    period = (alldata[0]['times'][iCrossRef_[1]] - alldata[0]['times'][iCrossRef_[0]]).total_seconds()
                    print('  --> Period of reference sat (minutes) : ', period/60.0)
                    fig = plt.figure(figsize = (7,10))
                    ax1 = fig.add_axes([0.13,0.03,0.86,0.29])
                    ax2 = fig.add_axes([0.13,0.35,0.86,0.29])
                    ax3 = fig.add_axes([0.13,0.70,0.86,0.29])
                    # should be concatenated...
                    alts = alldata[0]['alt(km)'][iCrossRef_]
                    timesRef = alldata[0]['times'][iCrossRef_]
                    ax1.plot(timesRef, alts, label = alldata[0]['sat'])
                    refSat = alldata[0]['sat']
                    refLats = alldata[0]['lat(deg)'][iCrossRef_]
                    refLons = alldata[0]['lon(deg)'][iCrossRef_]

                else:
                    satLats = alldata[0]['lat(deg)'][iCrossRef_]
                    satLons = alldata[0]['lon(deg)'][iCrossRef_]
                    sat = refSat + ' to ' + alldata[0]['sat']
                    dist = great_circle(refLats, refLons, satLats, satLons)
                    ax2.plot(timesRef, dist * period / 360.0 / 60.0, label = sat)

                    iCross_ = equatorward_crossing(alldata[0]['lat(deg)'])                    
                    alts = alldata[0]['alt(km)'][iCross_]
                    times = alldata[0]['times'][iCross_]
                    ax1.plot(times, alts, label = alldata[0]['sat'])
                    
                    satLats = alldata[0]['lat(deg)'][iCross_]
                    satLons = alldata[0]['lon(deg)'][iCross_]
                    dist = great_circle(refLats, refLons, satLats, satLons)
                    ax3.plot(timesRef, dist, label = sat)
                    
                    if (files == filelists[-1]):

                        ax1.set_ylabel('Altitude (km)')
                        ax2.set_ylabel('Total Seperation (minutes)')
                        ax3.set_ylabel('Long. Seperation (deg)')
                        ax1.legend()
                        ax2.legend()
                        ax3.legend()

                        
                        plotfile = 'test.png'
                        print('Writing plot : ', plotfile)
                        fig.savefig(plotfile)
                        plt.close()

            else:
                for data in alldata:
                    fig = plt.figure(figsize = (10,7))
                    ax = fig.add_axes([0.1,0.1,0.8,0.8])

                    iCross_ = equatorward_crossing(data['lat(deg)'])
                    
                    ax.plot(data['times'][iCross_], data['alt(km)'][iCross_])

                    sTime = data['times'][0].strftime('%Y%m%d_%H%M%S')
                    plotfile = data['sat'] + '_' + sTime + '.png'
                    print('Writing plot : ', plotfile)
                    fig.savefig(plotfile)
                    plt.close()
