#!/usr/bin/env python3

"""
kepler.py
Script to solve the orbit problem using the scipy.integrate.odeint().
http://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html
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

    for files in filelists:

        if (len(files) > 0):
        
            alldata = read_all_kepler_files(files, args.cat)

            for data in alldata:
                fig = plt.figure(figsize = (10,7))
                ax = fig.add_axes([0.1,0.1,0.8,0.8])

                ax.plot(data['times'], data['alt(km)'])

                sTime = data['times'][0].strftime('%Y%m%d_%H%M%S')
                plotfile = data['sat'] + '_' + sTime + '.png'
                print('Writing plot : ', plotfile)
                fig.savefig(plotfile)
                plt.close()
