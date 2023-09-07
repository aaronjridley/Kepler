#!/bin/sh

#./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcA -perigee=350 -apogee=350
#./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcB -perigee=349.9 -apogee=349.9
./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcC -perigee=350 -apogee=350 -dt=15

#./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcD -perigee=350 -apogee=350
#./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcE -perigee=349.9 -apogee=349.9
#./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcF -perigee=350.1 -apogee=350.1

# let's modify the altitude of gdcC to 340, making it go faster:

./kepler.py -restartin=.restart_gdcC_20030324_000000.csv -time=48 -modify -newalt=340 -dt=15

# then back up to 350:
./kepler.py -restartin=.restart_gdcC_20030326_000000.csv -time=168 -modify -newalt=350 -dt=15

./plot_kepler.py gdcC*csv 
