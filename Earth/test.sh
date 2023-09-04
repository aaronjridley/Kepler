#!/bin/sh

./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcA -perigee=350 -apogee=350
./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcB -perigee=349.9 -apogee=349.9
./kepler.py -starttime=20030317:00 -inc=82 -time=168 -sat=gdcC -perigee=350.1 -apogee=350.1

./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcD -perigee=350 -apogee=350
./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcE -perigee=349.9 -apogee=349.9
./kepler.py -starttime=20030317:00 -inc=82.5 -time=168 -sat=gdcF -perigee=350.1 -apogee=350.1

