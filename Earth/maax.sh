#!/bin/sh

DATE=20100216
ALT=37000

# baseline MAAX orbits:
# satellite A
./kepler.py -starttime=${DATE}:00 -inc=87 -time=336 -sat=maaxA_ltan00_inc80 -perigee=${ALT} -apogee=${ALT} -dt=15 -ltan=0
# satellite B phase 1
./kepler.py -starttime=${DATE}:00 -inc=87 -time=336 -sat=maaxB_ltan00_inc80_delay025 -delay=-0.25 -perigee=${ALT} -apogee=${ALT} -dt=15 -ltan=0
# satellite B phase 2
./kepler.py -starttime=${DATE}:00 -inc=87 -time=336 -sat=maaxB_ltan00_inc80_delay050 -delay=-0.5 -perigee=${ALT} -apogee=${ALT} -dt=15 -ltan=0




