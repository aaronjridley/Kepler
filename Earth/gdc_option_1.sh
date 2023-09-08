#!/bin/sh

# I created a spreadsheet to look at orbital periods and precession rates, here:
# https://docs.google.com/spreadsheets/d/1fyjCOK9WrXZvQXkkJqWU14FUMQWBYeYWcWiWIfEbKDU

rm -rf gdc*csv .restart_gdc*

# if we want to create a constellation with 3 sats in a pearls on a string, with spacing of
# say 1/8 and 1/8 orbits (45 deg) fixed, then we need to do the following:
# SC A - stays at 350 km (reference)
# SC B - moves to 350.32 km for 14 days, then back to 350 km - seperates by 1/32 orbit
# SC C - moves to 350.64 km for 14 days, then back to 350 km - seperates by 1/64 orbit
# We need to correct for drag after that...

DT=15.0

ALT1=350.0
ALT2=350.5
ALT3=351.0
ALT1b=350.25
ALT2b=350.75
ALT3b=351.25

#ALT1=350
#ALT2=350.32
#ALT3=350.64
#ALT1b=350.16
#ALT2b=350.48
#ALT3b=350.80

./kepler.py -starttime=20300101:00 -inc=82 -time=336 -sat=gdcA -alt=${ALT1} -dt=${DT}
./kepler.py -starttime=20300101:00 -inc=82.1 -time=336 -sat=gdcB -alt=${ALT2} -dt=${DT}
./kepler.py -starttime=20300101:00 -inc=82.2 -time=336 -sat=gdcC -alt=${ALT3} -dt=${DT}

./kepler.py -starttime=20300101:00 -inc=81.3 -time=336 -sat=gdcD -alt=${ALT1b} -dt=${DT}
./kepler.py -starttime=20300101:00 -inc=81.4 -time=336 -sat=gdcE -alt=${ALT2b} -dt=${DT}
./kepler.py -starttime=20300101:00 -inc=81.5 -time=336 -sat=gdcF -alt=${ALT3b} -dt=${DT}

# switch to a 2-week cycle, correcting for drag to being it back up to ${ALT1}:
./kepler.py -restartin=.restart_gdcA_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300115_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300129_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

# seperate more (1/16 largest):

./kepler.py -restartin=.restart_gdcA_20300212_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300212_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300212_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300212_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300212_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300212_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300226_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300312_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

# seperate more (3/32 largest):

./kepler.py -restartin=.restart_gdcA_20300326_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300326_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300326_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300326_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300326_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300326_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300409_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300423_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

# seperate more (1/8 largest):

./kepler.py -restartin=.restart_gdcA_20300507_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300507_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300507_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300507_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300507_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300507_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300521_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300604_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}

# Now let orbits move apart continuously

./kepler.py -restartin=.restart_gdcA_20300618_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300618_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300618_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300618_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300618_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300618_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300702_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300702_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300702_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300702_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300702_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300702_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcA_20300716_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcB_20300716_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcC_20300716_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./kepler.py -restartin=.restart_gdcD_20300716_000000.csv -time=336 -modify -newalt=${ALT1} -dt=${DT}
./kepler.py -restartin=.restart_gdcE_20300716_000000.csv -time=336 -modify -newalt=${ALT2} -dt=${DT}
./kepler.py -restartin=.restart_gdcF_20300716_000000.csv -time=336 -modify -newalt=${ALT3} -dt=${DT}

./plot_kepler.py -sats=gdcA,gdcB,gdcC,gdcD,gdcE,gdcF -cat -diff

./plot_kepler.py -sats=gdcA -cat -output
./plot_kepler.py -sats=gdcB -cat -output
./plot_kepler.py -sats=gdcC -cat -output
./plot_kepler.py -sats=gdcD -cat -output
./plot_kepler.py -sats=gdcE -cat -output
./plot_kepler.py -sats=gdcF -cat -output




