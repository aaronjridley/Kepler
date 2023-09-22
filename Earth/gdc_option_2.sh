#!/bin/sh

ALT1=350.0
DT=15.0

./kepler.py -starttime=20300101:0000 -inc=82 -time=336 -sat=gdcA_lt000 -alt=${ALT1} -dt=${DT}
./kepler.py -starttime=20300101:0015 -inc=82 -time=336 -sat=gdcB_lt000 -alt=${ALT1} -dt=${DT}

./kepler.py -starttime=20300101:0015 -inc=82 -time=336 -sat=gdcB_lt005 -ltan=0.5 -alt=${ALT1} -dt=${DT}
