#!/bin/bash

# Constants
ice=0
lead=1
iron=2

# Jobs

#
## 1 TeV in ice
#
#nohup iceblock -ne 100 -e 1000000 -t $ice -p e- > ice_1TeV.out &

#
## Moliere Radius jobs
#

#nohup iceblock -ne 100 -e 1000 -t $ice  > logs/ice_1GeV.out &
#nohup iceblock -ne 100 -e 1000 -t $iron > logs/iron_1GeV.out &
#nohup iceblock -ne 100 -e 1000 -t $lead > logs/lead_1GeV.out &

#nohup iceblock -ne 100 -e 10000 -t $ice  > logs/ice_10GeV.out &
#nohup iceblock -ne 100 -e 10000 -t $iron > logs/iron_10GeV.out &
#nohup iceblock -ne 100 -e 10000 -t $lead > logs/lead_10GeV.out &

#
## Plotting shower energy vs track length
#

#nohup iceblock -ne 100 -e 200000 -t $ice > logs/ice_200GeV.out &
#nohup iceblock -ne 100 -e 400000 -t $ice > logs/ice_400GeV.out &
#nohup iceblock -ne 100 -e 700000 -t $ice > logs/ice_700GeV.out &


#
## Want to count particles for higher energies
#

#nohup iceblock -ne 100 -e 100000 -t $ice > logs/ice_100GeV.out &
#nohup iceblock -ne 50  -e 500000 -t $ice > logs/ice_500GeV.out &
#nohup iceblock -ne 20  -e 1000000 -t $ice > logs/ice_1000GeV.out &

#nohup iceblock -ne 100 -e 100000 -t $ice -p gamma > logs/ice_100GeV_gam.out &
#nohup iceblock -ne 50  -e 500000 -t $ice -p gamma > logs/ice_500GeV_gam.out &
#nohup iceblock -ne 20  -e 1000000 -t $ice -p gamma > logs/ice_1000GeV_gam.out &

nohup iceblock -ne 20 -e 10000000 -t $ice -c 1 > logs/ice_10000GeV.out &

