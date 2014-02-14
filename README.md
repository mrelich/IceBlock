IceBlock
========

Simplified example shooting beam of particles into giant piece of ice.

In order to run, you will need to update setupfiles/setup such that
the environment variables are set specifically for your machine. This
means that you need to source the geant4make.sh script in your 
geant4-build directory.  You also need to make sure that the data
files containing cross-sectional info are pointing in the right 
area.  My $G4INSTALL points to geant4, not geant4-build, so I 
additionally had to link the data/ directory in geant4-build to 
my geant4 directory like so:

cd /path/to/untarred-geant4/
ln -s ../geant4-build/data data

To execute: 
iceblock -ne <number_of_events> -e <beam_energy> -t <enum_for_material_type>

<number_of_events> = number of events specified as an int
<beam_energy> = energy of the particle beam specified as an int
<enum_for_material_type> = an int for material type:
			 - 0 = ice
			 - 1 = lead
			 - 2 = iron

You can get more info from running:
iceblock -h

TODO: CLEAN THIS UP AND MAKE IT LOOK PRETTY!