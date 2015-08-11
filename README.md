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
iceblock -ne <number_of_events> -e <beam_energy> -np <number_of_primaries>

<number_of_events> = number of events specified as an int
<beam_energy> = energy of the particle beam specified as an int
<enum_for_material_type> = an int for material type:
			 - 0 = ice
			 - 1 = lead
			 - 2 = iron

You can get more info from running:
iceblock -h

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
=                      UPDATES AS OF AUGUST 2015                        =
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Since this README was written we have made the following additions:

1.) Beam Shape -- we now allow one to set the beam shape to either:
    	 * Single injection point (default)
	 * 2-D Gaussian with some width by using: --gauss 45 (for 45mm)
	 * Flat uniform distribution by using: --flat 45 (squre with dimension 45mm)

2.) Bunch Number -- This will allow one to specify the number of bunches. Please
    note that when this option is used, it creates bunches with the same number
    of primaries specified at np.  So for example -np 1000 --nbunch 5 will create
    5 bunches with 1000 primaries

3.) Timing Offset -- This allows one to control the timing difference between 
    bunches in the bunch train.  The default is 350ps, as the beam in Utah

4.) Beam Profile -- We allow one to specify a beam profile.  What that means is we
    allow one to imput a more complicated bunch structure instead of equal bunches.
    This allows control of the relative particle number in each bunch.  One thing
    to note: The total number of particles set using -np is distributed over the
    input bunch structure.  A simple case of 3 bunches, where 1st has 10% of 
    the charge, and 2nd has 50% and 3rd has 40%.  So for -np 100 you would 
    get: 10, 50, and 40 particles in the 1st, 2nd, and 3rd bunches respectively.

5.) Step Limit -- This will force Geant4 to take a step even time, regardless 
    if any physics interaction is taking place.  This is useful for the Endpoint
    method, I believe, and that is why it was added.  ZHS doesn't care about it.

6.) Antenna file -- The user can specify a text file for the Antenna positions.  
    Unfortunately I haven't had time to clean this up, but in reality the only 
    important numbers are x, y, z.  This should be cleaned up in the future
    and potentially merged with the lookup table inputs.


#-------------------------------------------------------------#
# Current Issues that could be cleaned up or dealt with
#-------------------------------------------------------------#

1.) Right now there is no way to easily run without refraction and with 
    refraction turned on.

2.) We may need a better way to handle lookup tables, or at least
    generate lookup tables for finer angular spacings. This will
    be important for characterizing in simulation the angular
    dependence of TR and Askaryan in our Delta Ice

3.) Maybe think of adding a command line option for changing index of 
    refraction.  This might be a source of systematic error in the end.

4.) We need a robust way to test the lookup tables.  Some parameters
    or code to verify that they are ok to be used for both the outside
    points and the inside points.
