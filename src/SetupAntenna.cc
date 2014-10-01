
#include "SetupAntenna.hh"

//----------------------------------------------------//
// Constructor
//----------------------------------------------------//
SetupAntenna::SetupAntenna(std::string infile)
{

  // If an antenna file is passed in, 
  // use that to configure.  Otherwise
  // use the default setup.
  if( !infile.empty() ){
    readAntennaFromFile(infile);
    return;
  }


  // Scale factor for the timing
  G4double sf = 1.78 / 3.e8 * 1e9;

  // Trying to match some ZHS results. 
  // So put the antenna at R = 100m and
  // for several angles.
  G4double R = 1000;
  //G4double R = 10;

  // Angles
  std::vector<G4double> angles;
  //angles.push_back(54.4);
  //angles.push_back(54.5);
  angles.push_back(55.829616);
  /*angles.push_back(54.6);
  angles.push_back(54.7);
  angles.push_back(54.8);
  angles.push_back(54.9);
  angles.push_back(55.0);
  angles.push_back(55.1);
  angles.push_back(55.2);
  angles.push_back(55.3);
  angles.push_back(55.4);
  angles.push_back(55.5);
  angles.push_back(55.6);
  angles.push_back(55.7);
  angles.push_back(55.8);
  */

  for(unsigned int i=0; i<angles.size(); ++i){
    G4double angle = angles.at(i) * 3.14159265358979312/180.;

    // Initialize some antennas here
    m_ants.push_back( new Antenna(R*sin(angle),
				  0,
				  R*cos(angle),
				  R * sf - 10,
				  2000,
				  0.01,
				  0,0,0));
  }// end loop over angles

  /*
  G4double R = 100; //115.602;
  G4double angle = 55.829616 * 3.14159265358979312/180.;
  m_ants.push_back( new Antenna(R*sin(angle),
				0,      
				R*cos(angle),
				R * sf - 10,    
				2000,
				0.01) ); 
  */

  // Testing old configuration
  //G4double x = 6;
  //G4double z = 4.04;
  //G4double R = sqrt(x*x+z*z);
  //m_ants.push_back( new Antenna(x,0,z,
  //R*sf-10,
  //				2000,
  //				0.01) );

}


//----------------------------------------------------//
// Destructor
//----------------------------------------------------//
SetupAntenna::~SetupAntenna()
{

  // Clean up detectors
  for(unsigned int i=0; i<m_ants.size(); ++i)
    delete m_ants.at(i);
  m_ants.clear();

}

//----------------------------------------------------//
// Read antenna configuration from an input file
//----------------------------------------------------//
void SetupAntenna::readAntennaFromFile(std::string infile)
{

  std::ifstream input(infile.c_str());

  // The file structure should be the following:
  // x y z

  // Scale factor for starting timing
  G4double sf = 1.78 / 2.99e8 * 1e9;

  // I will fix the number of points to 1000 
  // with 0.5 ns steps
  int np         = 2000;
  G4double stepSize = 0.01; // ns
  
  // Antenna positions
  G4double x = 0;
  G4double y = 0;
  G4double z = 0;
  G4double angle    = 0;
  G4double refAngle = 0;
  G4double zprime   = 0;

  // hack to not double count last line
  G4double prevX = 0;
  G4double prevY = 0;
  G4double prevZ = 0;
  while( input.good() ){

    input >> x >> y >> z >> angle >> refAngle >> zprime;
    
    // Not the best way (comparing doubles) but
    // try this for now...
    if( prevX == x && prevY == y && prevZ == z ) continue;
    prevX = x;
    prevY = y;
    prevZ = z;

    // Don't need to calculate anymore. Stored in file
    G4double R = sqrt(x*x+y*y+z*z);

    G4cout<<"Setting antenna: "<<x<<" "<<y<<" "<<z<<" "<<R
	  <<" "<<angle<<" "<<refAngle<<" "<<zprime<<G4endl;
    m_ants.push_back( new Antenna(x,y,z,
                                  R * sf - 10,
                                  np,
                                  stepSize,
				  angle,
				  refAngle,
				  zprime) );

  }
    
  G4cout<<"Setup antenna"<<G4endl;
  input.close();

}



