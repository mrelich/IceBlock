
#include "SetupAntenna.hh"
#include <vector>

//----------------------------------------------------//
// Constructor
//----------------------------------------------------//
SetupAntenna::SetupAntenna()
{

  // Scale factor for the timing
  G4double sf = 1.78 / 3.e8 * 1e9;

  // Trying to match some ZHS results. 
  // So put the antenna at R = 100m and
  // for several angles.
  G4double R = 100;
  
  // Angles
  std::vector<G4double> angles;
  angles.push_back(55.8197842754213767);
  //angles.push_back(55.829616);
  //angles.push_back(55.929616);
  //angles.push_back(56.029616);
  //angles.push_back(56.129616);
  //angles.push_back(56.229616);
  //angles.push_back(56.329616);
  //angles.push_back(56.429616);
  angles.push_back(56.829616);

  for(unsigned int i=0; i<angles.size(); ++i){
    G4double angle = angles.at(i) * 3.14159265358979312/180.;

    // Initialize some antennas here
    m_ants.push_back( new Antenna(R*sin(angle),
				  0,
				  R*cos(angle),
				  R * sf - 10,
				  1000,
				  0.05) );
  }// end loop over angles

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


