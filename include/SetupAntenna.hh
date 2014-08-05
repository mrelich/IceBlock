#ifndef SetupAntennas_h
#define SetupAntennas_h

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This will be a convenient way to setup the detectors. It leaves //
// the option later to read a set of coordinates from a text file  //
// to generate the detector array.                                 //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "Antenna.hh"
#include <vector>

class SetupAntenna
{

 public:
  
  // Constructor
  SetupAntenna();

  // Destructor
  ~SetupAntenna();

  // Method to retrived detectors
  std::vector<Antenna*> getAnts(){ return m_ants; };

 private:
  
  std::vector<Antenna*> m_ants;

};

#endif
