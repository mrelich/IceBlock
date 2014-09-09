#ifndef SetupAntennas_h
#define SetupAntennas_h

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This will be a convenient way to setup the detectors. It leaves //
// the option later to read a set of coordinates from a text file  //
// to generate the detector array.                                 //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "Antenna.hh"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

class SetupAntenna
{

 public:
  
  // Constructor
  SetupAntenna(std::string infile = "");

  // Destructor
  ~SetupAntenna();

  // Method to retrived detectors
  std::vector<Antenna*> getAnts(){ return m_ants; };

 private:
  
  std::vector<Antenna*> m_ants;

  // Method to read antenna config from files
  void readAntennaFromFile(std::string infile);

};

#endif
