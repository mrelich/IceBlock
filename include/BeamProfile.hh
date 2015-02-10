#ifndef BeamProfile_h
#define BeamProfile_h

#include "globals.hh"
#include <vector>
#include <string>
#include <fstream>

class BeamProfile
{

 public:
  
  // Constructor
  BeamProfile(){
    m_Qratios.clear();
    m_QratErrs.clear();
    m_initialized = false;
  };

  // Destructor
  ~BeamProfile(){
    m_Qratios.clear();
    m_QratErrs.clear();
  };
  
  // Init the beam profile
  void init(std::string inputfile){
    
    // Open input file
    std::ifstream infile(inputfile.c_str());

    // Set bool to initialized if file is good
    m_initialized = infile.good();

    // Loop over infile and get ratio
    G4int bunch     = 0;
    G4int prevBunch = -999;
    G4double ratio  = 0;
    G4double error  = 0;
    while( !infile.eof() ){
      infile >> bunch >> ratio >> error;
      if(prevBunch == bunch) continue;
      prevBunch = bunch;
      m_Qratios.push_back(ratio);
      m_QratErrs.push_back(error);
      G4cout<<"Adding: "<<ratio<<" "<<error<<G4endl;
    }

  };

  // Get Number of bunches
  unsigned int getN(){ return m_Qratios.size(); };

  // Get the charge ratio for given bunch
  G4double getQRatio(unsigned int i){ return m_Qratios.at(i); };
  
  // Get the charge ratio error
  G4double getQratErrs(unsigned int i){ return m_QratErrs.at(i); };

  // Get vector of Qratios
  std::vector<G4double> getQRatios(){ return m_Qratios; };

  // Get vector of Q ratio errors
  std::vector<G4double> getQratErrs(){ return m_QratErrs; };

  // Check if tool initialized
  G4bool isInit(){ return m_initialized; };

 private:
  
  std::vector<G4double> m_Qratios;   // Q ratio for structure
  std::vector<G4double> m_QratErrs;  // Error on Q ratio
  
  G4bool m_initialized;              // whether or not initialized

};

#endif
