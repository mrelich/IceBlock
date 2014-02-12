#ifndef SteppingVerbose_h
#define SteppingVerbose_h

#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include <fstream>

using namespace std;

class SteppingVerbose : public G4SteppingVerbose
{

 public:

  // Constructor / Destructor
  SteppingVerbose(string outputName);
  ~SteppingVerbose();

  // Access step info
  void StepInfo();
  
  // Access Tracking info
  void TrackingStarted();

 protected:

  ofstream outputFile;   // temporary output file

  G4int m_createdProc;
  G4int m_event;
  
};

#endif
