#ifndef SteppingVerbose_h
#define SteppingVerbose_h

#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

class SteppingVerbose : public G4SteppingVerbose
{

 public:

  // Constructor / Destructor
  SteppingVerbose();
  ~SteppingVerbose();

  // Access step info
  void StepInfo();
  
  // Access Tracking info
  void TrackingStarted();

};

#endif
