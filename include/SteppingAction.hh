#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "G4SteppingManager.hh"

class SteppingAction : public G4UserSteppingAction
{

 public:

  // Constructor / Destructor
  SteppingAction();
  ~SteppingAction();
  
  // User methods -- empty for me
  void UserSteppingAction(const G4Step*);
  

};

#endif
