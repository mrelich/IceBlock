#ifndef RunAction_h
#define RunAction_h

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"

class RunAction : public G4UserRunAction
{

 public:
  
  RunAction();
  ~RunAction();
  
  // Required methods
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

};

#endif
