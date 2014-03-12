#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//#include "MyTreeWriter.hh"
#include <fstream>

class SteppingAction : public G4UserSteppingAction
{

 public:

  // Constructor / Destructor
  SteppingAction(std::ofstream* output); //, MyTreeWriter* treeWriter);
  ~SteppingAction();
  
  // User methods -- empty for me
  void UserSteppingAction(const G4Step*);
  
 private:

  std::ofstream* m_output;  
  //MyTreeWriter* m_treeWriter;


};

#endif
