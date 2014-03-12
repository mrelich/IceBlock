#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include <fstream>
//#include "MyTreeWriter.hh"

class EventAction : public G4UserEventAction
{

 public:
  
  // Constructor/Destructor
  EventAction(std::ofstream* trkFile, std::ofstream* stepFile); //, MyTreeWriter* treeWriter);
  ~EventAction();

  // Required methods
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  std::ofstream* m_trkOutput;
  std::ofstream* m_stepOutput;
  //MyTreeWriter* m_treeWriter;

};

#endif
