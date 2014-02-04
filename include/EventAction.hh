#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

class EventAction : public G4UserEventAction
{

 public:
  
  // Constructor/Destructor
  EventAction();
  ~EventAction();

  // Required methods
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
};

#endif
