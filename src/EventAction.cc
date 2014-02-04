
#include "EventAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
EventAction::EventAction()
{
  // dummy
}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
EventAction::~EventAction()
{
  // dummy
}

//-----------------------------------------------------------------//
// Begin Event
//-----------------------------------------------------------------//
void EventAction::BeginOfEventAction(const G4Event*)
{
  // dummy
}

//-----------------------------------------------------------------//
// End Event
//-----------------------------------------------------------------//
void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
    G4cout << "    " << n_trajectories 
           << " trajectories stored in this event." << G4endl;
  }
}
