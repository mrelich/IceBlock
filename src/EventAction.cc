
#include "EventAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
EventAction::EventAction(std::ofstream* trkFile, 
			 std::ofstream* stepFile) :			 
  m_trkOutput(NULL),
  m_stepOutput(NULL)
  //m_treeWriter(NULL)
{
  m_trkOutput = trkFile;
  m_stepOutput = stepFile;
  //m_treeWriter = treeWriter;
}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
EventAction::~EventAction()
{
  // dummy
  m_trkOutput = NULL;
  m_stepOutput = NULL;
  //m_treeWriter = NULL;
}

//-----------------------------------------------------------------//
// Begin Event
//-----------------------------------------------------------------//
void EventAction::BeginOfEventAction(const G4Event* evt)
{
  // dummy
  (*m_trkOutput) << "Event: " << evt->GetEventID() << G4endl;
  (*m_stepOutput) << "Event: " << evt->GetEventID() << G4endl;
  //m_treeWriter->CreateEvent( evt->GetEventID() );
}

//-----------------------------------------------------------------//
// End Event
//-----------------------------------------------------------------//
void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  // periodic printing
  //if (event_id%10 == 0) {
  if (event_id%100000 == 0) {
    G4cout << "Events processed " << evt->GetEventID() << G4endl;
  }
  
  (*m_trkOutput) << "End: " << evt->GetEventID() << G4endl;
  (*m_stepOutput) << "End: " << evt->GetEventID() << G4endl;
  
  //G4cout<<"Writing event: "<<evt->GetEventID()<<G4endl;
  //m_treeWriter->WriteEvent();

}
