
#include "EventAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
EventAction::EventAction(std::ofstream* trkFile, 
			 std::ofstream* stepFile,
			 std::ofstream* vFile,
			 std::vector<Antenna*>* ants) :			 
  m_trkOutput(NULL),
  m_stepOutput(NULL),
  //m_treeWriter(NULL)
  m_VPotOut(NULL),
  m_ants(NULL)
{
  m_trkOutput  = trkFile;
  m_stepOutput = stepFile;
  m_VPotOut    = vFile;
  m_ants       = ants;
}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
EventAction::~EventAction()
{
  // dummy
  m_trkOutput  = NULL;
  m_stepOutput = NULL;
  m_VPotOut    = NULL;
  m_ants       = NULL;
}

//-----------------------------------------------------------------//
// Begin Event
//-----------------------------------------------------------------//
void EventAction::BeginOfEventAction(const G4Event* evt)
{
  // dummy
  (*m_trkOutput)  << "Event: " << evt->GetEventID() << G4endl;
  (*m_stepOutput) << "Event: " << evt->GetEventID() << G4endl;
  (*m_VPotOut)    << "Event: " << evt->GetEventID() << G4endl;
  for(unsigned int i=0; i<m_ants->size(); ++i)
    m_ants->at(i)->clear();
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

  // Write vector potential info
  G4double Ax=0, Ay=0, Az=0;
  G4double time = 0;
  for(unsigned int i=0; i<m_ants->size(); ++i){
    Antenna* ant = m_ants->at(i);

    (*m_VPotOut) << "Antenna: "<<i<<" Pos: "
		 <<ant->getX()<<" "
		 <<ant->getY()<<" "
		 <<ant->getZ()<<" "
		 <<ant->getAngle()<<" "
		 <<ant->getRefAngle()<<" "
		 <<ant->getZprime()<<G4endl;

    
    // Now save the data points
    unsigned int nP = ant->getN();
    //G4double step   = ant->getTStep() * 1e-9; // Put into [s]
    for(unsigned int ip=0; ip<nP; ++ip){
      ant->getPoint(ip,time,Ax,Ay,Az);

      (*m_VPotOut) << "\t" 
		   << time << " "
		   << Ax   << " " 
		   << Ay   << " " 
		   << Az   << G4endl;

    }// end loop over points

  }
  (*m_VPotOut)    << "End: " << evt->GetEventID() << G4endl;

  
}
