
#include "EventAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
EventAction::EventAction(std::ofstream* trkFile, 
			 std::ofstream* stepFile,
			 std::ofstream* vFile,
			 std::ofstream* eFile,
			 std::vector<Antenna*>* ants) :			 
  m_trkOutput(NULL),
  m_stepOutput(NULL),
  //m_treeWriter(NULL)
  m_VPotOut(NULL),
  m_EFOut(NULL),
  m_ants(NULL)
{
  m_trkOutput  = trkFile;
  m_stepOutput = stepFile;
  m_VPotOut    = vFile;
  m_EFOut      = eFile;
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

  // Write event info for track
  if(m_trkOutput)
    (*m_trkOutput)  << "Event: " << evt->GetEventID() << G4endl;
  
  // Write output for steps
  if(m_stepOutput)
    (*m_stepOutput) << "Event: " << evt->GetEventID() << G4endl;

  // Write output for Vector potential
  if(m_VPotOut){
    (*m_VPotOut)    << "Event: " << evt->GetEventID() << G4endl;
    for(unsigned int i=0; i<m_ants->size(); ++i)
      m_ants->at(i)->clear();
  }

  // Write output for Efield
  if(m_EFOut){
    (*m_EFOut)    << "Event: " << evt->GetEventID() << G4endl;
    for(unsigned int i=0; i<m_ants->size(); ++i)
      m_ants->at(i)->clear();
  }
  
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

  // Event end info for tracks
  if( m_trkOutput )
    (*m_trkOutput) << "End: " << evt->GetEventID() << G4endl;

  // Event end info for steps
  if( m_trkOutput )
    (*m_stepOutput) << "End: " << evt->GetEventID() << G4endl;
  
  //G4cout<<"Writing event: "<<evt->GetEventID()<<G4endl;
  //m_treeWriter->WriteEvent();

  // Write vector potential info
  if( m_VPotOut && m_EFOut){
    G4double Ax=0, Ay=0, Az=0;
    G4double Ex=0, Ey=0, Ez=0;
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

      (*m_EFOut) << "Antenna: "<<i<<" Pos: "
		 <<ant->getX()<<" "
		 <<ant->getY()<<" "
		 <<ant->getZ()<<" "
		 <<ant->getAngle()<<" "
		 <<ant->getRefAngle()<<" "
		 <<ant->getZprime()<<G4endl;

    
      // Now save the data points
      unsigned int nP = ant->getN();
      for(unsigned int ip=0; ip<nP; ++ip){
	ant->getPoint(ip,time,Ax,Ay,Az);
	
	(*m_VPotOut) << "\t" 
		     << time << " "
		     << Ax   << " " 
		     << Ay   << " " 
		     << Az   << G4endl;
	
      }// end loop over points

      // Save E-field points
      unsigned int nEP = ant->getEN();
      for(unsigned int ip=0; ip<nEP; ++ip){
	ant->getPoint(ip,time,Ex,Ey,Ez);
	
	(*m_EFOut) << "\t" 
		   << time << " "
		   << Ex   << " " 
		   << Ey   << " " 
		   << Ez   << G4endl;
	
      }// end loop over points

    }

    (*m_VPotOut)  << "End: " << evt->GetEventID() << G4endl;
    (*m_EFOut)    << "End: " << evt->GetEventID() << G4endl;

  }// end if VPotOut is opened file
  
}
