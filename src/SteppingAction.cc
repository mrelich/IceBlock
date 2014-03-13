
#include "SteppingAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingAction::SteppingAction(std::ofstream* output) :
  m_output(NULL)
  //m_treeWriter(NULL)
{
  m_output = output;
  //m_treeWriter = treeWriter;

}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
SteppingAction::~SteppingAction()
{

  m_output = NULL;
  //m_treeWriter = NULL;

}

//-----------------------------------------------------------------//
// Empty Action -- Using Stepping Verbos
//-----------------------------------------------------------------//
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  // Check the energy loss process
  int ProcessID(-100);
  if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
      == "Transportation" ) ProcessID=-1;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "msc" ) ProcessID=0;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "eIoni" ) ProcessID=1;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "eBrem" ) ProcessID=2;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "annihil" ) ProcessID=3;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "phot" ) ProcessID=4;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "compt" ) ProcessID=5;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "conv" ) ProcessID=6;
  else if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
          == "muIoni" ) ProcessID = 1;
  else if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	  == "UserSpecialCut" ) return;
  else ProcessID=-100;

  
  // Get step properties
  float x  = aStep->GetPreStepPoint()->GetPosition().x()/cm;
  float y  = aStep->GetPreStepPoint()->GetPosition().y()/cm;
  float z  = aStep->GetPreStepPoint()->GetPosition().z()/cm;
  float dE = aStep->GetTotalEnergyDeposit()/MeV;
  float dX = aStep->GetStepLength()/cm;
  
  // Get the current particle
  // add a step
  //Particle* part = m_treeWriter->GetEvt()->getLastPart();
  //part->addStep( Step(x,y,z,
  //dE,dX,
  //ProcessID) 
  //);   
  

  G4Track* track = aStep->GetTrack();
  
  (*m_output) <<
    x << " " <<
    y << " " <<
    z << " " <<
    dX << " " <<
    track->GetKineticEnergy() / MeV << " " <<
    dE << " " <<
    track->GetParticleDefinition()->GetPDGEncoding()<< " " <<
    track->GetTrackID() << " " <<
    track->GetParentID() << " " <<
    ProcessID << " " <<
    G4endl;

}
