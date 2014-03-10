
#include "SteppingVerbose.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingVerbose::SteppingVerbose(string outputName) :
  m_createdProc(-1),
  m_event(-1)
{
  
  string dir = string( getenv("PWD") );
  //outputFile.open((dir+"/output/"+outputName+".dat").c_str(),ofstream::out);
  //outputFile.open((dir+"/tmp/"+outputName+".dat").c_str(),ofstream::out);
  outputFile.open((dir+"/output/"+outputName+".dat").c_str(),ofstream::out);

}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
SteppingVerbose::~SteppingVerbose()
{

  outputFile.close();

}

//-----------------------------------------------------------------//
// Get step info
//-----------------------------------------------------------------//
void SteppingVerbose::StepInfo()
{

  // Necessary
  CopyState();

  if( fTrack->GetParticleDefinition()->GetPDGEncoding() == 0 ) return;

  //
  // Dump some basic info
  //

  int ProcessID(-100);
  if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
      == "Transportation" ) ProcessID=-1;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
           == "msc" ) ProcessID=0;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
           == "eIoni" ) ProcessID=1;      
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "eBrem" ) ProcessID=2;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "annihil" ) ProcessID=3;  
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "phot" ) ProcessID=4;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "compt" ) ProcessID=5;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "conv" ) ProcessID=6;
  else if(fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	  == "muIoni" ) ProcessID = 1;
  //else return; //  ProcessID=-100;
  else ProcessID=-100;

  outputFile <<
    //fStep->GetPostStepPoint()->GetPosition().x()/cm <<" "<<
    //fStep->GetPostStepPoint()->GetPosition().y()/cm <<" "<<
    //fStep->GetPostStepPoint()->GetPosition().z()/cm <<" "<<
    fStep->GetPreStepPoint()->GetPosition().x()/cm <<" "<<
    fStep->GetPreStepPoint()->GetPosition().y()/cm <<" "<<
    fStep->GetPreStepPoint()->GetPosition().z()/cm <<" "<<
    fStep->GetStepLength()/cm <<" "<<
    fTrack->GetKineticEnergy()/MeV<<" "<<
    fStep->GetTotalEnergyDeposit()/MeV<<" "<<
    //(fStep->GetTotalEnergyDeposit() - fStep->GetNonIonizingEnergyDeposit())/MeV<<" "<<
    fTrack->GetParticleDefinition()->GetPDGEncoding()<<" "<<
    fTrack->GetTrackID()<<" "<<
    fTrack->GetParentID()<<" "<<
    m_createdProc<<" "<<
    ProcessID<<" "<<
    //fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<" "<<
    G4endl;
}

//-----------------------------------------------------------------//
// Track info
//-----------------------------------------------------------------//
void SteppingVerbose::TrackingStarted()
{

  // Necessary
  CopyState();

  if( fTrack->GetParticleDefinition()->GetPDGEncoding() == 0 ) return;

  // The Created Process will crash for the intial particle.
  // So make sure we are grabbing something that exists
  if( fTrack->GetTrackID() != 1 ){

    if( fTrack->GetCreatorProcess()->GetProcessName()=="msc" )   m_createdProc = 0;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="eIoni" ) m_createdProc = 1;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="eBrem" ) m_createdProc = 2;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="phot" )  m_createdProc = 3;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="compt" ) m_createdProc = 4;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="conv" )  m_createdProc = 5;
  }
  else m_createdProc = -1;


  // Keep track of the event number such that we can 
  // read the file back in easily later to make tree
  // TODO: Maybe just make root tree now..?
  if( fTrack->GetTrackID() == 1 ){
    if( m_event >= 0 )
      outputFile << "End" << G4endl;
    m_event++;
    outputFile << "Event: " << m_event << G4endl;
  }

}

