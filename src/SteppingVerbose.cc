
#include "SteppingVerbose.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingVerbose::SteppingVerbose()
{

  outputFile.open("outputfile.dat",ofstream::out);

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

  //
  // Dump some basic info
  //

  outputFile <<
    fTrack->GetVolume()->GetName()<<" "<<
    fStep->GetPostStepPoint()->GetPosition().x()/cm <<" "<<
    fStep->GetPostStepPoint()->GetPosition().y()/cm <<" "<<
    fStep->GetPostStepPoint()->GetPosition().z()/cm <<" "<<
    fTrack->GetKineticEnergy()/MeV<<" "<<
    fTrack->GetParticleDefinition()->GetPDGEncoding()<<" "<<
    fTrack->GetTrackID()<<" "<<
    fTrack->GetParentID()<<" "<<
    G4endl;
  
}

//-----------------------------------------------------------------//
// Track info
//-----------------------------------------------------------------//
void SteppingVerbose::TrackingStarted()
{

  // Necessary
  CopyState();

}
