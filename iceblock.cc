
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Simple simulation to shoot particles into a block of some material.   //
// The goal is to study shower properties in Ice, but other materials    //
// may be used to compare with previoulsy established shower properties. //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

//---------------------------------------------------//
// Main
//---------------------------------------------------//
int main() // add arguments later
{

  // My own stepping
  G4VSteppingVerbose* verboseStep = new SteppingVerbose();
  G4VSteppingVerbose::SetInstance(verboseStep);

  // Default run manager
  G4RunManager* runManager = new G4RunManager;

  // Construct detector
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Set Physics list
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);
  
  // Set Primary action to be carried out
  G4VUserPrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(genAction);
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new SteppingAction);

  // Initialize G4 Kernel
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // Start the run
  G4int numberOfEvents = 1; // TODO: Make runtime option
  runManager->BeamOn(numberOfEvents);

  delete runManager;

  return 0;

}
