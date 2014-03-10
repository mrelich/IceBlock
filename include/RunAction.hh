#ifndef RunAction_h
#define RunAction_h

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

class RunAction : public G4UserRunAction
{

 public:
  
  //RunAction();
  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  ~RunAction();

  void GetCuts();
  void CriticalEnergy();
  
  // Required methods
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fPrimary;
  G4double  fRangeCut[3];
  G4double fEnergyCut[3];

};

#endif
