#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

 public:

  // Constructor / Destructor
  PrimaryGeneratorAction(DetectorConstruction*);
  ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);

 private:
  G4ParticleGun* particleGun;
  DetectorConstruction* myDetector;
  
};

#endif
