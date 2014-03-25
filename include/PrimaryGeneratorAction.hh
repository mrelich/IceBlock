#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include <string>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

 public:

  // Constructor / Destructor
  PrimaryGeneratorAction(DetectorConstruction*, G4float, std::string, G4int);
  ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun(){ return particleGun; };

 private:
  G4ParticleGun* particleGun;
  DetectorConstruction* myDetector;
};

#endif
