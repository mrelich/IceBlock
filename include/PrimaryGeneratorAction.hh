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
#include "Randomize.hh"

#include <fstream>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

 public:

  // Constructor / Destructor
  PrimaryGeneratorAction(DetectorConstruction*, G4float, std::string, G4int,
			 G4bool, G4bool, G4double);
  ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun(){ return particleGun; };

 private:
  G4ParticleGun* particleGun;
  DetectorConstruction* myDetector;
  G4int m_nParticles;

  G4bool m_flat;       // True for random flat distribution
  G4bool m_gauss;      // True for random gauss distribution
  G4double m_sigma;    // Describe width of two distributions

  //std::ofstream* f_test; // Used to get data on beam profile
};

#endif
