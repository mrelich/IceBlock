#ifndef PhysicsList_h
#define PhysicsList_h

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

// Following includes copied from Example02
// in novice examples from geant directory
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

class PhysicsList : public G4VUserPhysicsList
{

 public:
  
  PhysicsList();
  ~PhysicsList();

 protected:
  
  // Construct particles
  void ConstructParticle();
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructHadrons();

  // Construct Processes
  void ConstructProcess();
  void ConstructEM();

  // Set cuts
  void SetCuts();


};

#endif
