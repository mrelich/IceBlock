#ifndef PhysicsList_h
#define PhysicsList_h

#include "G4Decay.hh"

#include "G4KleinNishinaModel.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4UserSpecialCuts.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

// Following includes copied from Example02
// in novice examples from geant directory
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"


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

#include "G4EmStandardPhysics.hh"

class PhysicsList : public G4VUserPhysicsList
{

 public:
  
  PhysicsList(bool useThresh);
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
  void ConstructOp();
  void AddDecay();

  // Set cuts
  void SetCuts();

  G4EmStandardPhysics* emStandard;

  G4Cerenkov*          theCerenkovProcess;
  G4Scintillation*     theScintillationProcess;
  G4OpAbsorption*      theAbsorptionProcess;
  G4OpRayleigh*        theRayleighScatteringProcess;
  G4OpMieHG*           theMieHGScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;
  bool m_useThreshold;
  
};

#endif
