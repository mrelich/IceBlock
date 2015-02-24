
#include "PhysicsList.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
PhysicsList::PhysicsList(bool useThresh) :
  theCerenkovProcess(NULL),
  theScintillationProcess(NULL),
  theAbsorptionProcess(NULL),
  theRayleighScatteringProcess(NULL),
  theMieHGScatteringProcess(NULL),
  theBoundaryProcess(NULL),
  m_useThreshold(false)
{
  emStandard = new G4EmStandardPhysics();
  m_useThreshold = useThresh;
}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
PhysicsList::~PhysicsList()
{

  delete emStandard;

}

//-----------------------------------------------------------------//
// Construct particles
//-----------------------------------------------------------------//
void PhysicsList::ConstructParticle()
{

  // Need the geantino for transportation
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // Try to see if error in my EM def
  // UPDATE: Doesn't have any difference.
  //emStandard->ConstructParticle();

  // Initialize bosons
  //ConstructBosons();

  // Initialize leptons
  //ConstructLeptons();

  // Initialize hadrons 
  // Commented out for now, only looking at EM shower
  //ConstructHadrons();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

  // mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  // baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

}

//-----------------------------------------------------------------//
// Construct bosons
//-----------------------------------------------------------------//
void PhysicsList::ConstructBosons()
{

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  //G4OpticalPhoton::OpticalPhotonDefinition();

}

//-----------------------------------------------------------------//
// Construct Leptons
//-----------------------------------------------------------------//
void PhysicsList::ConstructLeptons()
{

  // Electron
  G4Electron::ElectronDefinition();

  // Positron
  G4Positron::PositronDefinition();

  // Muon
  // Needed for using Muon Beam
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

}

//-----------------------------------------------------------------//
// Construct Hadrons
//-----------------------------------------------------------------//
void PhysicsList::ConstructHadrons()
{

  // Proton
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  
  // Neutron  
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // Meson
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  // Ions
  G4Alpha::Alpha(); 
  G4Deuteron::DeuteronDefinition();
  G4GenericIon::GenericIonDefinition();
  G4He3::He3Definition();
  G4Triton::TritonDefinition();

}

//-----------------------------------------------------------------//
// Construct processes
//-----------------------------------------------------------------//
void PhysicsList::ConstructProcess()
{
  
  // Make sure to have transportation
  AddTransportation();

  // Turn on processes here.  Right now ONLY considering
  // the electromagnetic stuff.
  ConstructEM();
  AddDecay();

  // For a cross-check on my EM
  //emStandard->ConstructProcess();

}



//-----------------------------------------------------------------//
// Particle decay
//-----------------------------------------------------------------//
void PhysicsList::AddDecay()
{
  // Add Decay Process

  G4Decay* fDecayProcess = new G4Decay();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) { 

      pmanager ->AddProcess(fDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(fDecayProcess, idxAtRest);

    }
  }
}


//-----------------------------------------------------------------//
// Electromagnetic processes
//-----------------------------------------------------------------//
void PhysicsList::ConstructEM()
{


  // Copying from Example02 in novice geant4 examples
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  if( m_useThreshold ){
    ph->RegisterProcess(new G4UserSpecialCuts(), G4Gamma::GammaDefinition());
    ph->RegisterProcess(new G4UserSpecialCuts(), G4Electron::ElectronDefinition());
    ph->RegisterProcess(new G4UserSpecialCuts(), G4Positron::PositronDefinition());
  }

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();        
    G4String particleName = particle->GetParticleName();
    //G4ProcessManager* pmanager = particle->GetProcessManager();
    //pmanager->AddDiscreteProcess(new G4StepLimiter);
    //ph->RegisterProcess(new G4UserSpecialCuts(), particle);

    // Gamma
    if (particleName == "gamma" ) {
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      G4ComptonScattering* cs = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());
      ph->RegisterProcess(cs, particle);
      //ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);

      //particle->SetApplyCutsFlag(true);
      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
      //ph->RegisterProcess(new G4UserSpecialCuts(), particle);
    }
    // Electron
    else if (particleName == "e-") {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      // extra param for ionization...
      G4eIonisation* eIoni = new G4eIonisation();
      //eIoni->SetStepFunction(0.1,100*um);
      ph->RegisterProcess(eIoni, particle);
      //ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);      

      //particle->SetApplyCutsFlag(true);
      //ph->RegisterProcess(new G4UserSpecialCuts(), particle);
      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
    }
    // Positron
    else if (particleName == "e+" ) {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      G4eIonisation* eIoni = new G4eIonisation();
      //eIoni->SetStepFunction(0.1, 100*um);      
      ph->RegisterProcess(eIoni, particle);
      //ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);
      ph->RegisterProcess(new G4eplusAnnihilation,   particle);
      //particle->SetApplyCutsFlag(true);
      //ph->RegisterProcess(new G4UserSpecialCuts(), particle);
      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
    } 
    // Muons
    else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      ph->RegisterProcess(new G4MuMultipleScattering, particle);
      ph->RegisterProcess(new G4MuIonisation,         particle);
      ph->RegisterProcess(new G4MuBremsstrahlung,     particle);
      ph->RegisterProcess(new G4MuPairProduction,     particle);
      
    } 
    // Protons
    else if( particleName == "proton" || 
	     particleName == "pi-" ||
	     particleName == "pi+"    ) {
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);
      ph->RegisterProcess(new G4hBremsstrahlung,     particle);
      ph->RegisterProcess(new G4hPairProduction,     particle);       
      
    }
    // More complex stuff
    else if( particleName == "alpha" || 
               particleName == "He3"  )     {
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);
      
    }
    // Some ion
    else if( particleName == "GenericIon" ) { 
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);     
      
    }
    // otherwise
    else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino" )) {
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);        
    }     
  }
  
  
}


//-----------------------------------------------------------------//
// Set Cuts
//-----------------------------------------------------------------//
void PhysicsList::SetCuts()
{

  SetCutsWithDefault();

  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(0.611*MeV, 1*TeV);
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(5*MeV, 10*TeV);
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*MeV, 10*TeV);
  
  //G4double cutval = 7*um;
  //SetCutValue(cutval, "gamma");
  //SetCutValue(cutval, "e-");
  //SetCutValue(cutval, "e+");

  DumpCutValuesTable();



}
