
#include "PhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4UserSpecialCuts.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
PhysicsList::PhysicsList()
{

  emStandard = new G4EmStandardPhysics();

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
  emStandard->ConstructParticle();

  // Initialize bosons
  ConstructBosons();

  // Initialize leptons
  ConstructLeptons();

  // Initialize hadrons 
  // Commented out for now, only looking at EM shower
  //ConstructHadrons();
  
}

//-----------------------------------------------------------------//
// Construct bosons
//-----------------------------------------------------------------//
void PhysicsList::ConstructBosons()
{

  // gamma
  G4Gamma::GammaDefinition();

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
  
  // Make sure to have transportation : )
  AddTransportation();

  // Turn on processes here.  Right now ONLY considering
  // the electromagnetic stuff.
  ConstructEM();

  // For a cross-check on my EM
  //emStandard->ConstructProcess();

}

//-----------------------------------------------------------------//
// Electromagnetic processes
//-----------------------------------------------------------------//
void PhysicsList::ConstructEM()
{


  // Copying from Example02 in novice geant4 examples
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  ph->RegisterProcess(new G4UserSpecialCuts(), G4Gamma::GammaDefinition());
  ph->RegisterProcess(new G4UserSpecialCuts(), G4Electron::ElectronDefinition());
  ph->RegisterProcess(new G4UserSpecialCuts(), G4Positron::PositronDefinition());
  

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();        
    G4String particleName = particle->GetParticleName();
    //G4ProcessManager* pmanager = particle->GetProcessManager();

    // Gamma
    if (particleName == "gamma" ) {
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);
      particle->SetApplyCutsFlag(true);
      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
      //ph->RegisterProcess(new G4UserSpecialCuts(), particle);
    }
    // Electron
    else if (particleName == "e-") {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);      
      particle->SetApplyCutsFlag(true);
      //ph->RegisterProcess(new G4UserSpecialCuts(), particle);
      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
    }
    // Positron
    else if (particleName == "e+" ) {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);
      ph->RegisterProcess(new G4eplusAnnihilation,   particle);
      particle->SetApplyCutsFlag(true);
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

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(0.611*MeV, 1*TeV);
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(5*MeV, 10*TeV);
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*MeV, 10*TeV);
  
  G4double cutval = 0.1*mm;
  SetCutValue(cutval, "gamma");
  SetCutValue(cutval, "e-");
  SetCutValue(cutval, "e+");

  DumpCutValuesTable();



}
