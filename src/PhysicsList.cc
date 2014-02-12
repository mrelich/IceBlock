
#include "PhysicsList.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
PhysicsList::PhysicsList()
{

}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
PhysicsList::~PhysicsList()
{

}

//-----------------------------------------------------------------//
// Construct particles
//-----------------------------------------------------------------//
void PhysicsList::ConstructParticle()
{

  // Need the geantino for transportation
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // Initialize bosons
  ConstructBosons();

  // Initialize leptons
  ConstructLeptons();

  // Initialize hadrons (commented out for now)
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

}

//-----------------------------------------------------------------//
// Electromagnetic processes
//-----------------------------------------------------------------//
void PhysicsList::ConstructEM()
{

  // Copying from Example02 in novice geant4 examples
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    
    // NOTE:
    // This should only be set if we are running on ice
    // in order to match the papers results.
    // Update: Only turned on for gamma/e/p so as to avoid print out

    
    G4String particleName = particle->GetParticleName();

    // Gamma
    if (particleName == "gamma" ) {
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);
      particle->SetApplyCutsFlag(true);
    }
    // Electron
    else if (particleName == "e-") {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);      
      particle->SetApplyCutsFlag(true);
    }
    // Positron
    else if (particleName == "e+" ) {
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);
      ph->RegisterProcess(new G4eplusAnnihilation,   particle);
      particle->SetApplyCutsFlag(true);
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
  
  //SetCutValue(0.1*mm, "gamma");
  //SetCutValue(0.1*mm, "e+");
  //SetCutValue(0.1*mm, "e-");

  DumpCutValuesTable();



}
