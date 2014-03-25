
#include "PrimaryGeneratorAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC, 
					       G4float partEnergy,
					       std::string partType,
					       G4int n_particle) :
  myDetector(NULL)
{

  //
  // Specify detector
  //
  myDetector = myDC;

  //
  // Specify the number of particles to be simulated
  //
  //G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  //
  // specify the particle to be used
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(partType.c_str());

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); // z-direction
  particleGun->SetParticleEnergy(partEnergy);
  

}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{

  delete particleGun;

}

//-----------------------------------------------------------------//
// Generate the primaries
//-----------------------------------------------------------------//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Right now my block of ice is set such that it is in the
  // right side of the world volume.  World volume is 2*km cube
  // and the iceblock is a 1*km cube set such that ice is offset
  // So particles should be input at 0,0,0
  particleGun->SetParticlePosition(G4ThreeVector(0*mm,0*mm,0*mm));  
  particleGun->GeneratePrimaryVertex(anEvent);

}

