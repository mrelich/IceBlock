
#include "PrimaryGeneratorAction.hh"


//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC, 
					       G4float partEnergy,
					       std::string partType,
					       G4int n_particle,
					       G4bool b_flat,
					       G4bool b_gauss,
					       G4double sigma) :
  myDetector(NULL),
  m_flat(b_flat),
  m_gauss(b_gauss),
  m_sigma(sigma)
{

  //
  // Specify detector
  //
  myDetector = myDC;

  //
  // Specify the number of particles to be simulated
  //
  //G4int n_particle = 1;
  //particleGun = new G4ParticleGun(n_particle);
  particleGun  = new G4ParticleGun(1);
  m_nParticles = n_particle; 

  //
  // specify the particle to be used
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(partType.c_str());

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); // z-direction
  particleGun->SetParticleEnergy(partEnergy);

  //f_test = new std::ofstream("testingDist.txt",std::ofstream::out);
  //f_test = new std::ofstream("beamProfile/Gauss.txt",std::ofstream::out);
  //f_test = new std::ofstream("beamProfile/Flat.txt",std::ofstream::out);
  
  //
  // Specify the seeds
  //
  G4RandGauss::setTheSeed(1234567654);
  CLHEP::RandFlat::setTheSeed(1234567654);

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

  // We now have several ways to place the particle's vertex.
  // Implemented below are random gaussian distribution, random
  // flat distribution with some width, and the original way,
  // which is to place them all at (0,0,0).

  // Gaussian
  if( m_gauss ){
    
    // Setup random Gaussian with some seed.
    G4double x_rand = 0.0;
    G4double y_rand = 0.0;
    //for(G4int i=0; i<m_nParticles; ++i){
    for(G4int i=0; i<m_nParticles; ++i){
      x_rand = G4RandGauss::shoot(0,m_sigma);
      y_rand = G4RandGauss::shoot(0,m_sigma);
      particleGun->SetParticlePosition(G4ThreeVector(x_rand*mm,y_rand*mm,0*mm));
      particleGun->GeneratePrimaryVertex(anEvent);
      //(*f_test) << x_rand <<" " << y_rand << std::endl;
    }// end loop over # particles

  }// end if gauss

  // Flat distribution
  else if( m_flat ){

    // Get uniform random number
    G4double x_rand  = 0;
    G4double y_rand  = 0;
    G4double x_sign  = 1;
    G4double y_sign  = 1;
    for(G4int i=0; i<m_nParticles; ++i){
      
      x_rand = 999;
      y_rand = 999;
      x_sign = 1;
      y_sign = 1;
      if( G4UniformRand() < 0.5 ) x_sign = -1;
      if( G4UniformRand() < 0.5 ) y_sign = -1;
      
      // Make sure radius is <= m_sigma
      while( sqrt( x_rand*x_rand + y_rand*y_rand ) > m_sigma ){
	x_rand = m_sigma * G4UniformRand();
	y_rand = m_sigma * G4UniformRand();
      }
      x_rand *= x_sign;
      y_rand *= y_sign;

      particleGun->SetParticlePosition(G4ThreeVector(x_rand*mm,y_rand*mm,0*mm));
      particleGun->GeneratePrimaryVertex(anEvent);
      //(*f_test) << x_rand <<" " << y_rand << std::endl;

    }// end loop over # particles

  }// end if flat distribution
  
  // Nominal: all at (0,0,0)
  else{

    G4ThreeVector startPos = G4ThreeVector(0*mm,0*mm,0*mm);
    for(G4int i=0; i<m_nParticles; ++i){
      particleGun->SetParticlePosition(startPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }

  }// End default (0,0,0) starting point for all particles
  
}

