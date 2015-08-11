
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
					       G4double sigma,
					       G4int nbunch,
					       G4double tOffset,
					       BeamProfile* bp,
					       G4int theSeed) :
  myDetector(NULL),
  m_flat(b_flat),
  m_gauss(b_gauss),
  m_sigma(sigma),
  m_nbunch(nbunch),
  m_tOffset(tOffset),
  m_bp(NULL),
  m_theSeed(theSeed)
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
  G4RandGauss::setTheSeed(abs(theSeed + 54321));
  CLHEP::RandFlat::setTheSeed(abs(theSeed - 172634));  
  CLHEP::HepRandom::setTheSeed(m_theSeed);

  //
  // Set the beam profile object
  //
  m_bp = bp;

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

  // In addition we have the ability to create a bunch train. The
  // user should specify the bunch number and time delay between
  // bunches.  We can get more fancy and make it a time delay
  // distribution, but for right now that seems overkill.

  // Update Feb 9: Adding the option to read in beam profile that
  // is taken from the WC monitor. This is stored as a BeamProfile
  // object and is called m_bp in this class.

  // set z position of the beam
  G4double z_beam = -0.5*m;

  // Check if bp is initialized, if so use it.
  int nbunches = m_bp->isInit() ? m_bp->getN() : m_nbunch;

  // Loop over the bunches
  for(int nb = 0; nb < nbunches; ++nb){

    // Set the starting time and number of particles
    // per bunch. Again use bp info if initialized
    G4double tstart  = nb * m_tOffset * ns;
    G4int nparticles = (m_bp->isInit() ? m_bp->getQRatio(nb) : 1. ) * m_nParticles;
    
    //G4cout<<"Bunch: "<<nb<<" npart: "<<nparticles<<G4endl;
    
    // Gaussian
    if( m_gauss ){
      
      // Setup random Gaussian with some seed.
      G4double x_rand = 0.0;
      G4double y_rand = 0.0;
      //for(G4int i=0; i<nparticles; ++i){
      for(G4int i=0; i<nparticles; ++i){
	x_rand = G4RandGauss::shoot(0,m_sigma);
	y_rand = G4RandGauss::shoot(0,m_sigma);
	particleGun->SetParticlePosition(G4ThreeVector(x_rand*mm,y_rand*mm,z_beam));
	particleGun->SetParticleTime(tstart);
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
      for(G4int i=0; i<nparticles; ++i){
	
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
	
	particleGun->SetParticlePosition(G4ThreeVector(x_rand*mm,y_rand*mm,z_beam));
	particleGun->SetParticleTime(tstart);
	particleGun->GeneratePrimaryVertex(anEvent);
	//(*f_test) << x_rand <<" " << y_rand << std::endl;
	
      }// end loop over # particles
      
    }// end if flat distribution
    
    // Nominal: all at (0,0,0)
    else{
      
      G4ThreeVector startPos = G4ThreeVector(0*mm,0*mm,z_beam);
      for(G4int i=0; i<nparticles; ++i){
	particleGun->SetParticlePosition(startPos);
	particleGun->SetParticleTime(tstart);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
      
    }// End default (0,0,0) starting point for all particles
  
  }// end loop over bunch structure

}

