
#include "SteppingAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingAction::SteppingAction(std::ofstream* output,
			       std::vector<Antenna*> *ants) :
  m_output(NULL),
  //m_treeWriter(NULL)
  m_ants(NULL)
{
  m_output = output;
  //m_treeWriter = treeWriter;
  m_ants = ants;
}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
SteppingAction::~SteppingAction()
{

  m_output = NULL;
  //m_treeWriter = NULL;
  m_ants = NULL;

}

//-----------------------------------------------------------------//
// Turn on and off here what we want to run
//-----------------------------------------------------------------//
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  //WriteSteps(aStep);
  VPotentialZHSStyle(aStep);
}

//-----------------------------------------------------------------//
// User Stepping Action
//-----------------------------------------------------------------//
void SteppingAction::WriteSteps(const G4Step* aStep)
{

  // Check the energy loss process
  int ProcessID(-100);
  if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
      == "Transportation" ) ProcessID=-1;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "msc" ) ProcessID=0;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "eIoni" ) ProcessID=1;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "eBrem" ) ProcessID=2;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "annihil" ) ProcessID=3;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "phot" ) ProcessID=4;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "compt" ) ProcessID=5;
  else if( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
           == "conv" ) ProcessID=6;
  else if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
          == "muIoni" ) ProcessID = 1;
  else if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	  == "UserSpecialCut" ) return;
  else ProcessID=-100;

  
  // Get step properties
  float x  = aStep->GetPreStepPoint()->GetPosition().x()/cm;
  float y  = aStep->GetPreStepPoint()->GetPosition().y()/cm;
  float z  = aStep->GetPreStepPoint()->GetPosition().z()/cm;
  float dE = aStep->GetTotalEnergyDeposit()/MeV;
  float dX = aStep->GetStepLength()/cm;
  
  // Get the current particle
  // add a step
  //Particle* part = m_treeWriter->GetEvt()->getLastPart();
  //part->addStep( Step(x,y,z,
  //dE,dX,
  //ProcessID) 
  //);   
  

  G4Track* track = aStep->GetTrack();
  
  (*m_output) <<
    x << " " <<
    y << " " <<
    z << " " <<
    dX << " " <<
    track->GetKineticEnergy() / MeV << " " <<
    dE << " " <<
    track->GetParticleDefinition()->GetPDGEncoding()<< " " <<
    track->GetTrackID() << " " <<
    track->GetParentID() << " " <<
    ProcessID << " " <<
    G4endl;

}

//-----------------------------------------------------------------//
// Calculate the vector potential a la ZHS method
//-----------------------------------------------------------------//
void SteppingAction::VPotentialZHSStyle(const G4Step* aStep)
{

  // For now, only count electrons and positrons
  G4Track* track = aStep->GetTrack();
  if( abs(track->GetParticleDefinition()->GetPDGEncoding()) != 11) return;

  // Working in SI Units throughout in order to try to remove
  // any unit conversion errors.
  
  // Specify some constans... Consider moving these to some
  // file that is more easily changed among all scripts
  const G4double m_c  = 2.99792458e8;
  const G4double m_n  = 1.78;      // index of refraction for ice
  const G4double m_mu = 1e-7;      // permeability
  const G4double m_e  = 1.602e-19; // electric charge 
  const G4double m_pi = 3.141592653589;

  // Get the Pre and post step points and times
  G4ThreeVector P0 = aStep->GetPreStepPoint()->GetPosition()    / m;
  G4double      t0 = aStep->GetPreStepPoint()->GetGlobalTime()  / s;
  G4ThreeVector P1 = aStep->GetPostStepPoint()->GetPosition()   / m;
  G4double      t1 = aStep->GetPostStepPoint()->GetGlobalTime() / s;

  if(t1 == t0) return;

  // Now Get midpoint times and position
  //G4double      tm = (t0 + t1)/2.;
  G4ThreeVector Pm = G4ThreeVector( (P0.x() + P1.x())/2.,
				    (P0.y() + P1.y())/2.,
				    (P0.z() + P1.z())/2.);
  
  // Set the velocity vector in m/s
  G4ThreeVector V = getVelocity(P0,P1,t0,t1);
  
  // Loop over each antenna and calculate the vector
  // potential for that particular antenna
  for(unsigned int iA=0; iA<m_ants->size(); ++iA){
    Antenna* ant = m_ants->at(iA);
    
    // Antenna Position
    G4ThreeVector AntPos = G4ThreeVector( ant->getX(),
					  ant->getY(),
					  ant->getZ() );
    
    // Get R and unit vector
    G4double R      = 0;
    G4ThreeVector u = G4ThreeVector();
    setUnitVector(AntPos, Pm, u, R);
    
    // Get the detector time
    G4double tD0 = getTDetector(AntPos, P0, t0, m_n, m_c);				
    G4double tD1 = getTDetector(AntPos, P1, t1, m_n, m_c);				

    //G4acout<<"tD: "<<tD0<<" "<<tD1<<" t: "<<t0<<" "<<t1<<G4endl;

    // Load the antenna timing information
    G4double AntTmin  = ant->getTmin() * 1e-9;  // in seconds
    //G4double AntTmax  = ant->getTmax() * 1e-9;  // in seconds
    G4double AntTstep = ant->getTStep() * 1e-9; // in seconds
    
    // Locate what bin our particle falls in
    G4int iFirstBin = 0;
    G4int iLastBin  = int(ant->getN() - 1);
    G4int iStart    = int( (tD0 - AntTmin)/AntTstep );
    G4int iEnd      = int( (tD1 - AntTmin)/AntTstep );

    // Swap start times if there are issues
    // This really shouldn't happen
    if( iStart > iEnd ){
      G4int temp = iStart;
      iStart = iEnd;
      iEnd = temp;
    }

    // Don't count if out of bounds
    if(iEnd < iFirstBin)  continue;
    if(iStart > iLastBin) continue;
    if(iEnd > iLastBin)    iEnd = iLastBin;
    if(iStart < iFirstBin) iStart = iFirstBin;

    // Debug times
    //G4cout<<"\ttD0: "<<tD0<<" tD1: "<<" bins: "<<iStart<<" "<<iEnd
    //<<" first and last: "<<iFirstBin<<" "<<iLastBin<<G4endl;

    // Define the constant part of vector potential
    G4double charge = m_e * track->GetParticleDefinition()->GetPDGCharge();
    G4double constA = m_mu * charge / (4*m_pi*R);
    
    // We are ready to calculate the vector potential for
    // each vector component
    G4double Ax = constA * ( V.x() - V.dot(u)*u.x() );
    G4double Ay = constA * ( V.y() - V.dot(u)*u.y() );
    G4double Az = constA * ( V.z() - V.dot(u)*u.z() );
    
    // The binning seems to be the most difficult part. There
    // is a 'factor' that is attached to the binning determined
    // by a few limiting cases.  One of the parameters needed is
    // the time derivative of the t_detector
    G4double dtD_dt     = fabs(1. - (m_n/m_c)*(V.dot(u)));
    G4double dtD_dt_min = 1.e-15; // taken from Anne and ZHS code

    G4double factor = 1.;

    // Case 1 -- step doesn't cross bin
    if(iStart == iEnd){
      if( dtD_dt > dtD_dt_min ) factor = fabs((tD1-tD0)/AntTstep/dtD_dt);
      else                      factor = fabs((t1-t0)/AntTstep);
      
      // Save point
      ant->addPoint(iStart, Ax * factor, Ay * factor, Az * factor);
      
    }// end Case 1

    // Case 2 -- step crosses multiple bins
    else{
      
      // Loop over the bins
      for(int ibin = iStart; ibin < iEnd; ++ibin){
	factor = 1./dtD_dt;
	if( tD0 <= tD1 ){
	  if( ibin == iStart ) factor = ((iStart+1)*AntTstep - tD0)/AntTstep/dtD_dt;
	  if( ibin == iEnd )   factor = (tD1 - (iEnd)*AntTstep)/AntTstep/dtD_dt;
	}
	else{
	  if( ibin == iStart ) factor = ((iStart+1)*AntTstep - tD1)/AntTstep/dtD_dt;
	  if( ibin == iEnd )   factor = (tD0 - (iEnd)*AntTstep)/AntTstep/dtD_dt;
	}
	
	// Save the result
	ant->addPoint(ibin, Ax * factor, Ay * factor, Az * factor);
	
      }// end loop over bins

    }// end Case 2      
    
  }// end loop over antennas
  
  
}


//-----------------------------------------------------------------//
// Get the detector time
//-----------------------------------------------------------------//
G4double SteppingAction::getTDetector(G4ThreeVector ant,
				      G4ThreeVector point,
				      G4double time,
				      G4double n_,
				      G4double c_)
{

  //G4cout<<"\t\tGetting time: "<<time<<" "<<(n_/c_)
  //<<" "<<(ant - point).mag()<<" final result: "
  //<<(time + (n_/c_) * (ant - point).mag())<<G4endl;

  return time + (n_/c_) * (ant - point).mag();

}


//-----------------------------------------------------------------//
// Get Velocity vector
//-----------------------------------------------------------------//
G4ThreeVector SteppingAction::getVelocity(G4ThreeVector P0,
					  G4ThreeVector P1,
					  G4double t0,
					  G4double t1)
{

  // Specify dt
  G4double dt = (t1 - t0);
  
  // return the vector
  return G4ThreeVector( (P1.x() - P0.x())/dt,
			(P1.y() - P0.y())/dt,
			(P1.z() - P0.z())/dt
			);
			

}

//-----------------------------------------------------------------//
// Set the unit vector and distance from mid-point to antenna
//-----------------------------------------------------------------//
void SteppingAction::setUnitVector(G4ThreeVector v_ant,
				   G4ThreeVector v_midPoint,
				   G4ThreeVector &u,
				   G4double &R)
{

  // Set the vector from midpoint to antenna
  u = v_ant - v_midPoint;

  // Take R from the vector
  R = u.mag();
  
  // Now divide and make vector unit vector
  u /= R;


}
