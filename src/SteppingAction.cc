
#include "SteppingAction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingAction::SteppingAction(std::ofstream* output,
			       std::vector<Antenna*> *ants) :
  m_output(NULL),
  //m_treeWriter(NULL)
  m_ants(NULL),
  m_TRFirstPointFound(false),
  m_trtool(NULL)
{
  m_output = output;
  //m_treeWriter = treeWriter;
  m_ants = ants;
  //m_trtool = new TRTool(m_nAir, m_n);
  m_trtool = new TRTool(m_n, m_nAir);
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
  //VPotentialZHSStyle(aStep);
  //VPotentialEndpoint(aStep);
  TRFromZHS(aStep);

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

  if( m_output ){
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
  }// end if file is open

}

//-----------------------------------------------------------------//
// Calculate the vector potential a la ZHS method
//-----------------------------------------------------------------//
void SteppingAction::VPotentialZHSStyle(const G4Step* aStep)
{

  // For now, only count electrons and positrons
  G4Track* track = aStep->GetTrack();
  if( abs(track->GetParticleDefinition()->GetPDGEncoding()) != 11) return;

  // Placing a 1 MeV Energy threhold
  //if( aStep->GetPostStepPoint()->GetTotalEnergy()/MeV < 1 ) return;

  // Working in SI Units throughout in order to try to remove
  // any unit conversion errors.
  
  // Get the Pre and post step points and times
  G4ThreeVector P0;
  G4ThreeVector P1;
  G4double t0=0, t1=0;
  setInitialFinalPoint(aStep, P0, P1, t0, t1);


  // If track length or time are same do not use this track
  if( (P0-P1).mag() == 0 ) return;
  if(t1 == t0)             return;
  
  // Set the velocity vector in m/s
  G4ThreeVector V = getVelocity(P0,P1,t0,t1);
  //G4ThreeVector V = ((aStep->GetPostStepPoint()->GetVelocity()*s/m + 
  //aStep->GetPreStepPoint()->GetVelocity()*s/m) / 2.) * (P1-P0)/(P1-P0).mag();    
  
  // Now Get midpoint times and position
  G4ThreeVector Pm = G4ThreeVector( (P0.x() + P1.x())/2.,
				    (P0.y() + P1.y())/2.,
				    (P0.z() + P1.z())/2.);


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
    //if(iEnd > iLastBin)    continue;
    //if(iStart < iFirstBin) continue;

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

    //G4cout<<"A: "<<Ax<<" "<<Ay<<" "<<Az<<G4endl;
    //G4cout<<"v: "<<V<<G4endl;
    //G4cout<<"dot: "<<V.dot(u)*u.x()
    //	  <<" "<<V.dot(u)*u.y()
    //	  <<" "<<V.dot(u)*u.z()
    //	  <<G4endl;

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
      for(int ibin = iStart; ibin <= iEnd; ++ibin){
	factor = 1./dtD_dt;
	if( tD0 <= tD1 ){
	  if( ibin == iStart ) factor = ((AntTmin+(iStart+1)*AntTstep) - tD0)/AntTstep/dtD_dt;
	  if( ibin == iEnd )   factor = (tD1 - (AntTmin+(iEnd)*AntTstep))/AntTstep/dtD_dt;
	}
	else{
	  if( ibin == iStart ) factor = ((AntTmin+(iStart+1)*AntTstep) - tD1)/AntTstep/dtD_dt;
	  if( ibin == iEnd )   factor = (tD0 - (AntTmin+(iEnd)*AntTstep))/AntTstep/dtD_dt;
	}
	
	if(factor < 0) 
	  G4cout<<"Factor: "<<factor
	      <<" ibin: "<<ibin
	      <<" start: "<<iStart
	      <<" end: "<<iEnd<<G4endl;
	
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
// Get Beta ZHS style
//-----------------------------------------------------------------//
G4double SteppingAction::getBeta(const G4Step* step)
{

  G4double E0 = step->GetPreStepPoint()->GetTotalEnergy() / MeV;
  G4double E1 = step->GetPostStepPoint()->GetTotalEnergy() / MeV;
  
  G4double Eavg = (E1+E0)/2.;
  
  return sqrt(1 - 1/pow(Eavg/0.511,2));

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

//-----------------------------------------------------------------//
// Set the initial and final point
//-----------------------------------------------------------------//
void SteppingAction::setInitialFinalPoint(const G4Step* step,
                                          G4ThreeVector &P0,
                                          G4ThreeVector &P1,
                                          G4double &t0,
                                          G4double &t1)
{

  // Get the Pre and post step points and times
  // This is taken from Anne's code, and seems to give
  // velocities greater than light in vacuum for some steps.
  /*
  P0 = aStep->GetPreStepPoint()->GetPosition()    / m;
  P1 = aStep->GetPostStepPoint()->GetPosition()   / m;
  t0 = aStep->GetPreStepPoint()->GetGlobalTime()  / s;
  t1 = aStep->GetPostStepPoint()->GetGlobalTime() / s;
  */

  // Alternative calculation
  // This will use average velocity between steps. Much better
  // than taking pre and post step position and time where we 
  // can have speeds greater than light in vacuum due to neglecting
  // acceleration.
  P0 = step->GetPreStepPoint()->GetPosition()    / m;
  t0 = step->GetPreStepPoint()->GetGlobalTime()  / s;

  G4double Vi = (step->GetPreStepPoint()->GetVelocity() +
		 step->GetPostStepPoint()->GetVelocity())/ (m/s);
  Vi = Vi/2.;

  G4ThreeVector post = step->GetPostStepPoint()->GetPosition()    / m;
  G4ThreeVector uv   = (post-P0)/(post-P0).mag();

  t1 = step->GetPostStepPoint()->GetGlobalTime()  / s;
  P1   = P0 + uv * Vi * (t1-t0);
  
}


//-----------------------------------------------------------------//
// Transition radiation component
// This is an idea we had inspired from end-point method where we 
// calculate the 'transition radiation' compenent from the point
// just crossing the boundary.  As it turns out, Geant4 actually 
// stops the particle *right* before it crosses the boundary of
// a material.  So we will have take the material transition point,
// which occurs first, as our step outside the material, and then
// the next step will be considered the first step inside the
// material.  Luckily Geant4 walks through all steps for a single
// particle at a time.  So we can simply update a bool to keep track
// of the first and second steps for this calculation.
//-----------------------------------------------------------------//
void SteppingAction::TRFromZHS(const G4Step* aStep)
{

  // Temporary
  G4Track* trk = aStep->GetTrack();
  if( trk->GetParentID() != 0 ) return;

  // Get the previous and post step points. 
  G4StepPoint* SP0 = aStep->GetPreStepPoint();
  G4StepPoint* SP1 = aStep->GetPostStepPoint();

  // Right now we are only considering the transition from 
  // air to ice.  So only keep points where that step is
  // happening.
  G4Material* mat0 = SP0->GetMaterial();
  G4Material* mat1 = SP1->GetMaterial();

  // If we get to the boundary (meaning last point leaves the
  // world volume), then the material will not be defined. We
  // aren't interested in this point, so discard if either 
  // material is undefined
  if(!mat0 || !mat1) return;
  //if( !(mat1->GetName() == "ICE" && (mat0 != mat1 || m_TRFirstPointFound))) return;
  if( !(mat1->GetName() == "G4_AIR" && (mat0 != mat1 || m_TRFirstPointFound))) return;
  G4cout<<mat0->GetName()<<" "<<mat1->GetName()<<G4endl;  

  // There is also a contribution from backscatter electrons
  // meaning electrons produced in air above ice that scatter
  // back into the ice.  The contribution is small, but the code
  // is not designed to handle these guys.  This check is only
  // necessary for the first point
  // TODO: Add more thorough treatment to correctly remove
  // events that scatter back into the ice.
  G4ThreeVector P0 = SP0->GetPosition() / m;
  G4ThreeVector P1 = SP1->GetPosition() / m;
  if( !m_TRFirstPointFound ){
    if( !(P0.z() < P1.z() && P0.z() < 0 ) ) return;
  }

  // Now print out some info to test ice locations and everything
  G4cout<<"-------------------------------------------------------"<<G4endl;
  G4cout<<"P0: "<<P0<<G4endl;
  G4cout<<"P1: "<<P1<<G4endl;
  G4cout<<"Mat0: "<<mat0->GetName()<<G4endl;
  G4cout<<"Mat1: "<<mat1->GetName()<<G4endl;
  G4cout<<G4endl;
  
  // Now go through the calculation for this point. Two cases:
  //    1.) This is the first point -- we just get contribution
  //        from the ray passing through the Ice obeying Snell.
  //    2.) This is the second point -- we get ray from inside Ice
  //        and also the reflected one.
  
  // Now update whether or not we are working with the first point.
  // If we are with the first point, then TRFirstPoitnFound
  // should be false.  If it is the second point then it's true.
  // Update accordingly.
  // Check what medium we are in
  //G4double n   = m_nAir;
  //G4bool inice = false;
  //if( !m_TRFirstPointFound );
  //else{
  //  n = m_n;
  //  inice = true;
  // }
  //m_TRFirstPointFound = !m_TRFirstPointFound;

  G4double n   = m_n;
  G4bool inice = true;
  if( !m_TRFirstPointFound );
  else{
    n = m_nAir;
    inice = false;
  }
  m_TRFirstPointFound = !m_TRFirstPointFound;

  G4cout<<"index: "<<n<<G4endl;

  // For now, only count electrons and positrons
  G4Track* track = aStep->GetTrack();
  G4int pdg      = track->GetParticleDefinition()->GetPDGEncoding();
  G4int echarge  = track->GetParticleDefinition()->GetPDGCharge(); 
  if( abs(pdg) != 11) return;  

  // Set the initial and final vertex position
  G4double t0=0, t1=0;
  setInitialFinalPoint(aStep, P0, P1, t0, t1);

  // If track length or time are same do not use this track
  if( (P0-P1).mag() == 0 ) return;
  if(t1 == t0)             return;

  // Set the velocity vector in m/s
  G4ThreeVector V = getVelocity(P0,P1,t0,t1);

  // HACK for testing. The idea was to fix the single step
  // to +/- X cm and compare shape to ZHS paper.  This becomes
  // complicated though, as you need to swap the geometry 
  // and many other constants
  if(false){
    G4cout<<"Before "<<t0<<" "<<t1<<G4endl;
    V  =  0.99* m_c * G4ThreeVector(0,0,1);
    if( mat0->GetName() == "ICE" ){
      P0 = G4ThreeVector(0,0,-0.05);
      P1 = G4ThreeVector(0,0,0);
      t0 = 0; //t0 - (P1-P0).mag() / V.mag();
      t1 = t0 + (P1-P0).mag() / V.mag();
    }
    else{
      P0 = G4ThreeVector(0,0,0);
      P1 = G4ThreeVector(0,0,0.05);
      t0 = (P1-P0).mag() / V.mag();
      t1 = t0 + (P1-P0).mag() / V.mag();
    }
    G4cout<<"After "<<t0<<" "<<t1<<G4endl;
  }
  
  // Now Get midpoint times and position
  G4ThreeVector Pm = G4ThreeVector( (P0.x() + P1.x())/2.,
                                    (P0.y() + P1.y())/2.,
                                    (P0.z() + P1.z())/2.);

  // potential for that particular antenna
  for(unsigned int iA=0; iA<m_ants->size(); ++iA){
    Antenna* ant = m_ants->at(iA);
    
    // Get the vector that defines the direction we want
    // to calculate the electric field. In a homogenous medium,
    // this would be the antenna position. When dealing with refraction
    // we take a point on the ice.
    G4ThreeVector AntPos = G4ThreeVector(ant->getX(),
					 ant->getY(),
					 ant->getZ());
    G4double theta_i = 0; 
    G4double theta_r = 0;
    G4double R       = 0;
    G4double tD0     = 0;
    G4double tD1     = 0;
    G4ThreeVector u = G4ThreeVector();

    // Loop over the direct and reflected based on whether or
    // not we are in the ice.
    for(int it = 0; it<TRR_N; ++it){
      TRRay tr_it = (TRRay) it;

      //if(tr_it != TRR_refracted) continue;

      // Handle cases
      //if( inice && tr_it == TRR_refracted ) continue;
      //if( !inice && tr_it != TRR_refracted ) continue;
      if( inice && tr_it != TRR_refracted ) continue;
      if( !inice && tr_it == TRR_refracted ) continue;

      // Set the variables
      m_trtool->findPath(Pm, AntPos, t0, t1, theta_i, theta_r, R, u, tD0, tD1, tr_it);

      // Need to make sure distance makes sense
      if( R < 0 ){
	G4cout<<"Some issue probably in TRTool as distance is < 0"<<G4endl;
	G4cout<<"Distance: "<<R<<G4endl;
	continue;
      }

      // Define the constant part of vector potential
      G4double charge = m_e * echarge; //track->GetParticleDefinition()->GetPDGCharge();
      G4double constA = m_mu * charge / (4*m_pi*R);
      
      // We are ready to calculate the vector potential for
      // each vector component
      G4ThreeVector A = G4ThreeVector(constA * ( V.x() - V.dot(u)*u.x() ),
                                      constA * ( V.y() - V.dot(u)*u.y() ),
                                      constA * ( V.z() - V.dot(u)*u.z() ));

      // TODO: Add the transmission and reflection coeff to TRTool and
      // implement that here!!!!
      A = m_trtool->getFresnelCorrectedField(A, theta_i, theta_r, tr_it);

      // Get the time derivative which is needed for binning
      G4double dtD_dt     = fabs(1. - (n/m_c)*(V.dot(u)));
            
      fillForAntenna(ant, t0, t1, tD0, tD1, dtD_dt, A);
    
    }// end loop over antennas
  }// end loop over direct and reflected options
  
}

//-----------------------------------------------------------------//
// Method to fill the antenna for the ZHS TR calculation.
//-----------------------------------------------------------------//
void SteppingAction::fillForAntenna(Antenna* ant,
                                    G4double t0,
                                    G4double t1,
                                    G4double tD0,
                                    G4double tD1,
                                    G4double dtD_dt,
                                    G4ThreeVector A)
{

  // Load the antenna timing information
  G4double AntTmin  = ant->getTmin() * 1e-9;  // in seconds
  G4double AntTstep = ant->getTStep() * 1e-9; // in seconds

  //G4cout<<"AntTmin: "<<AntTmin<<" "<<(AntTmin + 4000*AntTstep)<<G4endl;
  
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
  if(iEnd < iFirstBin)  return;
  if(iStart > iLastBin) return;
  if(iEnd > iLastBin)    iEnd = iLastBin;
  if(iStart < iFirstBin) iStart = iFirstBin;
  //if(iEnd > iLastBin)    continue;
  //if(iStart < iFirstBin) continue;
  
  // The binning seems to be the most difficult part. There
  // is a 'factor' that is attached to the binning determined
  // by a few limiting cases.  One of the parameters needed is
  // the time derivative of the t_detector
  G4double dtD_dt_min = m_tolerance; //1.e-15; // taken from Anne and ZHS code
  
  G4double factor = 1.;

  // Case 1 -- step doesn't cross bin
  if(iStart == iEnd){
    if( dtD_dt > dtD_dt_min ) factor = fabs((tD1-tD0)/AntTstep/dtD_dt);
    else                      factor = fabs((t1-t0)/AntTstep);
    
    // Save point
    ant->addPoint(iStart, A.x() * factor, A.y() * factor, A.z() * factor);
    
  }// end Case 1
  
  // Case 2 -- step crosses multiple bins
  else{
    // Loop over the bins
    for(int ibin = iStart; ibin <= iEnd; ++ibin){
      factor = 1./dtD_dt;
      if( tD0 <= tD1 ){
        if( ibin == iStart ) factor = ((AntTmin+(iStart+1)*AntTstep) - tD0)/AntTstep/dtD_dt;
        if( ibin == iEnd )   factor = (tD1 - (AntTmin+(iEnd)*AntTstep))/AntTstep/dtD_dt;
      }
      else{
        if( ibin == iStart ) factor = ((AntTmin+(iStart+1)*AntTstep) - tD1)/AntTstep/dtD_dt;
        if( ibin == iEnd )   factor = (tD0 - (AntTmin+(iEnd)*AntTstep))/AntTstep/dtD_dt;
      }
      
      if(factor < 0) 
        G4cout<<"Factor: "<<factor
              <<" ibin: "<<ibin
              <<" start: "<<iStart
              <<" end: "<<iEnd<<G4endl;
      
      // Save the result
      ant->addPoint(ibin, A.x() * factor, A.y() * factor, A.z() * factor);
      
    }// end loop over bins
    
  }// end Case 2      
    
}

