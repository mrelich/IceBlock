
#include "SteppingAction.hh"
#include "Constants.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
SteppingAction::SteppingAction(std::ofstream* output,
			       std::vector<Antenna*> *ants,
  			       RefractionTool* refTool) :
  m_output(NULL),
  //m_treeWriter(NULL)
  m_ants(NULL),
  m_refTool(NULL),
  m_TRFirstPointFound(false)
{
  m_output = output;
  //m_treeWriter = treeWriter;
  m_ants = ants;
  m_refTool = refTool;
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
  //EFieldEndpointStyle(aStep);
  //TRFromZHS(aStep);
  
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

  // In order to get the material the pre/post-step points are in
  // one can use either GetPhysicalVolume() or GetMaterial()
  // This can then be used to assign index of refraction.

  // For now, only count electrons and positrons
  G4Track* track = aStep->GetTrack();
  G4int pdg      = track->GetParticleDefinition()->GetPDGEncoding();
  G4int echarge  = track->GetParticleDefinition()->GetPDGCharge(); 
  if( abs(pdg) != 11) return;  

  // Check what medium we are in
  // NOTE: Right now now only considering in ice stuff
  G4double n = m_nAir;
  if( track->GetVolume()->GetName() == "iceblock_phys")
    n = m_n;
  else return; // for now only track in ice particles
  
  // Set the initial and final vertex position
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

    // Save a copy for now. Clean this up later
    // but now want to test refraction quickly.
    G4ThreeVector AntPos_save = AntPos;

    // If the refraction tool has been initialized
    // then use the tool by updating position we are
    // calculating the field at (namely at the boundary)
    G4double theta_i = 0; 
    G4double theta_r = 0;
    G4double dummy   = 0;
    if( m_refTool->isInitialized() ){
      if(m_refTool->useLookup()){
	AntPos = m_refTool->getIntPoint(Pm, iA, theta_i, theta_r, dummy, true, true);	
      }
      else{
	AntPos = m_refTool->getIntPoint(Pm, AntPos, theta_i, theta_r);
	if( theta_i < 0 || theta_r < 0) continue;
      }
    }
    
    // check to make sure antenna position is valid
    if(AntPos.x() <= -9999 || AntPos.y() <= -9999 || AntPos.z() <= -9999)
      continue;

    // Get R and unit vector
    G4double R      = 0;
    G4ThreeVector u = G4ThreeVector();
    setUnitVector(AntPos, Pm, u, R);
    
    // Get the detector time
    G4double tD0 = getTDetector(AntPos, P0, t0, n, m_c);
    G4double tD1 = getTDetector(AntPos, P1, t1, n, m_c);

    // Now add back the propagation time in air
    if( m_refTool->isInitialized() ){
      tD0 = getTDetector(AntPos_save, AntPos, tD0, m_nAir, m_c);
      tD1 = getTDetector(AntPos_save, AntPos, tD1, m_nAir, m_c);
      //tD0 = getTDetector(AntPos_save, AntPos, tD0, m_n, m_c);
      //tD1 = getTDetector(AntPos_save, AntPos, tD1, m_n, m_c);
    }

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
    G4double R_pa = (AntPos_save - AntPos).mag();
    G4double charge = m_e * echarge; //track->GetParticleDefinition()->GetPDGCharge();
    G4double constA = m_mu * charge / (4*m_pi*(R+R_pa));

    // We are ready to calculate the vector potential for
    // each vector component
    G4double Ax = constA * ( V.x() - V.dot(u)*u.x() );
    G4double Ay = constA * ( V.y() - V.dot(u)*u.y() );
    G4double Az = constA * ( V.z() - V.dot(u)*u.z() );

    // Now get the refracted field using fresnel's coefficients
    // And also scale 1/R to account for distance from iceblock 
    // surface to antenna position
    if( m_refTool->isInitialized()){

      //G4double R_pa = (AntPos_save - AntPos).mag();
      //if(R_pa != 0){
      //	Ax /= R_pa;
      //Ay /= R_pa;
      //Az /= R_pa;
      //}

      if(m_refTool->useTool()){
	G4ThreeVector Aref = m_refTool->getTransmittedField(G4ThreeVector(Ax,Ay,Az),theta_i);
	Ax = Aref.x();
	Ay = Aref.y();
	Az = Aref.z();
      }
    }

    if( m_output ){
      G4ThreeVector vk = G4ThreeVector(0,0,1);
      (*m_output) << "\t" << iA
		  <<" "<< Az
		  <<" "<< V.angle((AntPos_save - Pm))
		  <<" "<< V.angle((AntPos - Pm))
		  <<" "<< Pm.x() << " " << Pm.y() << " " << Pm.z()
		  <<" "<< AntPos.x() << " " << AntPos.y() << " " << AntPos.z()
		  <<" "<< AntPos_save.x() << " " << AntPos_save.y() << " " << AntPos_save.z()
		  <<" "<< vk.angle(AntPos.unit())
		  <<" "<< vk.angle(AntPos_save.unit())
		  <<" "<< V.x() << " " << V.y() << " " << V.z()
		  <<G4endl;
    }

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
    G4double dtD_dt     = fabs(1. - (n/m_c)*(V.dot(u)));
    G4double dtD_dt_min = m_tolerance; //1.e-15; // taken from Anne and ZHS code
    
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
// Set intial and final points
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
  // <<(time + (n_/c_) * (ant - point).mag())<<G4endl;

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
				   G4ThreeVector v_Point,
				   G4ThreeVector &u,
				   G4double &R)
{

  // Set the vector from midpoint to antenna
  u = v_ant - v_Point;

  // Take R from the vector
  R = u.mag();
  
  // Now divide and make vector unit vector
  u /= R;


}

//-----------------------------------------------------------------//
// Get the Efield from Endpoint method
//-----------------------------------------------------------------//
void SteppingAction::EFieldEndpointStyle(const G4Step* aStep)
{

  // Work in CGS units
  G4double c = m_c / 1e-2;
  G4double q = m_e / 1e-1;
  G4double Econv = 1/1.e6; //1V/m = 10^6 abV/cm

  // For convenince later
  G4StepPoint* prestep  = aStep->GetPreStepPoint();
  G4StepPoint* poststep = aStep->GetPostStepPoint();
  G4Track* track = aStep->GetTrack();

  // only follow electrons
  if( abs(track->GetParticleDefinition()->GetPDGEncoding()) != 11) return;
  //G4cout<<"Track: "<<track->GetTrackLength()/cm<<" "<<track->GetTrackID()<<G4endl;

  // Alternative calculation for velocity
  // This will use average velocity between steps. Much better
  // than taking pre and post step position and time where we 
  // can have speeds greater than light in vacuum due to neglecting
  // acceleration.
  G4ThreeVector P0   = prestep->GetPosition()    / cm;
  G4double      t0   = prestep->GetGlobalTime()  / s;

  G4double      Vi   = (prestep->GetVelocity() +
			poststep->GetVelocity())/ (cm/s);
  Vi = Vi/2.;
  G4ThreeVector post = poststep->GetPosition()    / cm;
  G4ThreeVector uv   = (post-P0)/(post-P0).mag();

  G4double      t1   = poststep->GetGlobalTime()  / s;
  G4ThreeVector P1   = P0 + uv * Vi * (t1-t0);


  //G4cout<<"Pre:  "<<aStep->GetPreStepPoint()->GetVelocity()<<G4endl;
  //G4cout<<"Post: "<<aStep->GetPostStepPoint()->GetVelocity()<<G4endl;

  // If track length or time are same do not use this track
  if( (P0-P1).mag() == 0 ) return;
  if(t1 == t0)             return;

  // Set Beta
  G4ThreeVector Beta = getVelocity(P0,P1,t0,t1) / c;

  // Loop over the antennas and calculate the field 
  // recieved at each antenna.
  for(unsigned int iA=0; iA<m_ants->size(); ++iA){
    Antenna* ant = m_ants->at(iA);
    
    // Antenna Position -- put into cm
    G4ThreeVector AntPos = G4ThreeVector( ant->getX()/1e-2,
					  ant->getY()/1e-2,
					  ant->getZ()/1e-2 );
    
    // Antenna timing
    // Remember that for e-field the antenna 
    // is 0.5*step forward in time
    G4double AntTstep = ant->getTStep() * 1e-9;               // in seconds
    G4double AntTmin  = ant->getTmin() * 1e-9;                // in seconds
    //G4double AntTmin  = (ant->getTmin()+AntTstep/2.) * 1e-9;  // in seconds


    // Get Rs and unit vectors
    G4double R0 = 0, R1 = 0;
    G4ThreeVector u0 = G4ThreeVector(), u1 = G4ThreeVector();
    setUnitVector(AntPos, P0, u0, R0);
    setUnitVector(AntPos, P1, u1, R1);

    // Get the times at antenna
    G4double tD0 = getTDetector(AntPos,P0,t0,m_n,c);
    G4double tD1 = getTDetector(AntPos,P1,t1,m_n,c);
    G4double dt  = fabs(t0-t1); //fabs(tD0-tD1);

    // Get the charge
    G4double charge = q * track->GetParticleDefinition()->GetPDGCharge();

    // Locate what bin our particle falls into
    G4int ibin0       = int( (tD0 - AntTmin)/AntTstep );
    G4int ibin1       = int( (tD1 - AntTmin)/AntTstep );
    
    // Get the Efields
    //G4ThreeVector E0 = getEFieldEndpoint(Beta,u0,R0,dt,charge);
    //G4ThreeVector E1 = -1*getEFieldEndpoint(Beta,u1,R1,dt,charge);
    G4ThreeVector E0 = getEFieldEndpoint(Beta,u0,R0,AntTstep,charge,c) * Econv;
    G4ThreeVector E1 = -1*getEFieldEndpoint(Beta,u1,R1,AntTstep,charge,c) * Econv;

    if(false){
      G4cout<<"-----------------------------------------------"<<G4endl;
      G4cout<<"TrkID: "<<track->GetTrackID()<<G4endl;
      G4cout<<ibin0<<" "<<tD0<<" "<<E0<<G4endl;
      G4cout<<ibin1<<" "<<tD1<<" "<<E1<<G4endl;
      G4cout<<"Positions: "<<P0<<" "<<P1<<G4endl;
      G4cout<<"times: "<<t0<<" "<<t1<<" "<<t1-t0<<G4endl;
      G4cout<<"Delta: "<<(tD1-tD0)<<" "<<E0+E1<<G4endl;
      G4cout<<"Beta:  "<<Beta<<G4endl;
      G4cout<<"Distances: "<<R0<<" "<<R1<<" "<<u0<<" "<<u1<<G4endl;
    }

    // Store the information here
    // I think this should only work if we limit the step size
    // Step size limited in DetectorConstruction to something [mm]
    // Does this work? or does the binning need to be more complicated
    // like it was for ZHS??
    
    // Tried treating first stepping point as stopping point.
    // Doesn't change the result at all.
    //if(track->GetTrackID() == 1 && prestep->GetGlobalTime() ==0){
    //ant->addEPoint(ibin0, -E0.x(), -E0.y(), -E0.z());
    //}
    //else 
    ant->addEPoint(ibin0, E0.x(), E0.y(), E0.z());
    ant->addEPoint(ibin1, E1.x(), E1.y(), E1.z());
    

    //G4cout<<"Is first: "<<isFirstPoint<<G4endl;
    //G4cout<<"bin:      "<<ibin<<" time: "<<t<<" ret t: "<<tD<<G4endl;
    //G4cout<<Efield.z()<<G4endl;
    //G4cout<<"Ant bin:  "<<AntTstep<<G4endl;
    //G4cout<<"Ant Tmin: "<<AntTmin<<" "<<tD-AntTmin<<G4endl;

    if(false && dt > AntTstep){
      G4cout<<"Somehow have larger step? "<<dt<<" "<<AntTstep<<G4endl;
    }
		   
  }// end loop over antennas

}

//-----------------------------------------------------------------//
// Calculate the efield using equation 8 from end-point
//-----------------------------------------------------------------//
G4ThreeVector SteppingAction::getEFieldEndpoint(G4ThreeVector Beta,
						G4ThreeVector rhat,
						G4double R,
						G4double dt,
						G4double q,
						G4double c)
{

  // Get easy constant terms
  G4double C = q/(dt*R*c);
  
  // Get the numerator information
  G4ThreeVector num = rhat.cross(rhat.cross(Beta));
  
  // Get the denominator information
  G4double den = 1-m_n*Beta.dot(rhat);
  
  // Add some protections
  if( fabs(den) < m_tolerance ){
    G4cout<<"Below tolerance: "<<den<<G4endl;
    den = m_tolerance;
  }

  
  // Dump some information
  if(false){
    G4cout<<"---------------------------------------------"<<G4endl;
    G4cout<<"dt:   "<<dt<<G4endl;
    G4cout<<"Rhat: "<<rhat<<G4endl;
    G4cout<<"q:    "<<q<<G4endl;
    G4cout<<"R:    "<<R<<G4endl;
    G4cout<<"Beta: "<<Beta<<" "<<Beta.mag()<<G4endl;
    G4cout<<"Num:  "<<num<<G4endl;
    G4cout<<"den:  "<<den<<G4endl;
    G4cout<<"C:    "<<C<<G4endl;
    G4cout<<"Tot:  "<<C * num / den<<G4endl;
  }

  return C * num / den;
  
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

  if( !m_refTool->isInitialized() ){
    G4cout<<"In order to run TRFromZHSa the refraction tool"<<G4endl;
    G4cout<<"must be initialized and use the lookup tables"<<G4endl;
    G4cout<<"Doing nothing..."<<G4endl;
    return;
  }

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
  if( !(mat1->GetName() == "ICE" && (mat0 != mat1 || m_TRFirstPointFound))) return;
  
  // Temporary
  G4Track* trk = aStep->GetTrack();
  if( trk->GetParentID() != 0 ) return;

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
  G4double n   = m_nAir;
  G4bool inice = false;
  G4double shift = 0;
  if( !m_TRFirstPointFound );
  else{
    n = m_n;
    inice = true;
    shift = 0.01; // hack to handle the tables
  }
  m_TRFirstPointFound = !m_TRFirstPointFound;

  // For the in ice part, we need to also calculate
  // the reflected path.  Do this by looping over a vector
  // For now skipping reflected part...
  G4bool v_direct[2] = {true, false};
  int nDirectReflect = inice ? 1 : 1;

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
  
  // Now Get midpoint times and position
  G4ThreeVector Pm = G4ThreeVector( (P0.x() + P1.x())/2.,
				    (P0.y() + P1.y())/2.,
				    (P0.z() + P1.z())/2. + shift);
  
  // Loop over each antenna and calculate the vector
  // potential for that particular antenna
  for(unsigned int iA=0; iA<m_ants->size(); ++iA){
    Antenna* ant = m_ants->at(iA);
    
    // Get the vector that defines the direction we want
    // to calculate the electric field. In a homogenous medium,
    // this would be the antenna position. When dealing with refraction
    // we take a point on the ice.
    G4ThreeVector AntPos;
    G4double theta_i = 0; 
    G4double theta_r = 0;
    G4double R       = 0;

    // Loop over the direct and reflected based on whether or
    // not we are in the ice.
    for(int id = 0; id<nDirectReflect; ++id){
      
      // Update the antenna position
      AntPos = m_refTool->getIntPoint(Pm, iA, theta_i, theta_r, R, inice, v_direct[id]);	
      
      // check to make sure antenna position is valid
      if(AntPos.x() <= -9999 || AntPos.y() <= -9999 || AntPos.z() <= -9999){
	continue;
      }

      // Set the unit vector
      G4double dummy  = 0;
      G4ThreeVector u = G4ThreeVector();
      setUnitVector(AntPos, Pm, u, dummy);
      
      // Get the detector time 
      // NOTE: This is done differently here. Romain's lookup 
      // tables contain the optical path, which already accounts
      // for index of refraction. so t_detector = t_emission + R/c
      G4double tD0 = t0 + R/m_c;
      G4double tD1 = t1 + R/m_c;
      
      // Define the constant part of vector potential
      G4double charge = m_e * echarge; //track->GetParticleDefinition()->GetPDGCharge();
      G4double constA = m_mu * charge / (4*m_pi*R);
      
      // We are ready to calculate the vector potential for
      // each vector component
      G4ThreeVector A = G4ThreeVector(constA * ( V.x() - V.dot(u)*u.x() ),
				      constA * ( V.y() - V.dot(u)*u.y() ),
				      constA * ( V.z() - V.dot(u)*u.z() ));
      
      // Now get the refracted field using fresnel's coefficients.
      // Do this twice if we are going from Air->Ice->Air
      if(m_refTool->useTool() && !inice){
	//G4cout<<"Not in ice: "<<A<<G4endl;
	A = m_refTool->getTransmittedField(A,theta_i,false);
	//G4cout<<"\t"<<A<<G4endl;
	A = m_refTool->getTransmittedField(A,theta_r,true);
	//G4cout<<"\t"<<A<<G4endl;
      }
      if(m_refTool->useTool() && inice){
	//G4cout<<"InIce: "<<A<<G4endl;
	A = m_refTool->getTransmittedField(A,theta_i);
	//G4cout<<"\t"<<A<<G4endl;
      }

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
