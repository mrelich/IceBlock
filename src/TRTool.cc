
#include "TRTool.hh"

//-----------------------------------------------------------//
// Constructor
//-----------------------------------------------------------//
TRTool::TRTool(G4double n0, G4double n1) :
  m_n0(n0),
  m_n1(n1)
{

}

//-----------------------------------------------------------//
// Destructor
//-----------------------------------------------------------//
TRTool::~TRTool()
{
  //empty
}

//-----------------------------------------------------------//
// Find the path and angles for a given point
//-----------------------------------------------------------//
void TRTool::findPath(G4V3 Pm, 
		      G4V3 Pa, 
		      G4double t0,
		      G4double t1,
		      G4double &theta_i,
		      G4double &theta_r, 
		      G4double &distance,
		      G4V3 &u_r,
		      G4double &t0_ret,
		      G4double &t1_ret,
		      TRRay trr_type)
{
  
  if( trr_type == TRR_direct) 
    findPathDirect(Pm, Pa, t0, t1, theta_i, theta_r, distance, u_r, t0_ret, t1_ret);

  else if( trr_type == TRR_reflected) 
    findPathReflected(Pm, Pa, t0, t1, theta_i, theta_r, distance, u_r, t0_ret, t1_ret);

  else if( trr_type == TRR_refracted) 
    findPathRefracted(Pm, Pa, t0, t1, theta_i, theta_r, distance, u_r, t0_ret, t1_ret);
  
  else{
    G4cout<<"There is something wrong. You are trying to get"<<G4endl;
    G4cout<<"an unsupported type for findPath in TRTool"<<G4endl;
    G4cout<<"Returning garbage"<<G4endl;
    theta_i  = -999;
    theta_r  = -999;
    distance = -999;
    u_r.set(-999,-999,-999);
    t0_ret = -999;
    t1_ret = -999;
  }

}

//-----------------------------------------------------------//
// Direct path inside of ice to antenna
//-----------------------------------------------------------//
void TRTool::findPathDirect(G4V3 Pm, G4V3 Pa,
			    G4double t0, G4double t1,
			    G4double &theta_i, G4double &theta_r,
			    G4double &distance, G4V3 &u_r,
			    G4double &t0_ret, G4double &t1_ret)
{

  // This case is the standard and is rather trivial
  theta_i  = -999;
  theta_r  = -999;
  distance = (Pa - Pm).mag();
  u_r      = (Pa - Pm) / distance;
  
  t0_ret = t0 + m_n1 * distance / m_c;
  t1_ret = t1 + m_n1 * distance / m_c;     

  if(false){
    G4cout<<G4endl;
    G4cout<<"-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-="<<G4endl;
    G4cout<<"Final Direct results: "<<G4endl;
    G4cout<<"Distance: "<<distance<<" "<<u_r<<G4endl;
    G4cout<<"Times: "<<t0<<" "<<t1<<G4endl;
    G4cout<<"Ret:   "<<t0_ret<<" "<<t1_ret<<G4endl;
    G4cout<<G4endl;
  }

}

//-----------------------------------------------------------//
// Reflected path inside ice
// Only working in 2 dimensions. Assuming that antennas are
// in the x-z plane.
//-----------------------------------------------------------//
void TRTool::findPathReflected(G4V3 Pm, G4V3 Pa,
			       G4double t0, G4double t1,
			       G4double &theta_i, G4double &theta_r,
			       G4double &distance, G4V3 &u_r,
			       G4double &t0_ret, G4double &t1_ret)
{
  
  // x Point will lie between Pa and Pm. 
  // Enforce this
  G4double Pi_x = (Pm.z()*Pa.x() + Pm.x()*Pa.z()) / (Pm.z() + Pa.z());
  if( (Pi_x > Pm.x() && Pi_x > Pa.x()) || (Pi_x < Pm.x() && Pi_x < Pa.x()) ){
    G4cout<<"Something wrong... reflected point not in bounds"<<G4endl;
    G4cout<<"Pm: "<<Pm<<G4endl;
    G4cout<<"Pa: "<<Pa<<G4endl;
    G4cout<<"Pi_x: "<<Pi_x<<G4endl;
  }

  // Here assuming only 2-D case, so we don't solve for y
  G4ThreeVector Pi = G4ThreeVector(Pi_x, 0, 0);

  // Distance and unit vector
  u_r = (Pi - Pm) / (Pi - Pm).mag();
  distance = (Pa - Pi).mag() + (Pi - Pm).mag();
  
  // Angles
  theta_i = m_pi/2. - atan(Pm.z()/fabs(Pm.x() - Pi.x()));
  theta_r = -999;
  
  // Timing information can be done as in the direct case
  // since we are contained in one medium
  t0_ret = t0 + m_n1 * distance / m_c;
  t1_ret = t1 + m_n1 * distance / m_c;

  if(false){
    G4cout<<G4endl;
    G4cout<<"-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-="<<G4endl;
    G4cout<<"Final Reflected results: "<<G4endl;
    G4cout<<"Distance: "<<distance<<" "<<u_r<<" theta_i: "<<theta_i*180/m_pi<<G4endl;
    G4cout<<"Ret:   "<<t0_ret<<" "<<t1_ret<<G4endl;
    G4cout<<"Pm: "<<Pm<<G4endl;
    G4cout<<"Pa: "<<Pa<<G4endl;
    G4cout<<"Pi_x: "<<Pi_x<<G4endl;
    G4cout<<G4endl;
  }

}

//-----------------------------------------------------------//
// Refracted path from air --> ice
//-----------------------------------------------------------//
void TRTool::findPathRefracted(G4V3 Pm, G4V3 Pa,
			       G4double t0, G4double t1,
			       G4double &theta_i, G4double &theta_r,
			       G4double &distance, G4V3 &u_r,
			       G4double &t0_ret, G4double &t1_ret)
{

  // Find the interaction point
  G4V3 Pi = findIntPoint(Pm, Pa);
  
  // Now we have interaction point, calculate the two distances
  G4double dm = (Pi - Pm).mag();
  G4double da = (Pa - Pi).mag();
  
  // set the unit vector for the point outside the ice
  u_r = (Pi - Pm) / dm;
  
  // Set the distance
  distance = dm + da;

  // Set the angles
  theta_i = atan( fabs(Pm.x() - Pi.x()) / fabs(Pm.z()) ); 
  theta_r = atan( fabs(Pa.x() - Pi.x()) / fabs(Pa.z()) ); 
  
  // Set the times
  t0_ret = t0 + m_n0 * dm / m_c + m_n1 * da / m_c;
  t1_ret = t1 + m_n0 * dm / m_c + m_n1 * da / m_c;

  if( false ){
    G4cout<<G4endl;
    G4cout<<"-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-="<<G4endl;
    G4cout<<"Final Refracted results: "<<G4endl;
    G4cout<<"Angles: "<<theta_i*180/m_pi<<" "<<theta_r*180/m_pi<<G4endl;
    G4cout<<"Snell ok: "<<m_n0*sin(theta_i)<<" "<<m_n1*sin(theta_r)<<G4endl;
    G4cout<<"Ret:   "<<t0_ret<<" "<<t1_ret<<" "<<t1_ret-t0_ret<<G4endl;
    G4cout<<"distances: "<<dm<<" "<<da<<G4endl;
    G4cout<<Pm<<G4endl;
    G4cout<<Pi<<G4endl;
    G4cout<<Pa<<G4endl;
    G4cout<<G4endl;
  }

}

//-----------------------------------------------------------//
// Need to write a shitty minimizer...
//-----------------------------------------------------------//
G4V3 TRTool::findIntPoint(G4V3 Pm, G4V3 Pa)
{

  // Decide some tolerance
  G4double tolerance = 1e-6;
  G4double diff = 9999;
  bool prevStepRight = false;
  bool nextStep = false;
  G4double step = 2.0; // m

  // first guess defined by the line
  G4double xi = 0;
  if( fabs(Pa.x() - Pm.x()) > 1e-4 ) 
    xi = Pm.x() - Pm.z() / ((Pa.z()-Pm.z())/(Pa.x()-Pm.x()));

  while( diff > tolerance ){
    
    xi = findBestRefractedX(Pm.x(),Pm.z(),Pa.x(),Pa.z(),
			    xi, step, diff, nextStep);
    
    if(prevStepRight != nextStep){
      step /= 2.;
    }
    prevStepRight = nextStep;

    // for debugging
    if(false){
      G4cout<<"-------------------------------------------"<<G4endl;
      G4cout<<"Pm:  "<<Pm<<G4endl;
      G4cout<<"Pa:  "<<Pa<<G4endl;
      G4cout<<"Pix: "<<xi<<G4endl;
      G4cout<<"tol: "<<tolerance<<" where we are: "<<diff<<G4endl;
      G4cout<<"Step size: "<<step<<G4endl;
    }
  }    
  
  // Return te interaction point. Again ony in 2-D here
  return G4V3(xi, 0, 0);

}

//-----------------------------------------------------------//
// stupid minimizer
//-----------------------------------------------------------//
G4double TRTool::findBestRefractedX(G4double xm, G4double zm,
				    G4double xa, G4double za,
				    G4double xi, G4double step,
				    G4double &result,
				    bool &rightStep)
{
  
  G4double step_right = getdOP(xm,zm,xi+step,m_n0) + getdOP(xa,za,xi+step,m_n1);
  G4double step_left = getdOP(xm,zm,xi-step,m_n0) + getdOP(xa,za,xi-step,m_n1);

  if(false){
    G4cout<<"**********************************"<<G4endl;
    G4cout<<"Step right: "<<fabs(step_right)<<G4endl;
    G4cout<<"Step left: "<<fabs(step_left)<<G4endl;
    G4cout<<"xi: "<<xi<<" +/- "<<step<<G4endl;
    G4cout<<"m:  "<<xm<<" "<<zm<<G4endl;
    G4cout<<"a:  "<<xa<<" "<<za<<G4endl;
  }

  if( fabs(step_right) < fabs(step_left) ){
    rightStep = true;
    result    = fabs(step_right);
    return xi + step;
  }
  else{
    rightStep = false;
    result    = fabs(step_left);
    return xi - step;
  }

}

//-----------------------------------------------------------//
// Calculate dOp/dx
//-----------------------------------------------------------//
G4double TRTool::getdOP(G4double x0, G4double z0, G4double x1,
			G4double n)
{

  return n * (x0 - x1) / sqrt(pow(x1-x0,2)+z0*z0);

}

//-----------------------------------------------------------//
// Setup method to get the fresnel corrected field
//-----------------------------------------------------------//
G4V3 TRTool::getFresnelCorrectedField(G4V3 E, 
				      G4double theta_i,
				      G4double theta_r,
				      TRRay trr_type)
{

  if(trr_type == TRR_direct) return E;

  G4V3 Es = G4V3(0,0,E.z());
  G4V3 Ep = G4V3(E.x(),E.y(),0);
  
  if(trr_type == TRR_reflected)
    return getReflectedParallel(Ep, theta_i) + getReflectedPerp(Es, theta_i);
  
  else if(trr_type == TRR_refracted){
    G4V3 Enew = getRefractedParallel(Ep, theta_r) + getRefractedPerp(Es, theta_r);
    //G4cout<<theta_r * 180/m_pi<<" "<<Enew<<G4endl;
    return Enew;
  }
  else{
    G4cout<<"What the hell is going on! You are trying"<<G4endl;
    G4cout<<"to access something else besides direct, "<<G4endl;
    G4cout<<"reflected, or refracted"<<G4endl;
    return E;
  }


}

//-----------------------------------------------------------//
// Get refracted Efield parallel
// Note these contain the factor for divergence
//-----------------------------------------------------------//
G4V3 TRTool::getRefractedParallel(G4V3 Ep,
				  G4double theta_r)
{

  float snellterm =  m_n1/m_n0 * sin(theta_r);
  if( snellterm > 1 ) snellterm = 1;

  double num = 2*m_n1*cos(theta_r);
  double den = m_n1*sqrt(1-pow(snellterm,2)) + m_n0*cos(theta_r);
  
  return Ep.mag() * num/den * Ep.unit();

}

//-----------------------------------------------------------//
// Get refracted Efield perpendicular
// Note these contain the factor for divergence
//-----------------------------------------------------------//
G4V3 TRTool::getRefractedPerp(G4V3 Es,
			      G4double theta_r)
{

  //if( m_n1/m_n0 * sin(theta_r) > 1 ) return G4V3(0,0,0);
  float snellterm =  m_n1/m_n0 * sin(theta_r);
  if( snellterm > 1 ) snellterm = 1;

  double num = 2*m_n1*cos(theta_r);
  double den = m_n0*sqrt(1-pow(snellterm,2)) + m_n1*cos(theta_r);
  
  return Es.mag() * num/den * Es.unit();

}

//-----------------------------------------------------------//
// Get the reflected Efield parallel
//-----------------------------------------------------------//
G4V3 TRTool::getReflectedParallel(G4V3 Ep,
				  G4double theta_i)
{

  // Total reflection, forget Fresnel...
  if( m_n1/m_n0 * sin(theta_i) > 1 )  return Ep;  

  double num = m_n1 * cos(theta_i) - m_n0 * sqrt(1 - pow(m_n0/m_n1 * sin(theta_i),2));
  double den = m_n1 * cos(theta_i) + m_n0 * sqrt(1 - pow(m_n0/m_n1 * sin(theta_i),2));

  return Ep.mag() * num/den * Ep.unit();

}

//-----------------------------------------------------------//
// Get the reflected Efield perpendicular
//-----------------------------------------------------------//
G4V3 TRTool::getReflectedPerp(G4V3 Es,
			      G4double theta_i)
{

  // Total reflection, forget Fresnel...
  if( m_n1/m_n0 * sin(theta_i) > 1 )  return Es;
  
  double num = m_n0 * cos(theta_i) - m_n1 * sqrt(1 - pow(m_n0/m_n1 * sin(theta_i),2));
  double den = m_n0 * cos(theta_i) + m_n1 * sqrt(1 - pow(m_n0/m_n1 * sin(theta_i),2));

  return Es.mag() * num/den * Es.unit();

}
