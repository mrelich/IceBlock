
#include "RefractionTool.hh"
#include <fstream>
#include <sstream>

using namespace std;

//--------------------------------------------------------//
// Constructor
//--------------------------------------------------------//
RefractionTool::RefractionTool() :
  m_initialRot(NULL),
  m_backRot(NULL),
  m_identity(NULL),
  m_n0(1.78),
  m_n1(1.00),
  m_tolerance(1),
  m_initialized(false),
  m_useTool(false),
  m_useLookup(true)
{


}


//--------------------------------------------------------//
// Destructor
//--------------------------------------------------------//
RefractionTool::~RefractionTool()
{

  // This causes things to crash... wtf?
  // Root weirdness?
  //if(m_initialRot) m_initialRot->Delete();
  //if(m_backRot)    m_backRot->Delete();
  //if(m_identity)   m_identity->Delete();

}

//--------------------------------------------------------//
// Initialize the tool
//--------------------------------------------------------//
void RefractionTool::initialize(G4V3 planeNorm,
				G4V3 blockCenter,
				G4V3 blockDim,
				double index0,
				double index1,
				vector<Antenna*>* ants)
{

  // Set vectors and coordinates
  m_planeNorm   = G4V3(planeNorm.x(), planeNorm.y(), planeNorm.z());
  m_blockCenter = G4V3(blockCenter.x(), blockCenter.y(), blockCenter.z());
  m_blockDim    = G4V3(blockDim.x(), blockDim.y(), blockDim.z());

  // Set indices of refraction for two materials
  m_n0 = index0;
  m_n1 = index1;

  // Set up identity matrix
  m_identity = new G4RotationMatrix(G4V3(1,0,0),
				    G4V3(0,1,0),
				    G4V3(0,0,1));
				    

  //TMatrixT<double> identity = TMatrixT<double>(3,3);
  //identity(0,0) = 1.;
  //identity(1,1) = 1.;
  //identity(2,2) = 1.;
  //m_identity = new TMatrixT<double>( identity );
  
  // Set the rotation matrices
  setInitialRotation();
  setBackRotation();
  
  // Setup the z shift
  m_zshift.set(0,0, m_blockCenter.z() + m_blockDim.z()/2.);

  // Initialize the lookup tables for each antenna
  // TODO: Try to handle this better...
  G4int icetilt = G4int(asin(planeNorm.x())*180/m_pi);
  stringstream ss;
  for(unsigned int iA=0; iA<ants->size(); ++iA){
    G4double zpos = ants->at(iA)->getZ();
    ss.str("");
    ss << "data/lookup/lookup_ice_" << icetilt << "_zant_" << zpos << ".txt";
    ifstream infile(ss.str().c_str());
    if( !infile.good() ){
      G4cout<<"Lookup table doesn't exist: "<<ss.str()<<G4endl;
      G4cout<<"Please check path. Code may crash."<<G4endl;
      m_useLookup = false;
      break;
    }
    m_lookups.push_back( new LookUpTable(ss.str()) );
    
  }// end loop over antennas

  G4cout<<"Set use tool: "<<m_useLookup<<G4endl;

  // Set initialization flag
  m_initialized = true;

}

//--------------------------------------------------------//
// Get the Interaction Point using lookup tables
//--------------------------------------------------------//
G4V3 RefractionTool::getIntPoint(G4V3 pt, G4int iAnt, 
				 G4double &theta_i,
				 G4double &theta_r,
				 G4double &R,
				 G4bool inice,
				 G4bool direct)
{

  G4int lsize = m_lookups.size();
  if( iAnt > lsize ){
    G4cout<<"Lookup table doesn't exist for antenna!!!"<<G4endl;
    G4cout<<"Just returning the track middle point!!! "<<G4endl;
    return pt;
  }
  

  // Get the interaction Point
  //G4V3 intPoint = m_lookups.at(iAnt)->getG4vector(pt);
  G4V3 intPoint;
  if( inice ){
    intPoint = m_lookups.at(iAnt)->geticevector(pt, true, direct);
    G4V3 temp = m_lookups.at(iAnt)->geticevector(pt, false, direct);
    theta_i = temp.x();  
    theta_r = temp.y();  
    R       = temp.z();  
  }
  else{
    intPoint = m_lookups.at(iAnt)->getairvector(pt,true);
    G4V3 temp = m_lookups.at(iAnt)->getairvector(pt, false);
    theta_i = temp.x();
    theta_r = temp.y();
    R       = temp.z();
  }

  // In the lookup tables, if an
  // entry is not possible the point 
  // is set to -1,-1,-1
  if( intPoint.x() == -1 && intPoint.y() == -1 && intPoint.z() == -1 )
    intPoint.set(-9999,-9999,-9999);

  return intPoint;
}


//--------------------------------------------------------//
// Get the Interaction Point
//--------------------------------------------------------//
G4V3 RefractionTool::getIntPoint(G4V3 g4_pt,
				 G4V3 g4_pa,
				 G4double &theta_i,
				 G4double &theta_r)
{

  G4V3 pt = G4V3(g4_pt.x(), g4_pt.y(), g4_pt.z());
  G4V3 pa = G4V3(g4_pa.x(), g4_pa.y(), g4_pa.z());

  // Rotate the points to get pt' and pa'
  G4V3 ptp = *m_initialRot * pt;
  G4V3 pap = *m_initialRot * pa;

  // Translate the points such that the top face
  // of the iceblock sits in the x-y plane.
  // This gives pt'' and pa''
  G4V3 ptpp = ptp - m_zshift;
  G4V3 papp = pap - m_zshift;

  // Now find the interaction point
  G4V3 pipp = scanIntPoint(ptpp,papp);

  // Now apply a check that it complies with
  // Snell's law by checking incident angle
  // and refracted angle
  theta_i = atan(sqrt(pow(ptpp.x()-pipp.x(),2)+pow(ptpp.y()-pipp.y(),2))/fabs(ptpp.z()));
  theta_r = atan(sqrt(pow(papp.x()-pipp.x(),2)+pow(papp.y()-pipp.y(),2))/fabs(papp.z()));
  G4double expectedAngle = ( theta_i < asin(m_n1/m_n0) ) ? asin(m_n0*sin(theta_i)/m_n1) : 9999;
  
  //G4cout<<"Angles: "<<theta_i<<" "<<theta_r<<G4endl;
  //G4cout<<"Snell Check: "<<m_n0*sin(theta_i)<<" "<<m_n1*sin(theta_r)<<G4endl;
  //G4cout<<"\t"<<fabs(expectedAngle - theta_r) * 180/m_pi<<G4endl;

  if( theta_i >= asin(m_n1/m_n0) || fabs(expectedAngle - theta_r) * 180/m_pi > m_tolerance ){
    pipp.set(-9999,-9999,-9999);
    theta_i = -9999;
    theta_r = -9999;
    return G4V3(pipp.x(), pipp.y(), pipp.z());
  }

  // Now rotate and translate back along z axis
  pipp += m_zshift;
  pipp = *m_backRot * pipp;
  return G4V3(pipp.x(), pipp.y(), pipp.z());
  
}

//--------------------------------------------------------//
// Get the midpoint where the track cross the boundary
//--------------------------------------------------------//
G4V3 RefractionTool::getIntMidPoint(G4V3 P0, G4V3 P1)
{

  // Rotate the points to get iceblock parallel to z
  G4V3 P0p = *m_initialRot * P0;
  G4V3 P1p = *m_initialRot * P1;

  // Translate the points such that the bottom face
  // of the iceblock sits in the x-y plane.
  G4V3 zshift = m_zshift - G4V3(0,0, m_blockDim.z());
  G4V3 P0pp = P0p + zshift;
  G4V3 P1pp = P1p + zshift;

  // Now solve for x-y given z = 0
  // First get unit vector in direction of
  // track.
  G4V3 v     = (P1pp - P0pp).unit();

  // Now get the slope to z = 0
  G4double t = P0pp.z() / v.z();

  // Now we can calculate x and y
  G4V3 Pmpp = G4ThreeVector( P0pp.x() - t * v.x(),
			     P0pp.y() - t * v.y(),
			     0 );

  // Now rotate back to original coordinate system and return
  Pmpp -= zshift;  
  Pmpp = *m_backRot * Pmpp;
  return Pmpp;
		


}

//--------------------------------------------------------//
// Scan for intersection point
//--------------------------------------------------------//
G4V3 RefractionTool::scanIntPoint(G4V3 pt, G4V3 pa)
{

  // For this scan I assume that the pt and pa have been
  // rotated and translated such that it is sitting in
  // the x-y plane (ie. no scanning of z).

  G4V3 pi = G4V3(0,0,0); // The return value

  // TODO: This should be configurable!!!

  G4int nxsteps  = 400;
  double xmin  = m_blockCenter.x() - m_blockDim.x()/2.;
  double xmax  = m_blockCenter.x() + m_blockDim.x()/2.;
  double xstep = (xmax-xmin)/nxsteps;

  G4int nysteps  = 120; // only 30cm wide
  double ymin  = m_blockCenter.y() - m_blockDim.y()/2.;
  double ymax  = m_blockCenter.y() + m_blockDim.y()/2.;
  double ystep = (ymax-ymin)/nysteps;

  // Loop over coordinates and minimize the sum of the gradient
  // I am not sure if this is the best way, or maybe minimizing
  // each individually??
  double gradSum = 9999;
  G4V3 temp = G4V3(0,0,0);
  for(G4int ix=0; ix<nxsteps; ++ix){
    double xi = xmin + ix * xstep;
    
    for(G4int iy=0; iy<nysteps; ++iy){
      double yi = ymin + iy * ystep;

      temp.set(xi,yi,0);
      double gx = fabs( getGradient(pt,pa,temp,0) );
      double gy = fabs( getGradient(pt,pa,temp,1) );

      if( gx + gy <= gradSum ){
	gradSum = gx + gy;
	pi.set(xi,yi,0);
      }

    }// end loop over y points
  }// end loop over x points

  return pi;

}

//--------------------------------------------------------//
// Get the gradient for one coordinate
// needed for intersection scanning
// pt -- track point
// pa -- antenna position
// pi -- potential inter. point
// opt -- 0=x, 1=y, 1=z (axis)
//--------------------------------------------------------//
double RefractionTool::getGradient(G4V3 pt, G4V3 pa,
				   G4V3 pi, G4int opt)
{

  // If this method fails at any point, will return 9999
  double grad = 9999;

  // Decide coordinate we are looking at
  double ct = 0, ca = 0, ci = 0;
  if(opt == 0){
    ct = pt.x();
    ca = pa.x();
    ci = pi.x();
  }
  else if(opt == 1){
    ct = pt.y();
    ca = pa.y();
    ci = pi.y();
  }
  else if(opt == 2){
    ct = pt.z();
    ca = pa.z();
    ci = pi.z();
  }
  else return grad;

  // The point should lie in between the two points
  // otherwise there is no solution
  if( !((ct <= ci && ci <= ca) || (ca <= ci && ci <= ct)) )
    return grad;

  // Now we basically use the geometry to enforce snells law.
  // The gradient (if you work it out) is just the difference
  // of n0 * sin(theta0) - n1 * sin(theta1).

  // First the constant factors
  if( (pt - pi).mag() == 0 ) return grad;
  if( (pa - pi).mag() == 0 ) return grad;
  
  double t_const = m_n0 / (pt-pi).mag();
  double a_const = m_n1 / (pa-pi).mag();
  
  // Now calculate gradient
  grad = t_const * fabs(ct-ci) - a_const * fabs(ca-ci);

  return grad;


}

//--------------------------------------------------------//
// Set the intial rotation matrix
//--------------------------------------------------------//
void RefractionTool::setInitialRotation()
{

  //m_initialRot = getRotationMatrix(G4V3(0,0,1), m_planeNorm);
  m_initialRot = getRotationMatrix(m_planeNorm,G4V3(0,0,1));

}

//--------------------------------------------------------//
// Set the final rotation matrix (to go backwards)
//--------------------------------------------------------//
void RefractionTool::setBackRotation()
{

  //m_backRot = getRotationMatrix(m_planeNorm, G4V3(0,0,1));
  m_backRot = getRotationMatrix(G4V3(0,0,1), m_planeNorm);

}

//--------------------------------------------------------//
// Get generic rotation matrix
//--------------------------------------------------------//
//TMatrixT<double>* RefractionTool::getRotationMatrix(G4V3 v0,
//						    G4V3 v1)
G4RotationMatrix* RefractionTool::getRotationMatrix(G4V3 v0,
						    G4V3 v1)
{

  // To do this we first calculate some pieces. This is to take
  // v1 and rotate it into v0.
  // R = I + [v]_x + [v]_x^2 * (1-C)/S^2
  // v = v1 x v0
  // s = ||v|| 
  // c = v1 * v0
  // [v]_x = skew symmetrix matrix

  // Get V
  G4V3 v = v1.cross( v0 );
  
  // Get S
  G4double s = v.mag();
  
  // Get C
  G4double c = v1.dot( v0 );

  // Handle special case of v0 || v1 already
  if( c == 1 ) return m_identity;
  
  // Setup the skew symmetric matrix
  /*
  TMatrixT<double> v_x = TMatrixT<double>(3,3);
  v_x(0,1) = -v[2];
  v_x(0,2) = v[1];
  v_x(1,0) = v[2];
  v_x(1,2) = -v[0];
  v_x(2,0) = -v[1];
  v_x(2,1) = v[0];
  */
  
  // This sucks balls. The rotation matrix doesn't have the add operator
  // and I don't feel like writing my own class and including it.. so do
  // it manually here
  G4double cs = (1-c)/(s*s);
  G4V3 row0 = G4V3(1 + cs * (-pow(v[2],2)-pow(v[1],2)),
		   -v[2] + cs*v[0]*v[1],
		   v[1] + cs*v[0]*v[2]);
  G4V3 row1 = G4V3(v[2] - cs*v[0]*v[1],
		   1 + cs*(-pow(v[2],2)-pow(v[0],2)),
		   -v[0] + cs*v[1]*v[2]);
  G4V3 row2 = G4V3(-v[1] - cs*v[0]*v[2],
		   v[0] + cs*v[1]*v[2],
		   1 + cs*(-pow(v[0],2)-pow(v[1],2)));

  return new G4RotationMatrix(row0,row1,row2);

  // Setup the rotation matrix
  //return new TMatrixT<double>(*m_identity + v_x + v_x*v_x*((1-c)/(s*s)));
  
}

//--------------------------------------------------------//
// Get the rotated electric field following laws of
// refraction. NOTE: At the moment we are assuming plane
// wave condition, so only 1 angle is supported.
// To update, need to update getRefractedPerp and 
// getRefractedParallel, as they use Snell's condition
// to simplify the math
//--------------------------------------------------------//
G4V3 RefractionTool::getTransmittedField(G4V3 g4_E, 
					 G4double theta_i,
					 G4bool iceToAir)
{

  G4V3 E = G4V3(g4_E.x(), g4_E.y(), g4_E.z());

  // Get the inclinatin of the ice
  double tilt     = atan( m_planeNorm.x() / m_planeNorm.z() );

  // Get how much to rotate around y-axis
  double rotation = m_pi/2 - tilt;

  // Rotate the E-field
  G4V3 Eprime = getRotatedE(E, rotation, false);

  // Break into Es and Ep
  G4V3 Es = G4V3(0,0,Eprime.z());
  G4V3 Ep = G4V3(Eprime.x(),Eprime.y(),0);

  // Get adjusted Es and Ep
  G4V3 Es_prime = getRefractedPerp(Es, theta_i, iceToAir);
  G4V3 Ep_prime = getRefractedParallel(Ep, theta_i, iceToAir);

  // Now sum and rotate back
  G4V3 E_prime = getRotatedE(Es_prime+Ep_prime, rotation, true);

  // Return electric field
  return G4V3(E_prime.x(), E_prime.y(), E_prime.z());


}


//--------------------------------------------------------//
// Get the rotated E-field around y-axis
//--------------------------------------------------------//
G4V3 RefractionTool::getRotatedE(G4V3 E, 
				 double rotation,
				 bool back)
{

  int sign = back ? -1 : 1;
  //TMatrixT<double> rotMat = getRotationMatrixY(sign*rotation);
  G4RotationMatrix rotMat = getRotationMatrixY(sign*rotation);
  return rotMat * E;

}

//--------------------------------------------------------//
// Get Rotation matrix around y-axis
// This is counter clockwise rotation around y-axis
//--------------------------------------------------------//
//TMatrixT<double> RefractionTool::getRotationMatrixY(double beta)
G4RotationMatrix RefractionTool::getRotationMatrixY(double beta)
{

  /*TMatrixT<double> rotMat = TMatrixT<double>(3,3);
  rotMat(0,0) = cos(beta);
  rotMat(0,2) = -sin(beta);
  rotMat(1,1) = 1;
  rotMat(2,0) = sin(beta);
  rotMat(2,2) = cos(beta);
  */
  G4RotationMatrix rotMat = G4RotationMatrix( G4V3(cos(beta),0,sin(beta)),
					      G4V3(0,1,0),
					      G4V3(-sin(beta),0,cos(beta)) );
  return rotMat;

}

//--------------------------------------------------------//
// Get refracted E-field -- perpindicular component
//--------------------------------------------------------//
G4V3 RefractionTool::getRefractedPerp(G4V3 Es, 
				      G4double theta_i,
				      G4bool iceToAir)
{

  double n0 = m_n0;
  double n1 = m_n1;
  if(!iceToAir){
    n0 = m_n1;
    n1 = m_n0;
  };

  double Esmag = Es.mag();
  G4V3 Es_unit = Es.unit();
  double num   = 2*n0*cos(theta_i);
  double den   = n0*cos(theta_i) + sqrt(pow(n1,2)-pow(n0*sin(theta_i),2));

  return (Esmag * num / den) * Es_unit * 
    getCorrectionFactor(n0, n1, theta_i);

}

//--------------------------------------------------------//
// Get refracted E-field -- parallel component
//--------------------------------------------------------//
G4V3 RefractionTool::getRefractedParallel(G4V3 Ep, 
					  G4double theta_i,
					  G4bool iceToAir)
{

  double n0 = m_n0;
  double n1 = m_n1;
  if(!iceToAir){
    n0 = m_n1;
    n1 = m_n0;
  };

  double Epmag = Ep.mag();
  G4V3 Ep_unit = Ep.unit();
  double num   = 2*n0*n1*cos(theta_i);
  double den   = pow(n1,2)*cos(theta_i) + n0*sqrt(pow(n1,2)-pow(n0*sin(theta_i),2));

  return (Epmag * num / den) * Ep_unit *
    getCorrectionFactor(n0, n1, theta_i);;

}

//--------------------------------------------------------//
// Get correction factor
//--------------------------------------------------------//
G4double RefractionTool::getCorrectionFactor(G4double n0,
					     G4double n1, 
					     G4double theta_i)
{

  // The correction factor is only valid when Snell's law is safe.
  // So impose that condition here
  if( theta_i > asin(n1/n0) ) return 0;
  
  G4double factor = (n1/n0) * sqrt(1-pow(n0*sin(theta_i)/n1,2)) / cos(theta_i);
  return factor;

}
