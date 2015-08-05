#ifndef TRTool_hh
#define TRTool_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// A simplified version of the refraction tool to handle the current //
// transition radiation setup, where we are studying infinite plane. //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "G4ThreeVector.hh"
#include "Constants.hh"

typedef G4ThreeVector G4V3;

class TRTool
{
 public:

  // Constructor and destructor
  TRTool(G4double n0, G4double n1);
  ~TRTool();

  // Method to find the path 
  void findPath(G4V3 Pm,             // Midpoint of the track
		G4V3 Pa,             // Antenna/observer position
		G4double t0,         // initial time in particle frame
		G4double t1,         // initial time in particle frame
		G4double &theta_i,   // incident angle
		G4double &theta_r,   // refracted angle
		G4double &distance,  // optical path taking into account inde
		G4V3     &u_r,       // unit vector
		G4double &t0_ret,    // retarded time 
		G4double &t1_ret,    // retarded time		
		TRRay trr_type);     // The ray type

  // Method to get the overall refracted field
  G4V3 getFresnelCorrectedField(G4V3 E, 
				G4double theta_i,
				G4double theta_r,
				TRRay trr_type);
  

 private:
  
  // Find direct path
  void findPathDirect(G4V3 Pm, G4V3 Pa, 
		      G4double t0, G4double t1,
		      G4double &theta_i, G4double &theta_r,
		      G4double &distance, G4V3 &u_r,
		      G4double &t0_ret, G4double &t1_ret);

  // Find reflected
  void findPathReflected(G4V3 Pm, G4V3 Pa,
			 G4double t0, G4double t1,
			 G4double &theta_i, G4double &theta_r,
			 G4double &distance, G4V3 &u_r,
			 G4double &t0_ret, G4double &t1_ret);

  // Find refracted
  void findPathRefracted(G4V3 Pm, G4V3 Pa,
			 G4double t0, G4double t1,
			 G4double &theta_i, G4double &theta_r,
			 G4double &distance, G4V3 &u_r,
			 G4double &t0_ret, G4double &t1_ret);
  
  // Find the interaction point for refracted case
  G4V3 findIntPoint(G4V3 Pm, G4V3 Pa);

  // Stupid minimizer functions
  G4double findBestRefractedX(G4double xm, G4double zm,
			      G4double xa, G4double za,
			      G4double xi, G4double step,
			      G4double &result,
			      bool &rightStep);
  

  // Get the derivative of OP path for one piece
  G4double getdOP(G4double x0, G4double z0, G4double x1, G4double n);

  // Refracted parallel and perp field
  G4V3 getRefractedParallel(G4V3 Ep, G4double theta_r);
  G4V3 getRefractedPerp(G4V3 Es, G4double theta_r);

  // Reflected parallel and perp field
  G4V3 getReflectedParallel(G4V3 Ep, G4double theta_i);
  G4V3 getReflectedPerp(G4V3 Es, G4double theta_i);

  // Indices of refraction
  G4double m_n0;
  G4double m_n1;
  

};

#endif
