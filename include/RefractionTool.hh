#ifndef RefractionTool_hh
#define RefractionTool_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This class is meant to take as input the ice block geometry, the //
// antenna position, and the track midpoint position, then give the //
// interaction point on the ice surface.                            //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "Constants.hh"
//#include "TVector3.h"
//#include "TMatrixT.h"


typedef G4ThreeVector G4V3;
//typedef TVector3 G4V3;

class RefractionTool
{
 
 public:

  // Constructor
  RefractionTool();
  
  // Destructor
  ~RefractionTool();

  // initialize
  void initialize(G4ThreeVector planeNorm,    // vector normal to the top face
		  G4ThreeVector blockCenter,  // Position of the iceblock center
		  G4ThreeVector blockDim,     // Dimensions of the block    
		  double index0,              // index of refraction inside the block
		  double index1);             // index of refraction outside the block


  // Method to get the intercept/interaction point
  G4ThreeVector getIntPoint(G4ThreeVector g4_pt,      // track midpoint
			    G4ThreeVector g4_pa,      // antenna position
			    G4double &theta_i,        // incident angle
			    G4double &theta_r);       // refracted angle
  
  // Method to retrieve fresnel rotated electric field
  G4ThreeVector getTransmittedField(G4ThreeVector g4_E, 
				    G4double theta_i);
  
  bool isInitialized(){ return m_initialized; };

  void setUse(bool use){ m_useTool = use; };
  bool useTool(){ return m_useTool; };

 protected:

  // Scan for intercept point
  G4V3 scanIntPoint(G4V3 pt, G4V3 pa);

  // Gradient calc for specific point
  double getGradient(G4V3 pt, G4V3 pa, G4V3 pi, G4int opt);

  // Set the rotation matrix based on the normal vector
  void setInitialRotation();

  // Set rotation matrix to take you back
  void setBackRotation();

  // Get a generic rotation metrix basied on rotating
  // one vector into the other
  //G4RotationMatrix getRotationMatrix(G4V3 v0,   // Desired direction
  //G4V3 v1);  // Vector to be rotated
  //TMatrixT<double>* getRotationMatrix(G4V3 v0,   // Desired direction
  //				      G4V3 v1);  // Vector to be rotated
  G4RotationMatrix* getRotationMatrix(G4V3 v0,   // Desired direction
  				      G4V3 v1);  // Vector to be rotated


  
  // Get rotated electric field, which in this case 
  // is much simpler. We will rotate around the y-axis
  G4V3 getRotatedE(G4V3 E,              // E-field to be rotated
		   double rotation,     // rotation angle
		   bool back);          // initial rotation or rotating back

  
  // Get rotation matrix around y-axis
  //TMatrixT<double> getRotationMatrixY(double beta);
  G4RotationMatrix getRotationMatrixY(double beta);

  // Get Refracted field in perpendicular region
  G4V3 getRefractedPerp(G4V3 Es,      // E-field perpendicular to the plane
			double theta_i);  // angle of incidence
  
  // Get Refracted field in parallel region
  G4V3 getRefractedParallel(G4V3 Ep,         // E-field parallel to the plane
			    double theta_i); // angle of incidence


  G4RotationMatrix* m_initialRot;    // Initial rotation matrix
  G4RotationMatrix* m_backRot;       // Rotate back
  G4RotationMatrix* m_identity;      // Identity matrix
  //TMatrixT<double>* m_initialRot;    // Initial rotation matrix
  //TMatrixT<double>* m_backRot;       // Rotate back
  //TMatrixT<double>* m_identity;      // Identity matrix

  G4V3 m_planeNorm;                 // Vector for normal to plane
  G4V3 m_blockCenter;               // center coordinate of the iceblock
  G4V3 m_blockDim;                  // Dimensions of iceblock

  G4V3 m_zshift;                    // Shift along z-axis during transformation

  double m_n0;                    // index of refraction of block material
  double m_n1;                    // index of refraction outside block

  double m_tolerance;             // Tolerance for how far off refracted angle you allow

  bool m_initialized;             // Keep track if tool was initialized

  bool m_useTool;                 // Internally store whether or not to use

};

#endif
