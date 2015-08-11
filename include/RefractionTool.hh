#ifndef RefractionTool_hh
#define RefractionTool_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This class is meant to take as input the ice block geometry, the //
// antenna position, and the track midpoint position, then give the //
// interaction point on the ice surface.                            //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include <vector>

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "Constants.hh"
#include "Antenna.hh"
#include "LookUpTable.hh"

typedef G4ThreeVector G4V3;

class RefractionTool
{
 
 public:

  // Constructor
  RefractionTool();
  
  // Destructor
  ~RefractionTool();

  // initialize
  void initialize(G4V3 planeNorm,               // vector normal to the top face
		  G4V3 blockCenter,             // Position of the iceblock center
		  G4V3 blockDim,                // Dimensions of the block    
		  double index0,                // index of refraction inside the block
		  double index1,                // index of refraction outside the block
		  std::vector<Antenna*>* ants); // Antennas

  // Method to get the intercept/interaction point
  G4V3 getIntPoint(G4V3 g4_pt,      // track midpoint
		   G4V3 g4_pa,      // antenna position
		   G4double &theta_i,        // incident angle
		   G4double &theta_r);       // refracted angle

  // Get interaction point from lookup tables
  G4V3 getIntPoint(G4V3 pt,             // track midpoint
		   G4int iAnt,          // antenna number		   
		   G4double &theta_i,   // incident angle
		   G4double &theta_r,   // refracted angle
		   G4double &R,         // Total optical path
		   G4bool inice=true,   // is the point in ice?
		   G4bool direct=true); // do we want direct or reflected path?

  // Get the midpoint that crosses the surface boundary.
  // This is used in transition radiation calculation
  G4V3 getIntMidPoint(G4V3 P0,       // Start point for step
		      G4V3 P1);      // End point for step

  // Method to retrieve fresnel rotated electric field
  G4V3 getTransmittedField(G4V3 g4_E, 
			   G4double theta_i,
			   G4bool iceToAir=true);
  
  bool isInitialized(){ return m_initialized; };

  void setUse(bool use){ m_useTool = use; };
  bool useTool(){ return m_useTool; };
  bool useLookup(){ return m_useLookup; };

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
  G4V3 getRefractedPerp(G4V3 Es,          // E-field perpendicular to the plane
			G4double theta_i, // angle of incidence
			G4bool iceToAir); // to deal with ice->air or air->ice 

  // Get Refracted field in parallel region
  G4V3 getRefractedParallel(G4V3 Ep,          // E-field parallel to the plane
			    G4double theta_i, // angle of incidence
			    G4bool iceToAir); // to deal with ice->air or air->ice

  G4double getCorrectionFactor(G4double n0,        // index of material containing shower
			       G4double n1,        // index of material refracting into
			       G4double theta_i);  // incident angle

  G4RotationMatrix* m_initialRot;    // Initial rotation matrix
  G4RotationMatrix* m_backRot;       // Rotate back
  G4RotationMatrix* m_identity;      // Identity matrix
  //TMatrixT<double>* m_initialRot;  // Initial rotation matrix
  //TMatrixT<double>* m_backRot;     // Rotate back
  //TMatrixT<double>* m_identity;    // Identity matrix

  G4V3 m_planeNorm;                  // Vector for normal to plane
  G4V3 m_blockCenter;                // center coordinate of the iceblock
  G4V3 m_blockDim;                   // Dimensions of iceblock

  G4V3 m_zshift;                     // Shift along z-axis during transformation

  double m_n0;                       // index of refraction of block material
  double m_n1;                       // index of refraction outside block

  double m_tolerance;                // Tolerance for how far off refracted angle you allow

  bool m_initialized;                // Keep track if tool was initialized

  bool m_useTool;                    // Internally store whether or not to use

  std::vector<LookUpTable*> m_lookups;    // Lookup table to speed up scanning

  bool m_useLookup;

};

#endif
