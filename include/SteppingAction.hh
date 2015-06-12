#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//#include "MyTreeWriter.hh"
#include <fstream>
#include "Antenna.hh"
#include "RefractionTool.hh"

class SteppingAction : public G4UserSteppingAction
{

 public:

  // Constructor / Destructor
  SteppingAction(std::ofstream* output,
		 std::vector<Antenna*> *ants,
  		 RefractionTool *refTool = NULL);
  ~SteppingAction();
  
  // User Stepping action. This constrols
  // what is actually saved.
  void UserSteppingAction(const G4Step*);
  
  // Write the step information
  void WriteSteps(const G4Step*);

  // Write the vector potential info
  void VPotentialZHSStyle(const G4Step*);

  // Write the Efield from Endpoint
  void EFieldEndpointStyle(const G4Step*);

  // Find the intersection point on the surface
  // of the iceblock for a specific Antenna
  //G4ThreeVector FindIntersection(G4ThreeVector antPos,      // Antenna position
  //				 G4ThreeVector midPoint,    // Mid-point of the track
  //				 G4PhysicalVolume* volume); // IceBlock volume

 private:

  std::ofstream* m_output;  
  std::vector<Antenna*> *m_ants;
  //MyTreeWriter* m_treeWriter;

  // Get Velocity vector
  G4ThreeVector getVelocity(G4ThreeVector P0,
			    G4ThreeVector P1,
			    G4double t0,
			    G4double t1);
			
  // Set the unit vector to antenna
  // along with the distance to antenna
  void setUnitVector(G4ThreeVector v_ant,
		     G4ThreeVector v_Point,
		     G4ThreeVector &u,
		     G4double &R);

  // Get the detector time
  G4double getTDetector(G4ThreeVector ant,
			G4ThreeVector point,
			G4double time,
			G4double n_,
			G4double c_);

  // Testing
  G4double getBeta(const G4Step* step);

  // Calculate the e-field from paramters
  // for the endpoint method
  G4ThreeVector getEFieldEndpoint(G4ThreeVector Beta, // Velocity at step [dimensionless]
				  G4ThreeVector rhat, // unit vector towards antenna 
				  G4double R,         // distance to antenna [cm]
				  G4double dt,        // time-step [s]
				  G4double q,         // charge [esu]
				  G4double c);        // Speed of light [cm/s]

  // Tool for handling refraction
  RefractionTool* m_refTool;

};
#endif
