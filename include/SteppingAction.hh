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

  // Write the vector potential info from TR
  void TRFromZHS(const G4Step*);

 private:

  std::ofstream* m_output;  
  std::vector<Antenna*> *m_ants;
  //MyTreeWriter* m_treeWriter;

  // Set intial and final points for the step
  void setInitialFinalPoint(const G4Step* step,
			    G4ThreeVector &P0,
			    G4ThreeVector &P1,
			    G4double &t0,
			    G4double &t1);

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

  // Fill antennas from TR
  void fillForAntenna(Antenna* ant,
		      G4double t0,
		      G4double t1,
		      G4double tD0,
		      G4double tD1,
		      G4double dtD_dt,
		      G4ThreeVector A);


  // Tool for handling refraction
  RefractionTool* m_refTool;
  
  // Bool to tell whether or not outside point
  // has been used in TR calculation
  bool m_TRFirstPointFound;

};
#endif
