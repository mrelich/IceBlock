#ifndef Antenna_h
#define Antenna_h

#include "globals.hh"
#include <vector>

class Antenna
{
 public:
  
  // Constructor
  Antenna(G4double x, G4double y, G4double z,  // Antenna location
	  G4double tStart,                     // Start time for data taking
	  G4int nPoints,                       // Number of points to store
	  G4double stepSize,                   // Step size (bin width)
	  G4double angle,                      // Angle w/respect to beam direction
	  G4double refAngle,                   // Refracted Angle assuming some tilt in ice
	  G4double zprime);                     // Primed z coordinate
  // Destructor
  ~Antenna();

  
  // Methods to access some variables
  G4double getX(){ return m_x; };
  G4double getY(){ return m_y; };
  G4double getZ(){ return m_z; };
  G4double getAngle(){ return m_angle; };
  G4double getRefAngle(){ return m_refAngle; };
  G4double getZprime(){ return m_zprime; };

  // Method to calculate
  G4double getR(G4double x0, G4double y0, G4double z0);

  // Get timing info
  G4double getTmin()   { return m_tStart; };
  G4double getTmax()   { return m_tStart + m_nPoints*m_stepSize; };
  G4double getTStep()  { return m_stepSize; };
  G4double getNPoints(){ return m_nPoints; };

  // Add vector potential 
  void addPoint(unsigned int ip, 
		G4double Ax,
		G4double Ay, 
		G4double Az);

  // Get the number of points
  unsigned int getN(){ return m_Ax.size(); };
  
  // Accessor for the points
  void getPoint(unsigned int i,
		G4double &time,
		G4double &Ax,
		G4double &Ay,
		G4double &Az);

  // CLear out vector potential info
  void clear();
  
 private:
  
  G4double m_x;            // X position
  G4double m_y;            // Y position
  G4double m_z;            // Z position
  
  G4double m_angle;        // angle with respect to beam
  G4double m_refAngle;     // Angle refracted assuming some tilt in ice

  G4double m_zprime;       // Coordinate in refracted system. x is fixed.

  G4double m_tStart;       // Start time for data taking
  G4int    m_nPoints;      // Number of points to consider
  G4double m_stepSize;     // Resolution for detector

  std::vector<G4double> m_Ax; // Vector potential in X direction
  std::vector<G4double> m_Ay; // Vector potential in Y direction
  std::vector<G4double> m_Az; // Vector potential in Z direction

  
  
};

#endif
