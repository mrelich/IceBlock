
#include "Antenna.hh"

//--------------------------------------------------//
// Constructor
//--------------------------------------------------//
Antenna::Antenna(G4double x, G4double y, G4double z,
		 G4double tStart, G4int nPoints,
		 G4double stepSize)
{

  m_x = x;
  m_y = y;
  m_z = z;

  m_tStart = tStart;
  m_nPoints = nPoints;
  m_stepSize = stepSize;

  // Initialize the vector potential vectors
  m_Ax = std::vector<G4double> (m_nPoints, 0.);
  m_Ay = std::vector<G4double> (m_nPoints, 0.);
  m_Az = std::vector<G4double> (m_nPoints, 0.);

}

//--------------------------------------------------//
// Destructor
//--------------------------------------------------//
Antenna::~Antenna()
{
  // dummy
}

//--------------------------------------------------//
// Get distance to some arbitrary point
//--------------------------------------------------//
G4double Antenna::getR(G4double x0, G4double y0, G4double z0)
{

  return sqrt( pow(x0-m_x,2) + 
	       pow(y0-m_y,2) +
	       pow(z0-m_z,2) );

}

//--------------------------------------------------//
// Add a point
//--------------------------------------------------//
void Antenna::addPoint(unsigned int ip,
		       G4double Ax,
		       G4double Ay,
		       G4double Az)
{

  // Add the point to the list
  m_Ax[ip] += Ax;
  m_Ay[ip] += Ay;
  m_Az[ip] += Az;

}

//--------------------------------------------------//
// Accessor for points
//--------------------------------------------------//
void Antenna::getPoint(unsigned int i,
		       G4double &time,
		       G4double &Ax,
		       G4double &Ay,
		       G4double &Az)
{

  if( i >= m_Ax.size() ){
    G4cout<<"*** ERROR: Trying to access a point "
	 <<"that is outside range."<<G4endl;
    G4cout<<"Setting vector potential to 0."<<G4endl;
    Ax = Ay = Az = time = 0;
    return;
  }

  Ax = m_Ax[i];
  Ay = m_Ay[i];
  Az = m_Az[i];
  time = m_tStart + i*m_stepSize;

}

		       
//--------------------------------------------------//
// Clear out information
//--------------------------------------------------//
void Antenna::clear()
{
  
  for(unsigned int i=0; i<m_Ax.size(); ++i){
    m_Ax[i] = 0;
    m_Ay[i] = 0;
    m_Az[i] = 0;
  }

}
