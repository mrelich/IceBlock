
#include "Antenna.hh"

//--------------------------------------------------//
// Constructor
//--------------------------------------------------//
Antenna::Antenna(G4double x, G4double y, G4double z,
		 G4double tStart, G4int nPoints,
		 G4double stepSize,
		 G4double angle, G4double refAngle,
		 G4double zprime)
{

  m_x = x;
  m_y = y;
  m_z = z;

  m_angle    = angle;
  m_refAngle = refAngle;

  m_zprime = zprime;

  m_tStart   = tStart;
  m_nPoints  = nPoints;
  m_stepSize = stepSize;

  // Initialize the vector potential vectors
  m_Ax = std::vector<G4double> (m_nPoints, 0.);
  m_Ay = std::vector<G4double> (m_nPoints, 0.);
  m_Az = std::vector<G4double> (m_nPoints, 0.);

  // Initialize the vector for efield
  m_Ex = std::vector<G4double> (m_nPoints, 0.);
  m_Ey = std::vector<G4double> (m_nPoints, 0.);
  m_Ez = std::vector<G4double> (m_nPoints, 0.);

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
// Add a point
//--------------------------------------------------//
void Antenna::addEPoint(unsigned int ip,
			G4double Ex,
			G4double Ey,
			G4double Ez)
{

  // Add the point to the list
  m_Ex[ip] += Ex;
  m_Ey[ip] += Ey;
  m_Ez[ip] += Ez;

}

//--------------------------------------------------//
// Accessor for points
//--------------------------------------------------//
void Antenna::getEPoint(unsigned int i,
			G4double &time,
			G4double &Ex,
			G4double &Ey,
			G4double &Ez)
{

  if( i >= m_Ex.size() ){
    G4cout<<"*** ERROR: Trying to access a point "
	 <<"that is outside range."<<G4endl;
    G4cout<<"Setting vector potential to 0."<<G4endl;
    Ex = Ey = Ez = time = 0;
    return;
  }

  Ex = m_Ex[i];
  Ey = m_Ey[i];
  Ez = m_Ez[i];
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

  for(unsigned int i=0; i<m_Ex.size(); ++i){
    m_Ex[i] = 0;
    m_Ey[i] = 0;
    m_Ez[i] = 0;
  }

}
