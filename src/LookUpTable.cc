#include "LookUpTable.hh"
#include "Constants.hh"
#include <fstream>
//--------------------------------------------------//
// Constructor
//--------------------------------------------------//
LookUpTable::LookUpTable(std::string filename)
{
  std::ifstream lookupfile;
  lookupfile.open(filename.c_str(), std::ios::in);
  if (!lookupfile.good()){
    lookupfile.close();
    std::cout << " FILE IS NOT GOOD" << std::endl;
  }
  FillInfo(lookupfile);
  FillData(lookupfile);
  
  lookupfile.close();
}

//--------------------------------------------------//
// Destructor
//--------------------------------------------------//
LookUpTable::~LookUpTable()
{
  // dummy
}
//Fill the basic info in the header of the lookup file
//position  of the antenna, the tilt, the dimension of ice, the grid for each direction
void LookUpTable::FillInfo(std::ifstream& lookup){
  std::string buf = "";
  getline(lookup, buf);
  lookup >> m_xant >> m_yant >> m_zant >> m_icetilt >> m_xmin >> m_xmax >> m_ymin >> m_ymax >> m_zmin >> m_zmax >> m_gridx >> m_gridy >> m_gridz;
  float floatpointx = (m_xmax - m_xmin)/m_gridx;
  float floatpointy = (m_ymax - m_ymin)/m_gridy;
  float floatpointz = (m_zmax - m_zmin)/m_gridz;
  m_nrpointX  = floor(floatpointx);
  m_nrpointY  = floor(floatpointy);
  m_nrpointZ  = floor(floatpointz);
} 

//fill the big array of point structure (x,y,z)
//an array with z coord is filled first, then these vectors are filled for each y
// then these 2D vector are filled for each x
void LookUpTable::FillData(std::ifstream& lookup){
  for (int i = 0 ; i < m_nrpointX ; ++i){
    std::vector <std::vector <point> > yzice_v;
    for (int j = 0 ; j < m_nrpointY ; ++j){
      std::vector <point> zice_v;
      for (int k = 0 ; k < m_nrpointZ ; ++k){
	G4double x,y,z;
	lookup >> x >>  y  >> z ;
 	point thepoint;
	thepoint.x = x;
	thepoint.y = y;
	thepoint.z = z;
	zice_v.push_back(thepoint);
      }
      yzice_v.push_back(zice_v);
    }    
    m_points.push_back(yzice_v);
  }
  
}

//look up table were filled for an horizontal ice block 
//(i.e. the x,y,z were running over -0.5->0.5/-0.15->0.15/0->0.3)
//We need first to rotate the point in these coordinate 
//before looking at the index i.e. the position of the data in the file.
point LookUpTable::gethoripoint(G4double xice, G4double yice, G4double zice){
  point horipoint;
  horipoint.x = xice*cos(-m_icetilt*m_pi/180) + zice*sin(-m_icetilt*m_pi/180);
  horipoint.y = yice;
  horipoint.z = -xice*sin(-m_icetilt*m_pi/180) + zice*cos(-m_icetilt*m_pi/180);
  return horipoint;  
}

G4ThreeVector LookUpTable::getG4horipoint(G4double xice, G4double yice, G4double zice){
  point horipoint;
  horipoint.x = xice*cos(-m_icetilt*m_pi/180) + zice*sin(-m_icetilt*m_pi/180);
  horipoint.y = yice;
  horipoint.z = -xice*sin(-m_icetilt*m_pi/180) + zice*cos(-m_icetilt*m_pi/180);
  G4ThreeVector hori = G4ThreeVector(horipoint.x,horipoint.y,horipoint.z);
  return hori;  
}

//if the point in horizontal ice is larger than the limits,
//we will just return a dummy position 
bool LookUpTable::checkpoint(point pt){
  bool insidethebox = true;
  if (pt.x < m_xmin || pt.x > m_xmax || pt.y < m_ymin || pt.y > m_ymax || pt.z < m_zmin || pt.z > m_zmax)
    insidethebox = false;
  return insidethebox;
}

//after setting the new coord in horizontal position,
//we compute the index for each dimension
point LookUpTable::getpoint(G4double xice, G4double yice, G4double zice){
  point horipoint = gethoripoint(xice, yice, zice);
  point surfpoint ;  
  if (checkpoint(horipoint)){
    int indexx = (int)((horipoint.x - m_xmin)/m_gridx);
    int indexy = (int)((horipoint.y - m_ymin)/m_gridy);
    int indexz = (int)((horipoint.z - m_zmin)/m_gridz);
    surfpoint = m_points[indexx][indexy][indexz];
  }
  else{
    surfpoint.x = -1;
    surfpoint.y = -1;
    surfpoint.z = -1;
  } 
  return surfpoint;
}


G4ThreeVector LookUpTable::getG4vector(G4ThreeVector pos){
  point surfpoint = getpoint(pos.x(), pos.y(), pos.z());
  return  G4ThreeVector (surfpoint.x,surfpoint.y, surfpoint.z);
}

G4double LookUpTable::getX(G4double xice, G4double yice, G4double zice){
  point thepoint = getpoint(xice, yice, zice);
  return thepoint.x;
}

G4double LookUpTable::getY(G4double xice, G4double yice, G4double zice){
  point thepoint = getpoint(xice, yice, zice);
  return thepoint.y;
}

G4double LookUpTable::getZ(G4double xice, G4double yice, G4double zice){
  point thepoint = getpoint(xice, yice, zice);
  return thepoint.z;
}

void LookUpTable::printinfo(){
  std::cout << "antpos x,y,z = " << m_xant << " " << m_yant << " " << m_zant << std::endl;
  std::cout << "ice tilt = " << m_icetilt << std::endl;  
  std::cout << " xmin xmax = " << m_xmin << " " << m_xmax  << " ymin ymax = " << m_ymin << " " << m_ymax << " zmin zmax = " << m_zmin << " " << m_zmax << std::endl;
  std::cout << "nr of points x,y,z = " << m_nrpointX << " " << m_nrpointY << " " << m_nrpointZ << std::endl;
  std::cout << "grid x,y,z = " << m_gridx << " " << m_gridy << " " << m_gridz << std::endl;
}


void LookUpTable::printtable(){
  for (int i = 0 ; i < m_nrpointX ; ++i){
    for (int j = 0 ; j < m_nrpointY ; ++j){
      for (int k = 0 ; k < m_nrpointZ ; ++k){
	std::cout << "x = " << m_points[i][j][k].x << std::endl;
	std::cout << "y = " << m_points[i][j][k].y << std::endl;
	std::cout << "z = " << m_points[i][j][k].z << std::endl;
      }
    }    
  }
}
