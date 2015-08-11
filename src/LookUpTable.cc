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
  //  FillTestData(lookupfile);
  FillDataAirIce(lookupfile);
  
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
//position  of the antenna, the tilt, the grid for each direction
//dimension of the ice block
//dimension of the air block
void LookUpTable::FillInfo(std::ifstream& lookup){
  std::string buf = "";
  getline(lookup, buf);
  lookup >> m_xant >> m_yant >> m_zant >> m_icetilt >> m_gridx >> m_gridy >> m_gridz >> m_boundaryshift ;
  getline(lookup, buf);
  getline(lookup, buf);
  lookup >> m_icexmin >> m_icexmax >> m_iceymin >> m_iceymax >> m_icezmin >> m_icezmax ;
  getline(lookup, buf);
  getline(lookup, buf);
  lookup >> m_airxmin >> m_airxmax >>  m_airymin >> m_airymax >> m_airzmin >> m_airzmax;

  float icefloatpointx = (m_icexmax - m_icexmin - 2*m_boundaryshift)/m_gridx;
  float icefloatpointy = (m_iceymax - m_iceymin - 2*m_boundaryshift)/m_gridy;
  float icefloatpointz = (m_icezmax - m_icezmin - 2*m_boundaryshift)/m_gridz;
  //+1 come from the difference from interval and point !
  m_icenrpointX  = floor(icefloatpointx) +1;
  m_icenrpointY  = floor(icefloatpointy) +1;
  m_icenrpointZ  = floor(icefloatpointz) +1;
  //the number of reflected points is the same as direct ones
  m_ricenrpointX  = floor(icefloatpointx) +1;
  m_ricenrpointY  = floor(icefloatpointy) +1;
  m_ricenrpointZ  = floor(icefloatpointz) +1;

  float airfloatpointx = (m_airxmax - m_airxmin - 2*m_boundaryshift)/m_gridx;
  float airfloatpointy = (m_airymax - m_airymin - 2*m_boundaryshift)/m_gridy;
  float airfloatpointz = (m_airzmax - m_airzmin - 2*m_boundaryshift)/m_gridz;
  m_airnrpointX  = floor(airfloatpointx) + 1;
  m_airnrpointY  = floor(airfloatpointy) + 1;
  m_airnrpointZ  = floor(airfloatpointz) + 1;
} 

//fill the big array of point structure (x,y,z)
//an array with z coord is filled first, then these vectors are filled for each y
// then these 2D vector are filled for each x
//fill tables for only ice table. see next function for ice and air
void LookUpTable::FillData(std::ifstream& lookup){
  for (int i = 0 ; i < m_icenrpointX ; ++i){
    std::vector <std::vector <point> > yzice_v;
    for (int j = 0 ; j < m_icenrpointY ; ++j){
      std::vector <point> zice_v;
      for (int k = 0 ; k < m_icenrpointZ ; ++k){
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
    m_icepoints.push_back(yzice_v);
  }
  
}

void LookUpTable::FillDataAirIce(std::ifstream& lookup){
  // fill ice array
  // there are 2 "points" for ice: first is the direct point, second is the reflected
  // for the direct the angles are the incident and refracted
  // for the reflected path the angles are the first incident 
  // (which is the same as the first reflected and the same as the second incident)
  // and the refracted angle when the ray outs the ice.
  for (int i = 0 ; i < m_icenrpointX ; ++i){
    std::vector <std::vector <point> > yzice_v;
    std::vector <std::vector <point> > ryzice_v;
    for (int j = 0 ; j < m_icenrpointY ; ++j){
      std::vector <point> zice_v;
      std::vector <point> rzice_v;
      for (int k = 0 ; k < m_icenrpointZ ; ++k){
	G4double x,y,z,ii ,ir, op, rx,ry,rz,rii,rir, rop;
	lookup >> x >>  y  >> z >> ii >> ir >> op >> rx >> ry >>rz >>rii >> rir >>rop;
 	point thepoint;
	thepoint.x = x;
	thepoint.y = y;
	thepoint.z = z;
	thepoint.thetai = ii;
	thepoint.thetar = ir;
	thepoint.optpath = op;
	zice_v.push_back(thepoint);
	//reflected point
 	point rpoint;
	rpoint.x = rx;
	rpoint.y = ry;
	rpoint.z = rz;
	rpoint.thetai = rii;
	rpoint.thetar = rir;
	rpoint.optpath = rop;
	rzice_v.push_back(rpoint);
      }
      yzice_v.push_back(zice_v);
      ryzice_v.push_back(rzice_v);
    }    
    m_icepoints.push_back(yzice_v);
    m_ricepoints.push_back(ryzice_v);
  }
  //fill the air array
  for (int i = 0 ; i < m_airnrpointX ; ++i){
    std::vector <std::vector <point> > yzair_v;
    for (int j = 0 ; j < m_airnrpointY ; ++j){
      std::vector <point> zair_v;
      for (int k = 0 ; k < m_airnrpointZ ; ++k){
	G4double ax,ay,az,aii, air, aop;
	lookup >> ax >>  ay  >> az >> aii >> air >> aop;
 	point apoint;
	apoint.x = ax;
	apoint.y = ay;
	apoint.z = az;
	apoint.thetai = aii;
	apoint.thetar = air;
	apoint.optpath = aop;
	zair_v.push_back(apoint);
	//G4cout<<i<<" "<<j<<" "<<k<<" "<<ax<<" "<<ay<<" "<<az<<" "<<aii<<" "<<air<<G4endl;
      }
      yzair_v.push_back(zair_v);
    }    
    m_airpoints.push_back(yzair_v);
  }
  
}

//look up table were filled for an horizontal ice block 
//(i.e. the x,y,z were running over -0.5->0.5/-0.15->0.15/0->0.3)
//We need first to rotate the point in these coordinate 
//before looking at the index i.e. the position of the data in the file.
//note that we don't need such operation with the air point
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
bool LookUpTable::checkicepoint(point pt){
  bool insidethebox = true;
  if (pt.x < m_icexmin || pt.x > m_icexmax || pt.y < m_iceymin || pt.y > m_iceymax || pt.z < m_icezmin || pt.z > m_icezmax)
    insidethebox = false;
  return insidethebox;
}

bool LookUpTable::checkairpoint(point pt){
  bool insidethebox = true;
  if (pt.x < m_airxmin || pt.x > m_airxmax || pt.y < m_airymin || pt.y > m_airymax || pt.z < m_airzmin || pt.z > m_airzmax)
    insidethebox = false;
  return insidethebox;
}

//after setting the new coord in horizontal position,
//we compute the index for each dimension
point LookUpTable::geticepoint(G4double xice, G4double yice, G4double zice, bool direct){
  point horipoint = gethoripoint(xice, yice, zice);
  point surfpoint ;  
  if (checkicepoint(horipoint)){
    double nrofstepx = (horipoint.x - m_icexmin - m_boundaryshift)/m_gridx;
    double nrofstepy = (horipoint.y - m_iceymin - m_boundaryshift)/m_gridy;
    double nrofstepz = (horipoint.z - m_icezmin - m_boundaryshift)/m_gridz;
    int xfloor = floor(nrofstepx);
    int xceil = ceil(nrofstepx);
    int indexx = fabs(nrofstepx - xfloor) < fabs(nrofstepx - xceil) ? xfloor:xceil ;
    if (xceil >= m_icenrpointX )
      indexx = m_icenrpointX-1;
    if (xfloor < 0 )
      indexx = 0;
    int yfloor = floor(nrofstepy);
    int yceil = ceil(nrofstepy);
    int indexy = fabs(nrofstepy - yfloor) < fabs(nrofstepy - yceil) ? yfloor:yceil;
    if (yceil >= m_icenrpointY )
      indexy = m_icenrpointY-1;
    if (yfloor < 0 )
      indexy = 0;
    int zfloor = floor(nrofstepz);
    int zceil = ceil(nrofstepz);
    int indexz = fabs(nrofstepz - zfloor) < fabs(nrofstepz - zceil) ? zfloor:zceil;
    if (zceil >= m_icenrpointZ )
      indexz = m_icenrpointZ-1;
    if (zfloor < 0)
      indexz = 0;
    
    if (direct){
      surfpoint = m_icepoints[indexx][indexy][indexz];
    }
    else
      surfpoint = m_ricepoints[indexx][indexy][indexz];
  }
  else{
    surfpoint.x = -1;
    surfpoint.y = -1;
    surfpoint.z = -1;
  } 
  return surfpoint;
}

//we compute the index for each dimension
point LookUpTable::getairpoint(G4double xair, G4double yair, G4double zair){
  point airpoint;
  airpoint.x = xair;
  airpoint.y = yair;
  airpoint.z = zair;
  point surfpoint ;  
  //G4cout<<"\t\t\t\tIn get air point: "<<xair<<" "<<yair<<" "<<zair<<G4endl;
  if (checkairpoint(airpoint)){
    //G4cout<<"\t\t\t\t Air poinmt checks out"<<G4endl;
    double nrofstepx = (airpoint.x - m_airxmin - m_boundaryshift)/m_gridx;
    double nrofstepy = (airpoint.y - m_airymin - m_boundaryshift)/m_gridy;
    double nrofstepz = (airpoint.z - m_airzmin - m_boundaryshift)/m_gridz;
    int xfloor = floor(nrofstepx);
    int xceil = ceil(nrofstepx);
    int indexx = fabs(nrofstepx - xfloor) < fabs(nrofstepx - xceil) ? xfloor:xceil ;
    if (xceil >= m_airnrpointX )
      indexx = m_airnrpointX-1;
    if (xfloor < 0 )
      indexx = 0;
    int yfloor = floor(nrofstepy);
    int yceil = ceil(nrofstepy);
    int indexy = fabs(nrofstepy - yfloor) < fabs(nrofstepy - yceil) ? yfloor:yceil;
    if (yceil >= m_airnrpointY )
      indexy = m_airnrpointY-1;
    if (yfloor < 0 )
      indexy = 0;
    int zfloor = floor(nrofstepz);
    int zceil = ceil(nrofstepz);
    int indexz = fabs(nrofstepz - zfloor) < fabs(nrofstepz - zceil) ? zfloor:zceil;
    if (zceil >= m_airnrpointZ )
      indexz = m_airnrpointZ-1;
    if (zfloor < 0)
      indexz = 0;
    surfpoint = m_airpoints[indexx][indexy][indexz];
    //G4cout<<"\t\t\t\t\tindices: "<<indexx<<" "<<indexy<<" "<<indexz<<G4endl;
    //G4cout<<"\t\t\t\t\t"<<surfpoint.x<<" "<<surfpoint.y<<" "<<surfpoint.z<<G4endl;
  }
  else{
    //G4cout<<"Air point doesn't check out..."<<G4endl;
    surfpoint.x = -1;
    surfpoint.y = -1;
    surfpoint.z = -1;
    surfpoint.thetai = -1;
    surfpoint.thetar = -1;
    surfpoint.optpath = -1;
  } 
  return surfpoint;
}

bool LookUpTable::checkpath(point thepoint){
  bool pathexist = true;
  if (thepoint.x == -1 || thepoint.y == -1)
    pathexist = false;
  return pathexist;
}

G4ThreeVector LookUpTable::getairvector(G4ThreeVector pos, bool coord){
  point surfpoint;
  point position;
  position.x = pos.x();
  position.y = pos.y();
  position.z = pos.z();
  //check first if the point is in the defined block of air
  if (checkairpoint(position)){
    surfpoint = getairpoint(pos.x(), pos.y(), pos.z());
    //check if the point can reach the antenna
    if (checkpath(surfpoint))
      if (coord)
	return  G4ThreeVector (surfpoint.x,surfpoint.y, surfpoint.z);
      else
	return  G4ThreeVector (surfpoint.thetai,surfpoint.thetar, surfpoint.optpath);
    else
      return G4ThreeVector (-1,-1,-1);
  }
  else 
    return G4ThreeVector (-1,-1,-1);
}

G4ThreeVector LookUpTable::geticevector(G4ThreeVector pos, bool coord, bool direct){
  point surfpoint;
  point position;
  position.x = pos.x();
  position.y = pos.y();
  position.z = pos.z();
  if (checkicepoint(position)){
    surfpoint = geticepoint(pos.x(), pos.y(), pos.z(), direct);
    if (checkpath(surfpoint))
      if (coord)
	return  G4ThreeVector (surfpoint.x,surfpoint.y, surfpoint.z);
      else
	return  G4ThreeVector (surfpoint.thetai,surfpoint.thetar, surfpoint.optpath);
    else
      return G4ThreeVector (-1,-1,-1);
  }
  else 
    return G4ThreeVector (-1,-1,-1);
}



void LookUpTable::printinfo(){
  std::cout << "antpos x,y,z = " << m_xant << " " << m_yant << " " << m_zant << std::endl;
  std::cout << "ice tilt = " << m_icetilt << std::endl;  
  std::cout << " ice xmin xmax = " << m_icexmin << " " << m_icexmax  << " ymin ymax = " << m_iceymin << " " << m_iceymax << " zmin zmax = " << m_icezmin << " " << m_icezmax << std::endl;
  std::cout << " air xmin xmax = " << m_airxmin << " " << m_airxmax  << " ymin ymax = " << m_airymin << " " << m_airymax << " zmin zmax = " << m_airzmin << " " << m_airzmax << std::endl;
  std::cout << "nr of points x,y,z for ice = " << m_icenrpointX << " " << m_icenrpointY << " " << m_icenrpointZ << std::endl;
  std::cout << "nr of points x,y,z for air = " << m_airnrpointX << " " << m_airnrpointY << " " << m_airnrpointZ << std::endl;
  std::cout << "grid x,y,z = " << m_gridx << " " << m_gridy << " " << m_gridz << std::endl;
}


void LookUpTable::printtable(){
  for (int i = 0 ; i < m_icenrpointX ; ++i){
    for (int j = 0 ; j < m_icenrpointY ; ++j){
      for (int k = 0 ; k < m_icenrpointZ ; ++k){
	std::cout << "x = " << m_icepoints[i][j][k].x << std::endl;
	std::cout << "y = " << m_icepoints[i][j][k].y << std::endl;
	std::cout << "z = " << m_icepoints[i][j][k].z << std::endl;
      }
    }    
  }
}


void LookUpTable::checktable(){
  double step = 0.05;
  std::ofstream outicedirect("tablecheck/icecheck.txt", std::ofstream::out);
  std::ofstream outicerefl("tablecheck/icereflectedcheck.txt", std::ofstream::out);
  std::ofstream outair("tablecheck/aircheck.txt", std::ofstream::out);
  for (double itx =  m_icexmin ; itx < m_icexmax ; itx = itx + step){
    for (double ity =  m_iceymin ; ity < m_iceymax ; ity = ity + step){
      for (double itz =  m_icezmin ; itz < m_iceymax ; itz = itz + step){
	point tiltpoint;
	tiltpoint.x = itx*cos(m_icetilt*m_pi/180) + itz*sin(m_icetilt*m_pi/180);
	tiltpoint.y = ity;
	tiltpoint.z = -itx*sin(m_icetilt*m_pi/180) + itz*cos(m_icetilt*m_pi/180);
	G4ThreeVector source = G4ThreeVector(tiltpoint.x,tiltpoint.y,tiltpoint.z);
	G4ThreeVector directpoint = geticevector(source, true, true);
	G4ThreeVector reflectedpoint = geticevector(source, true, false);
	outicedirect << tiltpoint.x << " " << tiltpoint.y << " " << tiltpoint.z <<  " " << directpoint.x() << " " << directpoint.y() << " " << directpoint.z() << std::endl;
	outicerefl << tiltpoint.x << " " << tiltpoint.y << " " << tiltpoint.z <<  " " << directpoint.x() << " " << directpoint.y() << " " << directpoint.z() << std::endl;
      }
    }    
  }

  for (double itx =  m_airxmin ; itx < m_airxmax ; itx = itx + step){
    for (double ity =  m_airymin ; ity < m_airymax ; ity = ity + step){
      for (double itz =  m_airzmin ; itz < m_airymax ; itz = itz + step){
	G4ThreeVector source = G4ThreeVector(itx,ity,itz);
	G4ThreeVector directpoint = getairvector(source, true);
	outair << source.x() << " " << source.y() << " " << source.z() <<  " " << directpoint.x() << " " << directpoint.y() << " " << directpoint.z() << std::endl;
      }
    }    
  }
  outicedirect.close();
  outicerefl.close();
  outair.close();
}
