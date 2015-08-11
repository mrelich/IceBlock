#ifndef LookUpTable_h
#define LookUpTable_h

#include "globals.hh"
#include "G4Track.hh"
#include <vector>

//structure included in the array used to look up
//thetai is the incident angle
//thetar is either the refracted or reflected angle
struct point{
  G4double x;
  G4double y;
  G4double z;
  G4double thetai;
  G4double thetar;
  G4double optpath;
};

class LookUpTable
{
public:

  // Constructor
  //  LookUpTable(); //grid for each dimension in m
  LookUpTable(std::string lookupfile); 
  // Destructor
  ~LookUpTable();

  void FillInfo(std::ifstream& file);
  void FillData(std::ifstream& file);
  void FillDataAirIce(std::ifstream& file);

  void FillTestData(std::ifstream& file);
  void printtable();
  void printinfo();
  //function to set back in the coordinate of the flat ice
  //need to compute the index of the table
  point gethoripoint(G4double x, G4double y, G4double z);
  G4ThreeVector getG4horipoint(G4double x, G4double y, G4double z);
  bool checkicepoint(point p);
  bool checkairpoint(point p);
  bool checkpath(point p);
  
  // Methods to access some variables
  point geticepoint(G4double x, G4double y, G4double z, bool direct);
  point getairpoint(G4double x, G4double y, G4double z);
  //  G4ThreeVector getG4vector(G4ThreeVector position);
  G4ThreeVector getairvector(G4ThreeVector position, bool coord);
  //directpath bool used to specify if you want the direct path or the reflected one
  //coord bool is to specify if you want the coordinate vector or the angle and total path vector
  G4ThreeVector geticevector(G4ThreeVector position, bool coord, bool directpath);
  void checktable();
  
  // CLear out vector potential info
  void clear(){
    m_icepoints.clear();
    m_ricepoints.clear();
    m_airpoints.clear();
  };
  
private:
  //////common variables for air and ice point/////
  G4double m_xant;
  G4double m_yant;
  G4double m_zant;

  G4double m_gridx;
  G4double m_gridy;
  G4double m_gridz;

  G4double m_boundaryshift;

  G4double m_icetilt;
  
  ///////variables for ice/////
  //ice block dimension
  G4double m_icexmin;
  G4double m_icexmax;
  G4double m_iceymin;
  G4double m_iceymax;
  G4double m_icezmin;
  G4double m_icezmax;
  //nr of point for ice table
  int m_icenrpointX;
  int m_icenrpointY;
  int m_icenrpointZ;
  //nr of reflection point for ice table
  int m_ricenrpointX;
  int m_ricenrpointY;
  int m_ricenrpointZ;

  ///////variables for ice/////
  //air block dimension
  G4double m_airxmin;
  G4double m_airxmax;
  G4double m_airymin;
  G4double m_airymax;
  G4double m_airzmin;
  G4double m_airzmax;
  //nr of points for the air table
  int m_airnrpointX;
  int m_airnrpointY;
  int m_airnrpointZ;


  // point in ice
  std::vector < std::vector< std::vector<point> > > m_icepoints;
  // reflected point in ice
  std::vector < std::vector< std::vector<point> > > m_ricepoints;
  // point in air
  std::vector < std::vector< std::vector<point> > > m_airpoints;
  //for test purpose
  //  std::vector < std::vector< std::vector<point> > > m_sourcepoints;
  
  

};

#endif
