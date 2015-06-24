#ifndef LookUpTable_h
#define LookUpTable_h

#include "globals.hh"
#include "G4Track.hh"
#include <vector>

struct point{
  G4double x;
  G4double y;
  G4double z;
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

  void FillTestData(std::ifstream& file);
  void printtable();
  void printinfo();
  point gethoripoint(G4double x, G4double y, G4double z);
  G4ThreeVector getG4horipoint(G4double x, G4double y, G4double z);
  bool checkpoint(point p);
  
  // Methods to access some variables
  point getpoint(G4double x, G4double y, G4double z);
  G4ThreeVector getG4vector(G4ThreeVector position);
  G4double getX(G4double x, G4double y, G4double z);
  G4double getY(G4double x, G4double y, G4double z);
  G4double getZ(G4double x, G4double y, G4double z);

  // CLear out vector potential info
  void clear();
  
private:
  G4double m_xant;
  G4double m_yant;
  G4double m_zant;

  G4double m_gridx;
  G4double m_gridy;
  G4double m_gridz;

  G4double m_icetilt;
  
  G4double m_xmin;
  G4double m_xmax;
  G4double m_ymin;
  G4double m_ymax;
  G4double m_zmin;
  G4double m_zmax;

  int m_nrpointX;
  int m_nrpointY;
  int m_nrpointZ;

    

  std::vector < std::vector< std::vector<point> > > m_points;
  
  

};

#endif
