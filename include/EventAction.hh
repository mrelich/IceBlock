#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include <fstream>
//#include "MyTreeWriter.hh"
#include "Antenna.hh"

class EventAction : public G4UserEventAction
{

 public:
  
  // Constructor/Destructor
  EventAction(std::ofstream* trkFile, std::ofstream* stepFile,
	      std::ofstream* eFile, std::vector<Antenna*>* ants);
  ~EventAction();

  // Required methods
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  std::ofstream* m_trkOutput;        // File to store the trk output
  std::ofstream* m_stepOutput;       // File to store step output
  //MyTreeWriter* m_treeWriter; 
  std::ofstream* m_EfieldOut;         // File to store E-field vars
  std::vector<Antenna*>* m_ants;      // Pointer to antenna

};

#endif
