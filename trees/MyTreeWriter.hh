#ifndef MyTreeWriter_hh
#define MyTreeWriter_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This is a helper class to fill the event information so that we can pass   //
// basic quantities that are output by Geant4.  There will be an Event class  // 
// which will contain all the particles.  The particle information can be     //
// updated as we figure out what we need.                                     //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "Analysis.h"

class MyTreeWriter
{

 public:
  
  // Constructor and destructor
  MyTreeWriter(TString fname);
  virtual ~MyTreeWriter();

  // Make a new event
  void CreateEvent(int evtNum);

  // Get Event
  Event* GetEvt(){ return m_event; };

  // Write Event
  void WriteEvent();

  // Write Tree and close file
  void Finalize();

 protected:
  
  TTree* m_tree;
  TFile* m_file;
  Event* m_event;  

};

#endif
