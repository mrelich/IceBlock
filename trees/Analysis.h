#ifndef Analysis_h
#define Analysis_h

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Object class for analysis objects.  This will be for   //
// the event, particles, and particle properties.         //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "TObject.h"
#include <vector>
#include "TVector3.h"

//--------------------------------//
// Class for step information
//--------------------------------//
class Step : public TObject
{

 public:

  // Constructor and destructor
  Step(float x, float y, float z, 
       float dE, float dX, 
       int eLossProc){

    // Position
    m_position = TVector3(x, y, z);

    m_dE        = dE;         // Energy lost in step
    m_dX        = dX;         // Step size
    m_eLossProc = eLossProc;  // process caused E-loss
  };
  Step(){};
  virtual ~Step(){ clear(); };

  // Methods to retrive step information
  TVector3  getPos(){ return m_position; };
  float     getDE() { return m_dE; };
  float     getDX() { return m_dX; };
  int getELossProc(){ return m_eLossProc; };

  // CLean up
  void clear(){
    m_position.Clear();
    m_dE = 0;
    m_dX = 0;
    m_eLossProc = -999;
  };

  ClassDef(Step,1);

 protected:
  
  TVector3 m_position;
  float    m_dE;
  float    m_dX;
  int      m_eLossProc;

};

//-----------------------------------//
// Track Class
//-----------------------------------//
class Track : public TObject
{

 public:

  // Constructor and destructor
  Track(){};
  virtual ~Track(){
    clear();
  };

  // Assignement operator
  Track& operator=(const Track &rhs){
    if( this != &rhs ){
      m_initPos = rhs.m_initPos;
      m_finalPos = rhs.m_finalPos;
      m_E = rhs.m_E;
      m_trkID = rhs.m_trkID;
      m_parentID = rhs.m_parentID;
    }
    return *this;
  };
  
  // Add properties
  void setInitialPosition(float x, float y, float z){ m_initPos = TVector3(x,y,z); };
  void setFinalPosition(float x, float y, float z){ m_finalPos = TVector3(x,y,z); };
  void setInitialE(float energy){ m_E = energy; };
  void setTrackID(int trkid){ m_trkID = trkid; };
  void setParentID(int id){ m_parentID = id; };
  
  // Get properties
  TVector3 initPos(){ return m_initPos; };
  TVector3 finalPos(){ return m_finalPos; };
  float E(){ return m_E; };
  int parentID(){ return m_parentID; };
  int trkID(){ return m_trkID; };

  // Clear
  void clear(){
    m_initPos.Clear();
    m_finalPos.Clear();
    m_E = 0;
    m_parentID = -999;
    m_trkID = -999;
  };

  ClassDef(Track,1);

 protected:
  
  TVector3 m_initPos;
  TVector3 m_finalPos;
  float    m_E;
  int      m_trkID;
  int      m_parentID;
  

};

//-----------------------------------//
// Particle Class
//-----------------------------------//
class Particle : public TObject
{

 public:
  
  // Constructor and Destructor
  Particle(){};
  virtual ~Particle(){
    clear();
  };

  // Assignment operator
  Particle& operator=(const Particle &rhs){
    if( this != &rhs ){
      m_pdg = rhs.m_pdg;
      m_trk = rhs.m_trk;
      m_steps = rhs.m_steps;
    }
    return *this;
  };
  
  // Add properties
  void setPDG(int pdgid){ m_pdg = pdgid; };
  void setTrack(Track trk){ m_trk = trk; };
  void addStep(Step step){ m_steps.push_back(step); };

  // Retrieve properties
  int pdg(){ return m_pdg; };
  Track* track(){ return &m_trk; };
  std::vector<Step> getSteps(){ return m_steps; };
  Step* getStep(int pos){ return &(m_steps.at(pos)); };
  Step* getLastStep(){ return &(m_steps.back()); };

  // Clear
  void clear(){
    m_pdg = 0;
    m_trk.clear();
    m_steps.clear();
  };
  
  ClassDef(Particle,1);
  
 protected:
  
  int m_pdg;
  Track m_trk;
  std::vector<Step> m_steps;

};

//-----------------------------------//
// Event Class
//-----------------------------------//
class Event : public TObject
{

 public:
  
  // Constructor and Destructor
  Event(){};
  virtual ~Event(){
    clear();
  };
  
  // Set event properties
  void setEvtNumber(int evtNum){ m_evtNum = evtNum; };
  void addParticle(Particle part){ m_particles.push_back(part); };

  // Get event properties
  int getEvtNumber(){ return m_evtNum; };
  std::vector<Particle> getParts(){ return m_particles; };
  Particle* getPart(int pos){ return &(m_particles.at(pos)); };
  Particle* getLastPart(){ return &(m_particles.back()); };

  // Reset Properites
  void clear(){
    m_evtNum = -999;
    m_particles.clear();
  };

  ClassDef(Event,1);

 protected:
  
  int m_evtNum;                 // The event number
  std::vector<Particle> m_particles; // Store particles for event

};

#endif
