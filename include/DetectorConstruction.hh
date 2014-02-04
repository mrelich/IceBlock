#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4NistManager.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{

 public:
  
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  
 private:
  
  // Logical volumes
  G4LogicalVolume* m_world_log;
  G4LogicalVolume* m_iceblock_log;
  
  // Physical volumes
  G4VPhysicalVolume* m_world_phys;
  G4VPhysicalVolume* m_iceblock_phys;
  
};

#endif
