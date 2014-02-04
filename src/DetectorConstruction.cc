
#include "DetectorConstruction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
DetectorConstruction::DetectorConstruction() :
  m_world_log(NULL),
  m_iceblock_log(NULL),
  m_world_phys(NULL),
  m_iceblock_phys(NULL)
{

}

//-----------------------------------------------------------------//
// Destructor
//-----------------------------------------------------------------//
DetectorConstruction::~DetectorConstruction()
{

}

//-----------------------------------------------------------------//
// Construct the volumes
//-----------------------------------------------------------------//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //
  // Get NIST manager and load AIR for world
  //

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");

  // For keepting track of what is initialized
  G4double a, z;
  G4double density;
  G4int nel;

  // Add some elements
  G4Element* H   = new G4Element("Hydrogen", "H", z=1., a=1.00794*g/mole );
  G4Element* O   = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole  );

  //
  // Create Material for the detector.
  //

  G4Material* ICE =  new G4Material("ICE", density=0.920*g/cm3, nel=2,
  kStateSolid, 216.15*kelvin); 
  ICE->AddElement(H, 2 );
  ICE->AddElement(O, 1 );

  //
  // Create the world volume
  //
  
  G4double world_x = 2000.0 * m;
  G4double world_y = 2000.0 * m;
  G4double world_z = 2000.0 * m;
  
  // World box
  G4Box* world_box = new G4Box("WORLD", 
			       0.5*world_x, 
			       0.5*world_y, 
			       0.5*world_z);
  // Logical volume
  m_world_log = new G4LogicalVolume(world_box, 
				    Air, 
				    "world_log");
  
  // Physical volume
  m_world_phys = new G4PVPlacement(0,                // no rotation
				   G4ThreeVector(),  // at (0,0,0)
				   m_world_log,      // the logical volume  
				   "world_phys",     // name
				   0,                // Mother Volume 
				   false,            // no boolean operator
				   0);               // copy number
  
  //
  // Create the giant iceblock
  //
  
  G4double iceblock_x = 1000.0 * m;
  G4double iceblock_y = 1000.0 * m;
  G4double iceblock_z = 1000.0 * m;

  // Iceblock is a box
  G4Box* iceblock_box = new G4Box("ICEBLOCK",
				  0.5*iceblock_x,
				  0.5*iceblock_y,
				  0.5*iceblock_z);

  // Logical volume
  m_iceblock_log = new G4LogicalVolume(iceblock_box,
				      ICE,
				      "iceblock_log");
 

  // Physical volume
  m_iceblock_phys = new G4PVPlacement(0,                             // no rotation
				      G4ThreeVector(0,0,1000.0 * m), // shifted z->1km
				      m_iceblock_log,                 // its logical volume
				      "iceblock_phys",                // name
				      m_world_log,                   // Mother Volume
				      false,                         // no boolean operator
				      0);                            // copy number
  
  
  //
  // Return the physical world volume
  //

  return m_world_phys;

}
