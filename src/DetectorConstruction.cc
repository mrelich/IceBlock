
#include "DetectorConstruction.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
DetectorConstruction::DetectorConstruction(G4int detMat) :
  m_world_log(NULL),
  m_iceblock_log(NULL),
  m_world_phys(NULL),
  m_iceblock_phys(NULL),
  m_detMaterial(0)
{

  // Set detector material
  m_detMaterial = detMat;

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
  G4Element* H   = new G4Element("Hydrogen", "H" , z=1. , a=1.00794*g/mole );
  G4Element* O   = new G4Element("Oxygen"  , "O" , z=8. , a= 16.00*g/mole  );
  G4Element* Pb  = new G4Element("Lead"    , "Pb", z=82., a= 207.19*g/mole    );
  G4Element* Fe  = new G4Element("Iron"    , "Fe", z=26., a= 55.845*g/mole    );
  
  //
  // Create Material for the detector.
  //

  G4Material* det;
  if( m_detMaterial == Mat_LEAD ){
    det =  new G4Material("LEAD", density=11.3*g/cm3, nel=1,
			  kStateSolid, 290.*kelvin); 
    det->AddElement(Pb, 1 );
    
  }
  else if( m_detMaterial == Mat_IRON ){
    det =  new G4Material("IRON", density=7.874*g/cm3, nel=1,
			  kStateSolid, 290.*kelvin); 
    det->AddElement(Fe, 1 );
    
  }
  else{ // default is ice
    det =  new G4Material("ICE", density=0.920*g/cm3, nel=2,
			  kStateSolid, 216.15*kelvin); 
    det->AddElement(H, 2 );
    det->AddElement(O, 1 );

  }

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
  
  G4double iceblock_x = 1500.0 * m;
  G4double iceblock_y = 1500.0 * m;
  G4double iceblock_z = 1500.0 * m;

  // Iceblock is a box
  G4Box* iceblock_box = new G4Box("ICEBLOCK",
				  0.5*iceblock_x,
				  0.5*iceblock_y,
				  0.5*iceblock_z);

  // Logical volume
  m_iceblock_log = new G4LogicalVolume(iceblock_box,
				      det,
				      "iceblock_log");

  // Try adding User limits.
  // UPDATE: Doesn't seem to impact results at all...
  //m_iceblock_log->SetUserLimits(new G4UserLimits(DBL_MAX,   // stepMax
  //DBL_MAX,   // trackMax
  //DBL_MAX,   // timeMax
  //0.611*MeV, // kinMin
  ////100.0*MeV, // kinMin
  //0)         // rangeMin
  //);

  // Physical volume
  m_iceblock_phys = new G4PVPlacement(0,                             // no rotation
				      G4ThreeVector(),
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
