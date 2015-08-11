

#include "DetectorConstruction.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

//-----------------------------------------------------------------//
// Constructor
//-----------------------------------------------------------------//
DetectorConstruction::DetectorConstruction(G4int detMat, G4double EThresh, 
					   bool useThresh, G4double stepLimit,
					   RefractionTool* refTool,
					   std::vector<Antenna*>* ants,
					   G4int rotation) :
  m_world_log(NULL),
  m_iceblock_log(NULL),
  m_world_phys(NULL),
  m_iceblock_phys(NULL),
  m_detMaterial(0),
  m_material(NULL),
  m_threshold(0),
  m_useThreshold(false),
  m_stepLimit(NULL),
  m_stepLimitValue(0),
  m_rotIce(NULL),
  m_refTool(NULL),
  m_ants(ants),
  m_icetilt(rotation)
{

  // Set detector material
  m_detMaterial = detMat;

  // Set threshold
  m_threshold    = EThresh;
  m_useThreshold = useThresh;
 
  // Set step limit
  m_stepLimitValue = stepLimit;

  // Tool for refraction to be initialized here
  m_refTool = refTool;

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
  //G4Element* Pb  = new G4Element("Lead"    , "Pb", z=82., a= 207.19*g/mole    );
  //G4Element* Fe  = new G4Element("Iron"    , "Fe", z=26., a= 55.845*g/mole    );

  //
  // Create Material for the detector.
  //
  
  //G4Material* det;

  if( m_detMaterial == Mat_LEAD ){
    //det =  new G4Material("LEAD", density=11.3*g/cm3, nel=1,
    //kStateSolid, 290.*kelvin); 
    //det->AddElement(Pb, 1 );
    m_material = nist->FindOrBuildMaterial("G4_Pb");
  }
  else if( m_detMaterial == Mat_IRON ){
    //det =  new G4Material("IRON", density=7.874*g/cm3, nel=1,
    //kStateSolid, 290.*kelvin); 
    //det->AddElement(Fe, 1 );
    m_material = nist->FindOrBuildMaterial("G4_Fe");    
  }
  else{ // default is ice
    m_material =  new G4Material("ICE", density=0.924*g/cm3, nel=2);
    m_material->AddElement(H, 2 );
    m_material->AddElement(O, 1 );
    m_material->GetIonisation()->SetMeanExcitationEnergy(78*eV);
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

  // Dimensions of ice at Delta
  G4double iceblock_x = 100.0 * cm;
  G4double iceblock_y = 30.0 * cm;
  G4double iceblock_z = 30.0 * cm;
  //if( !m_refTool->useTool() ){
  //iceblock_x = 1500.0 * m;
  //iceblock_y = 1500.0 * m;
  //xiceblock_z = 1500.0 * m;
  //}


  // Set block position
  G4ThreeVector block_pos = G4ThreeVector(0,0,iceblock_z*0.5);
  //if( !m_refTool->useTool() ) block_pos.set(0,0,0);

  // Define rotation matrix for block. Here rotate
  // in the y-direction only.
  G4double tilt = m_icetilt * m_pi / 180.;
  m_rotIce = new G4RotationMatrix(G4ThreeVector(cos(tilt),0,sin(tilt)),
				  G4ThreeVector(0,1,0),
				  G4ThreeVector(-sin(tilt),0,cos(tilt)));

  G4cout<<"Rotation Matrix: "<<G4endl;
  G4cout<<(*m_rotIce)<<G4endl;

  // UPDATE: So the way G4 handles the the placement of the object
  // requires that the translation vector be defined in the rotated
  // frame.  Therefore we apply the inverse rotation for the translation
  // vector to insure that the vector is with respect to the mother volume.
  block_pos *= (m_rotIce->inverse());
  G4cout<<block_pos<<G4endl;

  // Iceblock is a box
  G4Box* iceblock_box = new G4Box("ICEBLOCK",
				  0.5*iceblock_x,
				  0.5*iceblock_y,
				  0.5*iceblock_z);

  // Logical volume
  m_iceblock_log = new G4LogicalVolume(iceblock_box,
				      m_material,
				      "iceblock_log");

  // ADDING: Set the minimum cuts here
  if( m_useThreshold){
    m_iceblock_log->SetUserLimits( new G4UserLimits(DBL_MAX,           // Step length Max
						    DBL_MAX,           // Track length Max
						    DBL_MAX,           // Time Max for track
						    m_threshold * MeV, // Minimum Kinetic energy
						    0.)                // Range Min for track
				   );
  }
  
  // Physical volume
  m_iceblock_phys = new G4PVPlacement(m_rotIce,                      // Rotation Matrix
  				      block_pos,                     //iceblock_z*0.49),
  				      m_iceblock_log,                // its logical volume
  				      "iceblock_phys",               // name
  				      m_world_log,                   // Mother Volume
  				      false,                         // no boolean operator
  				      0);                            // copy number
  
  // Initialize the refraction tool with the geometry
  // settings of the ice
  //if( m_refTool->useTool() )
  m_refTool->initialize(G4ThreeVector(sin(tilt),0,cos(tilt)),
			block_pos/1000,
			G4ThreeVector(iceblock_x/1000,iceblock_y/1000,iceblock_z/1000),
			m_n,
			m_nAir,
			m_ants);
  
  //
  // Set user limits for step size
  //
  if(m_stepLimitValue > 0){
    //G4double maxStep = 0.5*mm; //1*mm; //1*cm;
    G4double maxStep = m_stepLimitValue*mm; //1*mm; //1*cm;
    m_stepLimit = new G4UserLimits(maxStep);
    m_iceblock_log->SetUserLimits(m_stepLimit);
    m_world_log->SetUserLimits(m_stepLimit);
  }
  
  //
  // Return the physical world volume
  //  

  return m_world_phys;

}
