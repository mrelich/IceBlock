
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Full Simulation using Geant4 to calculate the electric field for an //
// EM shower in Ice.  It is to be used for the ARA TA-ELS experiment   //
// conducted in Utah.                                                  // 
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//


// Geant4 Classes
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "TrackingAction.hh"

// Some additional classes for convenience
#include "BeamProfile.hh"
#include "SetupAntenna.hh"

// Standard
#include <sstream>
#include <fstream>
#include "time.h"

//---------------------------------------------------//
// Help menu for all the options
//---------------------------------------------------//
void help()
{

  cout << endl;
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "-ne <int> " << endl;
  cout << "\t Specify the number of events " << endl;
  cout << "\t default is 1" << endl;
  cout << "-np <int> " << endl;
  cout << "\t Specify the number of particles " << endl;
  cout << "\t default is 1 " << endl;
  cout << "-e <int> " << endl;
  cout << "\t Specify the beam energy in MeV" << endl;
  cout << "\t Default is 1000 MeV" << endl;
  cout << "-t <int> " << endl; 
  cout << "\t Specify the interaction material " << endl;
  cout << "\t Default is 0 (ice) " << endl;
  cout << "\t\t 0 -- ice " << endl;
  cout << "\t\t 1 -- lead " << endl;
  cout << "\t\t 2 -- iron " << endl;
  cout << "-p <char> " << endl;
  cout << "\t Specify the particle type" << endl;
  cout << "\t Default is e-" << endl;
  cout << "\t\t e- -- Electron" << endl;
  cout << "\t\t e+ -- Positron" << endl;
  cout << "\t\t mu- -- Muon" << endl;
  cout << "\t\t gamma -- Photon" << endl;
  cout << "-c <float>" << endl;
  cout << "\t Specify the energy threshold to use in MeV" << endl;
  cout << "\t Default is 0" << endl;
  cout << "-i <file.txt> " << endl;
  cout << "\t Defaule is NULL" << endl;
  cout << "\t\tSpecify input file to read antenna config from"<<endl;
  cout << "--flat <float>" << endl;
  cout << "\t Inject random flat distribution with width <float> [mm]" << endl;
  cout << "--gauss <float>" << endl;
  cout << "\t Inject Gaussian distribution with <float> sigma [mm]" << endl;
  cout << "--nbunch <int>" << endl;
  cout << "\t Create N bunches" << endl;
  cout << "--offset <float>" << endl;
  cout << "\t Timing offset for bunches (default 0.35 ns)"<<endl;
  cout << "--beam <file.txt> " << endl;
  cout << "\t Input beam profile file" << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << endl;
  cout << endl;

}

//---------------------------------------------------//
// Main
//---------------------------------------------------//
int main(int argc, char** argv)
{

  // Store time
  clock_t tStart = clock();

  //============================================//
  // List all the options here with some small
  // description of what they are used for
  //============================================//
  G4int nEvents          = 1;     // number of events
  G4int nParticles       = 1;     // number of particles
  G4int beamEnergy       = 1000;  // beam energy
  G4int detMaterial      = 0;     // detector material
  std::string partType   = "e-";  // primary particle
  G4double threshold     = 0;     // energy threshold for testing -- leave 0
  bool useThreshold      = false; // Whether or not to use this threshold -- don't
  std::string antFile    = "";    // File path for antenna positions
  G4bool b_flat          = false; // Use uniform beam
  G4bool b_gauss         = false; // Use Gaussian distributed beam
  G4double sigma         = 0;     // width of the beam
  G4int nbunches         = 1;     // Number of bunches in beam
  G4double tOffset       = 0.350; // Timing between bunches [ns]
  std::string beamFile   = "";    // Beam profile from txt file

  //============================================//
  // Load the options
  //============================================//
  for(int i=1; i<argc; ++i){
    if( strcmp(argv[i], "-ne") == 0 )
      nEvents = atoi( argv[++i] );
    else if( strcmp(argv[i], "-np") == 0 )
      nParticles = atoi( argv[++i] );
    else if( strcmp(argv[i], "-e") == 0 )
      beamEnergy = atoi( argv[++i] );
    else if( strcmp(argv[i], "-t") == 0 )
      detMaterial = atoi( argv[++i] );
    else if( strcmp(argv[i], "-p") == 0 )
      partType = string(argv[++i]);
    else if( strcmp(argv[i], "-c") == 0 ){
      threshold = atof( argv[++i] );
      useThreshold = true;
    }
    else if( strcmp(argv[i], "-i") == 0 )
      antFile = argv[++i];
    else if( strcmp(argv[i], "--flat") == 0 ){
      b_flat = true;
      sigma = atof( argv[++i] );
    }
    else if( strcmp(argv[i], "--gauss") == 0 ){
      b_gauss = true;
      sigma = atof( argv[++i] );
    }
    else if( strcmp(argv[i], "--nbunch") == 0 )
      nbunches = atoi( argv[++i] );
    else if( strcmp(argv[i], "--offset") == 0 )
      tOffset = atof( argv[++i] );
    else if( strcmp(argv[i], "--beam") == 0 )
      beamFile = argv[++i];
    else{
      help();
      return 0;
    }
  }//end loop over arguments

  // Make sure det material makes sense
  if( detMaterial > 2 ){
    cout<<"Error: Det material not recognized: "<<detMaterial<<endl;
    help();
    return 0;
  }

  //============================================//
  // Make an output name including the run 
  // information. Perhaps could come up with 
  // something better.
  // THIS IS GETTING TOO LONG
  //============================================//
  stringstream ss;
  ss << "output_" << nEvents << "_" << beamEnergy << "_";

  // Check det material
  if(detMaterial == 1)      ss << "lead";
  else if(detMaterial == 2) ss << "iron";
  else                      ss << "ice";

  // Check particle type
  if(partType == "e-")         ss << "_eBeam";    
  else if(partType == "mu-")   ss << "_muBeam";    
  else if(partType == "e+")    ss << "_pBeam";    
  else if(partType == "gamma") ss << "_gBeam";
  else                         ss << "_unkBeam";
  
  // Add the number of particles
  ss << "_np" << nParticles;

  // Decide if antenna file specified
  if( !antFile.empty() ){
    string s_ant = antFile.substr(0,antFile.find(".txt"));
    s_ant        = s_ant.substr(antFile.find("/")+1,-1);
    cout<<s_ant<<endl;
    ss << "_" << s_ant;
  }
  else{ // uses a single antenna a theta_c 100m away
    ss << "_HardCodedAntenna_R100m";
  }

  // Add some more info about input distribution
  BeamProfile* bp = new BeamProfile();
  if(b_flat)  ss << "_RandFlat" << sigma;
  else if(b_gauss) ss << "_RandGauss" << sigma;
  else             ss << "_singlePos";

  // If threshold is set, append to file
  if(useThreshold) ss << "_thresh" << threshold << "MeV";

  // Add number of bunches information
  if(beamFile.empty()) ss << "_bunches" << nbunches;
  else{
    bp->init(beamFile);
    if( !bp->isInit() ){
      cout<<"************************************"<<endl;
      cout<<"Error: Beam Profile not initialized"<<endl;
      cout<<"************************************"<<endl;
      return 0.;
    }
    string s_beam = beamFile.substr(0,beamFile.find(".txt"));
    s_beam        = s_beam.substr(beamFile.find("/")+1,-1);
    cout<<"Loading beam: "<<s_beam<<endl;
    ss << "_" << s_beam;
    //nbunches = bp->getN(); 
  }

  //============================================//
  // Configure Geant Below
  //============================================//

  // Setup the Antennas
  SetupAntenna* m_AntSetup = new SetupAntenna(antFile);
  std::vector<Antenna*> m_Ants = m_AntSetup->getAnts();
  
  // Default run manager
  G4RunManager* runManager = new G4RunManager;

  // Construct detector
  DetectorConstruction* detector = new DetectorConstruction(detMaterial,
							    threshold,
							    useThreshold);
  runManager->SetUserInitialization(detector);

  // Set Physics list
  G4VUserPhysicsList* physics = new PhysicsList(useThreshold);
  runManager->SetUserInitialization(physics);
  
  // Open up a file for output
  // This is not currently used anymore since we have the
  // e-field directly calculated
  //std::ofstream trackOutput(("tracks/"+ss.str()+".dat").c_str(), std::ofstream::out);
  //std::ofstream stepOutput(("steps/"+ss.str()+".dat").c_str(), std::ofstream::out);
  std::ofstream APotOut(("efield/A_"+ss.str()+".dat").c_str(), std::ofstream::out);

  // Output some meta data info for the vector potential plots
  APotOut << "# " << nEvents << " " 
	  << nParticles << " "
	  << beamEnergy << " "
	  << m_Ants.size() << " "
	  << m_Ants.at(0)->getNPoints() << " "
	  << m_Ants.at(0)->getTStep()   << " "
	  << nbunches                   << " "
	  << tOffset
	  << G4endl;

  // Set Primary action generator which will inject particles
  PrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(detector,
								 beamEnergy,
								 partType,
								 nParticles,
								 b_flat, b_gauss, sigma,
								 nbunches, tOffset,
								 bp);  
  
  runManager->SetUserAction(genAction);

  // Setup the actions to carry out  
  runManager->SetUserAction(new RunAction(detector,genAction));
  runManager->SetUserAction(new EventAction(NULL,NULL, //&trackOutput, &stepOutput,
					    &APotOut,
					    &m_Ants));
  runManager->SetUserAction(new SteppingAction(NULL /*&stepOutput*/, &m_Ants));
  //runManager->SetUserAction(new TrackingAction(&trackOutput));

  // Initialize G4 Kernel
  runManager->Initialize();

  // Start the run
  runManager->BeamOn(nEvents);

  // Close up the files
  //trackOutput.close();
  //stepOutput.close();
  APotOut.close();
    
  // Dump some info
  cout<<"Number of particles "<<nParticles<<endl;
  cout<<"Runtime is: "<<(double)(clock() - tStart)/CLOCKS_PER_SEC<<" seconds "<<endl;

  // Clean up and finish
  delete m_AntSetup;
  delete runManager;
  return 0;

}
