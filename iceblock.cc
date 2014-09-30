
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Simple simulation to shoot particles into a block of some material.   //
// The goal is to study shower properties in Ice, but other materials    //
// may be used to compare with previoulsy established shower properties. //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

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
//#include "MyTreeWriter.hh"

#include "SetupAntenna.hh"

#include <sstream>
#include <fstream>
//#include <ctime>
#include "time.h"

//---------------------------------------------------//
// Help menu
//---------------------------------------------------//
void help()
{

  cout << endl;
  cout << endl;
  cout << "-------------------------------------------" << endl;
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
  cout << "-------------------------------------------" << endl;
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
  
  // Options
  G4int nEvents          = 1;     // number of events
  G4int nParticles       = 1;     // number of particles
  G4int beamEnergy       = 1000; // MeV
  G4int detMaterial      = 0;     // detector material
  std::string partType   = "e-";   
  G4double threshold     = 0;
  bool useThreshold      = false;
  std::string antFile    = "";

  // Options
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

  // Make an output name including the run information
  stringstream ss;
  ss << "output_" << nEvents << "_" << beamEnergy << "_";

  if(detMaterial == 1)      ss << "lead";
  else if(detMaterial == 2) ss << "iron";
  else                      ss << "ice";

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
  else{
    ss << "_HardCodedAntenna_R1000m";
    //ss << "_HardCodedAntenna_R1000m_oldV";
    //ss << "_HardCodedAntenna_R100m_newV";
    //ss << "_HardCodedAntenna_BetaCalc";
    //ss << "_HardCodedAntenna_R10m";
  }

  // If threshold is set, append to file
  if(useThreshold) ss << "_thresh" << threshold << "MeV";

  // Setup the Antennas
  SetupAntenna* m_AntSetup = new SetupAntenna(antFile);
  std::vector<Antenna*> m_Ants = m_AntSetup->getAnts();
  

  // My own stepping
  //G4VSteppingVerbose* verboseStep = new SteppingVerbose(ss.str());
  //G4VSteppingVerbose::SetInstance(verboseStep);

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
  std::ofstream trackOutput(("tracks/"+ss.str()+".dat").c_str(), std::ofstream::out);
  std::ofstream stepOutput(("steps/"+ss.str()+".dat").c_str(), std::ofstream::out);
  std::ofstream APotOut(("efield/A_"+ss.str()+".dat").c_str(), std::ofstream::out);

  // My output tree
  // This seems to take too long!!
  //MyTreeWriter* treeWriter = new MyTreeWriter("trees/test.root");
  //MyTreeWriter* treeWriter = new MyTreeWriter(("trees/"+TString(ss.str().c_str())+".root"));
  
  // Set Primary action to be carried out
  PrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(detector,
								 beamEnergy,
								 partType,
								 nParticles);  
  runManager->SetUserAction(genAction);
  runManager->SetUserAction(new RunAction(detector,genAction));
  runManager->SetUserAction(new EventAction(&trackOutput, &stepOutput,
					    &APotOut,
					    &m_Ants));
  runManager->SetUserAction(new SteppingAction(&stepOutput, &m_Ants));
  runManager->SetUserAction(new TrackingAction(&trackOutput));

  // Initialize G4 Kernel
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  //G4UImanager* UI = G4UImanager::GetUIpointer();
  //UI->ApplyCommand("/run/verbose 1");
  //UI->ApplyCommand("/event/verbose 1");
  //UI->ApplyCommand("/tracking/verbose 1");

  // Start the run
  runManager->BeamOn(nEvents);

  // Write tree and end
  //treeWriter->Finalize();
  //delete treeWriter;

  trackOutput.close();
  stepOutput.close();
  APotOut.close();
    

  delete runManager;
  cout<<"Number of particles"<<nParticles<<endl;
  cout<<"Runtime is: "<<(double)(clock() - tStart)/CLOCKS_PER_SEC<<" seconds "<<endl;

  delete m_AntSetup;

  return 0;

}
