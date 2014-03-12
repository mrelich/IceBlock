
#include "MyTreeWriter.hh"

//------------------------------------------------------------//
// Constructor
//------------------------------------------------------------//
MyTreeWriter::MyTreeWriter(TString fname)
{

  // Open a file
  m_file = new TFile(fname.Data(),"recreate");

  // Associate the tree with the file
  m_tree = new TTree("geant4","Geant4 Output Tree");
  m_tree->SetAutoSave(1000000);
  m_tree->SetMaxTreeSize(3000000000u);
  
  // Create an event object
  m_event = new Event();

  // Associate event to tree
  m_tree->Branch("event", &m_event);
  
}

//------------------------------------------------------------//
// Destructor
//------------------------------------------------------------//
MyTreeWriter::~MyTreeWriter()
{

  delete m_event;

}

//------------------------------------------------------------//
// Create an event
//------------------------------------------------------------//
void MyTreeWriter::CreateEvent(int evtNum)
{
  
  m_event->Clear();
  m_event->setEvtNumber(evtNum);

}

//------------------------------------------------------------//
// Write Event
//------------------------------------------------------------//
void MyTreeWriter::WriteEvent()
{

  m_tree->Fill();

}

//------------------------------------------------------------//
// Finalize
//------------------------------------------------------------//
void MyTreeWriter::Finalize()
{

  m_file->Write();
  m_file->Close();

}
