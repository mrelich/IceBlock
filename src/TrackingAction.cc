//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
// $Id: TrackingAction.cc 76464 2013-11-11 10:22:56Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

//#include "DetectorConstruction.hh"
#include "G4Track.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//TrackingAction::TrackingAction(DetectorConstruction* DET, EventAction* EA)
TrackingAction::TrackingAction(std::ofstream* file) 
  : G4UserTrackingAction(), m_output(file)
  //m_treeWriter(treeWriter)
  //fDetector(DET), fEventAction(EA)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  // few initialisations
  //
  //if (aTrack->GetTrackID() == 1) {
  //fXstartAbs = fDetector->GetxstartAbs();
  //fXendAbs   = fDetector->GetxendAbs();
  //fPrimaryCharge = aTrack->GetDefinition()->GetPDGCharge();
  //}

  G4ThreeVector position = aTrack->GetPosition();
  G4double      energy   = aTrack->GetKineticEnergy();
  G4double      total_e  = aTrack->GetTotalEnergy();
  //G4int         trk_id   = aTrack->GetTrackID();
  //G4int         parent_id= aTrack->GetParentID();
  //G4int         pdg_id   = aTrack->GetParticleDefinition()->GetPDGEncoding();

  // Save
  (*m_output) <<
    position.x()/cm << " " <<
    position.y()/cm << " " <<
    position.z()/cm << " " <<
    energy/MeV      << " " <<
    total_e/MeV      << " ";
  
  // Get the event
  //Event* event  = m_treeWriter->GetEvt();

  // Create Particle
  //event->addParticle( Particle() );
  //Particle* part = event->getLastPart();
  //part->setPDG( aTrack->GetParticleDefinition()->GetPDGEncoding() );

  // Create Track object
  //part->setTrack( Track() );
  //Track* trk = part->track();
  //trk->setInitialPosition( position.x()/cm, position.y()/cm, position.z()/cm );
  //trk->setTrackID(aTrack->GetTrackID());
  //trk->setParentID(aTrack->GetParentID());
  //trk->setInitialE(aTrack->GetKineticEnergy()/MeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4ThreeVector position = aTrack->GetPosition();
  G4int         trk_id   = aTrack->GetTrackID();
  G4int         parent_id= aTrack->GetParentID();
  G4int         pdg_id   = aTrack->GetParticleDefinition()->GetPDGEncoding();

  (*m_output) <<
    position.x()/cm << " " <<
    position.y()/cm << " " <<
    position.z()/cm << " " <<
    trk_id << " " <<
    parent_id << " " <<
    pdg_id << " " <<
    G4endl;

  // Save tracks final position
  //Track* trk = m_treeWriter->GetEvt()->getLastPart()->track();
  //trk->setFinalPosition(position.x()/cm, position.y()/cm, position.z()/cm);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

