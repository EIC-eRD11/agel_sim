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
// Rich advanced example for Geant4
// RichTbSteppingAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "globals.hh"
#include "ExN06SteppingAction.hh"
#include "G4SteppingManager.hh"
//#include "RichTbAnalysisManager.hh"
//#include "RichTbMaterial.hh"
//#include "RichTbGeometryParameters.hh"
//#include "RichTbMaterialParameters.hh"
//#include "RichTbRunConfig.hh"
//#include "RichTbPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4ThreeVector.hh"
#include "G4OpticalPhoton.hh"
#include "G4PionMinus.hh"

#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;


#ifdef G4ANALYSIS_USE
  #include "AIDA/AIDA.h"
#endif

ExN06SteppingAction::
ExN06SteppingAction(  ExN06PrimaryGeneratorAction* RPrimGenAction)
{
//  exn06RunConfig= rConfig;
  rPrimGenAction = RPrimGenAction;
//  HpdPhElectronKE=rConfig->getHpdPhElectronEnergy();
  uParticleChange=new G4VParticleChange();  

  G4String fout = "photons.txt";
  fOut.open(fout.c_str(), std::ios_base::out);

}

ExN06SteppingAction::~ExN06SteppingAction()
{
  fOut.close();
}

void ExN06SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	//cout <<"Entry"<<endl;

	G4StepPoint* pPostStepPoint = aStep ->GetPostStepPoint();
	G4StepPoint* pPreStepPoint  = aStep ->GetPreStepPoint();

    if ( pPostStepPoint->GetPhysicalVolume() != NULL) {
	  assert (pPostStepPoint->GetPhysicalVolume());

      G4String tpostVol = pPostStepPoint->GetPhysicalVolume()->GetName();
      G4double ke = pPostStepPoint->GetKineticEnergy();    // in eV
      G4double wavel = ke*1240*1e5;                        // wavelength in nm
      if (tpostVol == "phodet") {
    	G4double pcharge = pPostStepPoint->GetCharge();
    	G4double g_x = pPostStepPoint->GetPosition().getX();
	    G4double g_y = pPostStepPoint->GetPosition().getY();
	    G4double g_z = pPostStepPoint->GetPosition().getZ();


	    G4ThreeVector g_mom = pPostStepPoint->GetMomentumDirection();
        fOut <<iev<<" "<< g_x <<" "<< g_y <<" "<< g_z << " "<<pcharge<<" "<<ke<<" "<<wavel<<" "<< g_mom[0] <<" "<< g_mom[1] <<" "<< g_mom[2] << endl;
    	if (pcharge != 0.0) {iev++;}
      }
    }
  //cout <<"Out"<<endl;
}


