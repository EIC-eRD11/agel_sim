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
// RichTbSteppingAction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef ExN06SteppingAction_h
#define ExN06SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include <iostream>

class ExN06AnalysisManager;
//class ExN06RunConfig;
class ExN06PrimaryGeneratorAction;

class ExN06SteppingAction : public G4UserSteppingAction {
	G4int iev;

 public:

  ExN06SteppingAction( ExN06PrimaryGeneratorAction* );
  virtual ~ExN06SteppingAction();
  void UserSteppingAction(const G4Step* aStep);
  ExN06AnalysisManager* getRAnalysisManager()
  {return  ranalysisManager; }
//  ExN06RunConfig* getExN06RunConfig()
//  {return exn06RunConfig; }
  G4double getHpdPhElectronKE() 
  {return HpdPhElectronKE; }
  
private:
  std::fstream fOut;

  ExN06AnalysisManager* ranalysisManager;
//  ExN06RunConfig* exn06RunConfig;
  ExN06PrimaryGeneratorAction* rPrimGenAction;
  G4double HpdPhElectronKE;
  G4VParticleChange* uParticleChange;

};
#endif
