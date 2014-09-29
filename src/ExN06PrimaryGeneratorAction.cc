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
//
// $Id: ExN06PrimaryGeneratorAction.cc,v 1.6 2006-06-29 17:54:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN06PrimaryGeneratorAction.hh"
#include "ExN06PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06PrimaryGeneratorAction::ExN06PrimaryGeneratorAction()
{
  particleGun = new G4ParticleGun();
  
  //create a messenger for this class
  gunMessenger = new ExN06PrimaryGeneratorMessenger(this);
  
  //default kinematic
  // handled in GeneratePrimaries
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06PrimaryGeneratorAction::~ExN06PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//===========================================================================

void ExN06PrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//===========================================================================

void ExN06PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }
     	       
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton); 
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}

//===========================================================================

void ExN06PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //particleTable ->DumpTable("ALL");


  particleGun->SetParticleTime(0.0*ns);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  G4int ev_type = 4;  // 3=parallel rays (turn off agel).   4-ring (turn agel on).
  G4double Ekin;

  if (ev_type == 1) {            // 2 muons, 1 GeV  ======================================================================
    particle = particleTable->FindParticle("mu+");
	particleGun->SetParticleDefinition(particle);

	G4double momentum = 1000.*MeV;
    G4double sigmaMomentum = 0.*MeV;
    G4double sigmaAngle = 0.*deg;

    G4double pp = momentum + (G4UniformRand()-0.5)*sigmaMomentum;
    G4double mass = particle->GetPDGMass();
    Ekin = std::sqrt(pp*pp+mass*mass)-mass;
    particleGun->SetParticleEnergy(Ekin);

    G4double angle = (G4UniformRand()-0.5)*sigmaAngle;
    particleGun->SetParticleMomentumDirection(G4ThreeVector(std::sin(angle),0.,std::cos(angle)));

    G4int nparticles=2;
    G4double sigmayy = 6.0*cm;
    G4double yy;
    for (G4int i=1; i<=nparticles; i++) {
      yy = 2.0*cm + (G4UniformRand()-0.5)*sigmayy;
      particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, yy ,-11.0*cm));
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }    // ev_type=1

  if (ev_type == 2) {     // x,y=0, sweep over wavelengths 200-600 nm  ====================================================
		G4int nbins=10000;    //                                                  study spectra
		particle = particleTable->FindParticle("opticalphoton");
		particleGun->SetParticleDefinition(particle);
	    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1.0));
	    particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm ,-11.0*cm));
		for (G4int i=1; i<nbins; i++) {
			Ekin = ( 1.7+(4.8-1.7)*G4UniformRand() )*eV;   // 1.7-4.8eV ~= 200-600 nm
			particleGun->SetParticleEnergy(Ekin);
		    particleGun->GeneratePrimaryVertex(anEvent);
		}   // end loop over wavelengths at fixed position
  }     // end ev_type=2

  if (ev_type == 3) {     // fixed wavelength, x=0, sweep over y =========================================================
	G4int nbins=100;    //                                                  study focusing
	particle = particleTable->FindParticle("opticalphoton");
	particleGun->SetParticleDefinition(particle);
//    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1.0));
//    particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0.1737, 0.9848));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    Ekin = 3.2*eV;   // ~400 nm
	particleGun->SetParticleEnergy(Ekin);

	for (G4int i=1; i<nbins; i++) {
		G4double angle = G4UniformRand()*360.0*deg;
		SetOptPhotonPolar(angle);
        G4double yrand = (G4UniformRand()-0.5)*8.0*cm;
		particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, yrand ,-9.5*cm));
	    particleGun->GeneratePrimaryVertex(anEvent);
	}   // end loop over wavelengths at fixed position
  }

  if (ev_type == 4) {            // 10 muons, 5 GeV at x,y = (0,0) ===================================================================
    particle = particleTable->FindParticle("mu+");
	particleGun->SetParticleDefinition(particle);

	G4double momentum = 5000.*MeV;
	G4double sigmaMomentum = 0.*MeV;
	G4double sigmaAngle = 0.*deg;
	G4double pp = momentum + (G4UniformRand()-0.5)*sigmaMomentum;
    G4double mass = particle->GetPDGMass();
    Ekin = std::sqrt(pp*pp+mass*mass)-mass;
    particleGun->SetParticleEnergy(Ekin);
	    G4double angle = (G4UniformRand()-0.5)*sigmaAngle;
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1.0));     // straight
	//particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0.139,0.990));    // 8 degrees
    //particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0.208,0.978));    // 12 degrees
	//particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0.309,0.951));    // 18 degrees
    G4int npart = 1;
    for (G4int i=1; i<=npart; i++) {
      if (i<=npart)particleGun->SetParticlePosition(G4ThreeVector(0.01*cm, 3.01*cm ,-11.0*cm));
      if (i>npart) particleGun->SetParticlePosition(G4ThreeVector(0.01*cm, 5.0*cm ,-11.0*cm));
      particleGun->GeneratePrimaryVertex(anEvent);
    }

  }    // ev_type=4

}
