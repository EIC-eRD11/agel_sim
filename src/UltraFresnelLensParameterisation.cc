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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraFresnelLensParameterisation.cc
//    ****************************************************
//
//    Class derived from G4VPVParameterisation and used to define a Fresnel lens geometry
//    through a parameterised replication of G4Cons volumes. These volumes are frustra 
//    of cones describing the lens grooves.
//    An  UltraFresnelLensParameterisation object is created in the UltraFresnelLens class
//
#include <cmath>
#include "UltraFresnelLensParameterisation.hh"
#include "UltraFresnelLens.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraFresnelLensParameterisation::UltraFresnelLensParameterisation(UltraFresnelLens* Lens)  
{
   
   FresnelLens     = Lens ;
   GrooveWidth     = Lens->GetGrooveWidth() ;
   NumberOfGrooves = Lens->GetNumberOfGrooves() ;

   dZOffset = Lens->GetSagita((NumberOfGrooves-0)*(GrooveWidth)) - 
              Lens->GetSagita((NumberOfGrooves-1)*(GrooveWidth)) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraFresnelLensParameterisation::~UltraFresnelLensParameterisation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLensParameterisation::ComputeTransformation
// ComputeTransformation defines where each groove is placed
(const G4int GrooveNo, G4VPhysicalVolume* physVol) const
{
  G4double  Rmin1 = (GrooveNo+0)*(GrooveWidth) ;
  G4double  Rmax1 = (GrooveNo+1)*(GrooveWidth) ;
  //G4cout <<" groovewidth = "<< GrooveWidth/cm<<G4endl;
  G4double     dZ = FresnelLens->GetSagita(Rmax1) - FresnelLens->GetSagita(Rmin1) ;
  //dZ = 1.0*mm;    //xxx fixed for polycone

  if (dZ <= 0.0){
     G4Exception("UltraFresnelLensParameterisation::ComputeTransformation: Groove depth<0 !");
  }
 
  G4double LensSize = 5.525 * cm;
  //G4ThreeVector origin(LensSize/2.0 -mystery_offset, -LensSize/2.0 +mystery_offset, (dZ-dZOffset)/2.);  // why the extra 1.1 cm? xxx
  G4ThreeVector origin(LensSize/2.0, -LensSize/2.0, 0);
  //G4ThreeVector origin(0,0, (dZ-dZOffset)/2.);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLensParameterisation::ComputeDimensions
(G4Polycone& Groove, const G4int GrooveNo, const G4VPhysicalVolume*) const
{
  G4double  Rmin1 = (GrooveNo+0)*(GrooveWidth) ;
  G4double  Rmax1 = (GrooveNo+1)*(GrooveWidth) ;

  G4double  Rmin2 = Rmin1 ;
  G4double  Rmax2 = Rmin2+0.0001*mm ;
  
  G4double  dZ = FresnelLens->GetSagita(Rmax1) - FresnelLens->GetSagita(Rmin1) ;

  if (dZ <= 0.0){
     G4Exception("UltraFresnelLensParameterisation::ComputeDimensions: Groove depth<0 !");
  }

  G4double z_poly[3], rmin_poly[3], rmax_poly[3];
  z_poly[0] = -FresnelLens->LensThickness/2.0;
  z_poly[1] =  FresnelLens->LensThickness/2.0 - dZ;
  z_poly[2] =  FresnelLens->LensThickness/2.0;

  rmin_poly[0] = Rmin1;
  rmin_poly[1] = Rmin1;
  rmin_poly[2] = Rmin2;

  rmax_poly[0] = Rmax1;
  rmax_poly[1] = Rmax1;
  rmax_poly[2] = Rmax2;

 // Taper the outer part of the lens to a square shape, using PhiAngle (HvH, 2014)
    G4double circle_end = G4double(GrooveNo) / G4double(NumberOfGrooves);
    circle_end = circle_end*(NumberOfGrooves+1)/(NumberOfGrooves);    // small correction
    G4double phi0, phidelta;
    if (circle_end>99.6*sqrt(2.0)) {          // test marker at the tips 99.6->0.6
      phi0     = pi/2.0;
      phidelta = 0.1;
    }
    else if (circle_end<0.5*sqrt(2.0)) {      // quarter-circle, as before
      phi0     = pi/2.0;
      phidelta = pi/2.0;
    }
    else {                                    // taper to a point
      G4double start_angle = 0.5*asin( pow(circle_end,-2) -1.0 );
      phi0     = -start_angle+0.75*pi;
      phidelta =  2.0*start_angle;
    }
  G4Polycone foo("LensPV", phi0, phidelta, 3, z_poly, rmin_poly, rmax_poly);
  Groove = foo;



#ifdef ULTRA_VERBOSE

G4cout << "UltraFresnelLensParameterisation: GrooveNo " << GrooveNo+1 << 
" Rmin1, Rmax1(mm): " << Rmin1/mm <<" "<<  Rmax1/mm << " dZ(mm) " << dZ/mm << G4endl ;
#endif



}

