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
//    *      UltraFresnelLens.cc
//    ****************************************************
//
//    Class for definition of the Ultra Fresnel Lens. 
//    An  UltraFresnelLens object is created in the UltraDetectorConstruction class.
//    This class makes use of the UltraFresnelLensParameterisation class.
//    The lens profile is define through the GetSagita method.   
//
#include <cmath>
#include "UltraFresnelLens.hh"
#include "UltraFresnelLensParameterisation.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4IntersectionSolid.hh"

UltraFresnelLens::UltraFresnelLens( G4double Diameter, G4int nGrooves, G4Material* Material, G4VPhysicalVolume * MotherPV, G4ThreeVector Pos)
{

 LensMaterial    = Material ;
 G4cout<<" ***** Material *** "<<LensMaterial<<" *****"<<G4endl;
 LensDiameter    = Diameter ;
 NumberOfGrooves = nGrooves;
 GrooveWidth     = (Diameter/2.0)/nGrooves ;
 LensPosition    = Pos ;

 if( GrooveWidth <= 0 ){
       G4Exception("UltraFresnelLens constructor: GrooveWidth<=0");
   }

 G4double  Rmin1 = (NumberOfGrooves-1)*(GrooveWidth) ;
 G4double  Rmax1 = (NumberOfGrooves-0)*(GrooveWidth) ;
 LensThickness   = GetSagita(Rmax1)-GetSagita(Rmin1) ; // Height of the highest groove
 LensThickness +=1.0*mm;          // minimum acrylic thickness at the deepest part of the outermost groove

BuildLens(MotherPV) ;

}


UltraFresnelLens::~UltraFresnelLens(  )
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLens::BuildLens(G4VPhysicalVolume *MotherPV){

G4double StartPhi = pi/2.0 ;     // modified to describe a quarter
G4double DeltaPhi = pi/2.0 ;     // sector of the lens

G4double  LensMotherDz = LensThickness/2.0 + 0.25*mm ;

G4Material *LensMotherMaterial = MotherPV->GetLogicalVolume()->GetMaterial() ;
G4cout<<"\n*** in BuildLens: mother vol made of ***"<< LensMotherMaterial<<" *************\n";
G4cout<<"\n*** in BuildLens: lens made of ***"<< LensMaterial<<" *************\n";

G4VPVParameterisation *FresnelLensParam = new UltraFresnelLensParameterisation(this);

//////////////////////////////////////////////  1  //////////////////////////////////////////////////////////////////
// The lens is built in four 90-degree sections.
// first quarter, holder box:
G4double size2 = sqrt(2.0)*LensDiameter/8.0;
G4Box* LensMotherBox1 = new G4Box("LensMotherBox1",size2, size2, LensMotherDz);
G4LogicalVolume *LensMotherBox1_log = new G4LogicalVolume(LensMotherBox1,LensMotherMaterial,"LensMotherBox1_log",0,0,0);
G4RotationMatrix *rot1 = new G4RotationMatrix();
rot1->rotateZ(pi/2.0);

G4VPhysicalVolume *LensMotherBox1_phys
    = new G4PVPlacement(rot1,G4ThreeVector(size2, size2, LensPosition[2]),"LensMotherBox1_phys",LensMotherBox1_log,MotherPV,false,0);
G4VisAttributes* SurfaceVisAtt1 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
SurfaceVisAtt1->SetVisibility(true);
SurfaceVisAtt1->SetForceWireframe(true);
LensMotherBox1_log->SetVisAttributes(SurfaceVisAtt1);

// first quarter, the lens POLYCONE:
G4double poly_z[2] = {0, 1.*mm};
G4double poly_rmin[2] = {40.0*mm, 40.0*mm};
G4double poly_rmax[2] = {50.0*mm, 40.001*mm};
G4Polycone *solidGroove1 = new G4Polycone("Groove1",StartPhi+pi,DeltaPhi, 2, poly_z, poly_rmin, poly_rmax);
G4LogicalVolume *logicalGroove1 = new G4LogicalVolume(solidGroove1,LensMaterial,"Groove_log1",0,0,0);
LensPhysicalVolume = new G4PVParameterised("LensPV2",logicalGroove1,LensMotherBox1_phys,kZAxis,NumberOfGrooves,FresnelLensParam) ;

//////////////////////////////////////////////  2  //////////////////////////////////////////////////////////////////
// second quarter, holder box:
G4cout<< "ultrafresnel: " << LensDiameter<< G4endl;

G4Box* LensMotherBox2 = new G4Box("LensMotherBox2",size2, size2, LensMotherDz);
G4LogicalVolume *LensMotherBox2_log = new G4LogicalVolume(LensMotherBox2,LensMotherMaterial,"LensMotherBox2_log",0,0,0);
G4VPhysicalVolume *LensMotherBox2_phys
    = new G4PVPlacement(0,G4ThreeVector(-size2, size2, LensPosition[2]),"LensMotherBox2_phys",LensMotherBox2_log,MotherPV,false,0);
//      = new G4PVPlacement(0,G4ThreeVector(0,0,0),"LensMotherBox2_phys",LensMotherBox2_log,MotherPV,false,0);

G4VisAttributes* SurfaceVisAtt2 = new G4VisAttributes(G4Colour(0.0,1.0,.0));  // 010 green
SurfaceVisAtt2->SetVisibility(true);
SurfaceVisAtt2->SetForceWireframe(true);
LensMotherBox2_log->SetVisAttributes(SurfaceVisAtt2);

// second quarter, the lens POLYCONE:
G4Polycone *solidGroove2 = new G4Polycone("Groove2",StartPhi+pi,DeltaPhi, 2, poly_z, poly_rmin, poly_rmax);
G4LogicalVolume *logicalGroove2 = new G4LogicalVolume(solidGroove2,LensMaterial,"Groove_log2",0,0,0);
LensPhysicalVolume = new G4PVParameterised("LensPV2",logicalGroove2,LensMotherBox2_phys,kZAxis,NumberOfGrooves,FresnelLensParam) ;

//////////////////////////////////////////////  3  //////////////////////////////////////////////////////////////////
// third quarter, holder box:
G4Box* LensMotherBox3 = new G4Box("LensMotherBox3",size2, size2, LensMotherDz);
G4LogicalVolume *LensMotherBox3_log = new G4LogicalVolume(LensMotherBox3,LensMotherMaterial,"LensMotherBox3_log",0,0,0);
G4RotationMatrix *rot3 = new G4RotationMatrix();
rot3->rotateZ(-pi*0.5);

G4VPhysicalVolume *LensMotherBox3_phys
    = new G4PVPlacement(rot3,G4ThreeVector(-size2, -size2, LensPosition[2]),"LensMotherBox3_phys",LensMotherBox3_log,MotherPV,false,0);
G4VisAttributes* SurfaceVisAtt3 = new G4VisAttributes(G4Colour(1.0,0.0,.0));
SurfaceVisAtt3->SetVisibility(true);
SurfaceVisAtt3->SetForceWireframe(true);
LensMotherBox3_log->SetVisAttributes(SurfaceVisAtt3);

// third quarter, the lens POLYCONE:
G4Polycone *solidGroove3 = new G4Polycone("Groove3",StartPhi+pi,DeltaPhi, 2, poly_z, poly_rmin, poly_rmax);
G4LogicalVolume *logicalGroove3 = new G4LogicalVolume(solidGroove3,LensMaterial,"Groove_log3",0,0,0);
LensPhysicalVolume = new G4PVParameterised("LensPV2",logicalGroove3,LensMotherBox3_phys,kZAxis,NumberOfGrooves,FresnelLensParam) ;

//////////////////////////////////////////////  4  //////////////////////////////////////////////////////////////////
// fourth quarter, holder box:
//G4double size2 = sqrt(2.0)*LensDiameter/8.0;
G4Box* LensMotherBox4 = new G4Box("LensMotherBox4",size2, size2, LensMotherDz);
G4LogicalVolume *LensMotherBox4_log = new G4LogicalVolume(LensMotherBox4,LensMotherMaterial,"LensMotherBox4_log",0,0,0);
G4RotationMatrix *rot4 = new G4RotationMatrix();
rot4->rotateZ(-pi);

G4VPhysicalVolume *LensMotherBox4_phys
    = new G4PVPlacement(rot4,G4ThreeVector(size2, -size2, LensPosition[2]),"LensMotherBox4_phys",LensMotherBox4_log,MotherPV,false,0);
G4VisAttributes* SurfaceVisAtt4 = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
SurfaceVisAtt4->SetVisibility(true);
SurfaceVisAtt4->SetForceWireframe(true);
LensMotherBox4_log->SetVisAttributes(SurfaceVisAtt4);

// fourth quarter, the lens POLYCONE:
G4Polycone *solidGroove4 = new G4Polycone("Groove1",StartPhi+pi,DeltaPhi, 2, poly_z, poly_rmin, poly_rmax);
G4LogicalVolume *logicalGroove4 = new G4LogicalVolume(solidGroove4,LensMaterial,"Groove_log0",0,0,0);
LensPhysicalVolume = new G4PVParameterised("LensPV2",logicalGroove4,LensMotherBox4_phys,kZAxis,NumberOfGrooves,FresnelLensParam) ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//==============================================================================

G4double UltraFresnelLens::GetSagita(G4double radius) 
{

  G4double Conic = -1.0;                                  // original:
  G4int lens_type = 3;
  G4double Curvature;
  G4double Aspher[8] = {0,0,0,0,0,0,0,0};

  if (lens_type == 1) {
	Curvature = 0.00437636761488/mm;
	Aspher[0] = 4.206739256e-05/(mm);
	Aspher[1] = 9.6440152e-10/(mm3);
	Aspher[2] = -1.4884317e-15/(mm2*mm3);
  }

if (lens_type == 2) {
  Curvature = 0.0132/mm;                // r=77mm, f~14cm
  Aspher[0] = 32.0e-05/(mm);
  Aspher[1] = -2.0e-7/(mm3);
  Aspher[2] =  1.2e-13/(mm2*mm3);
}

if (lens_type == 3) {
  Curvature = 0.0150/mm;         // r=77mm, f~12.5cm
  Aspher[0] = 42.0e-05/(mm);
  Aspher[1] = -3.0e-7/(mm3);
  Aspher[2] =  1.2e-13/(mm2*mm3);
}
if (lens_type == 4) {
  Curvature = 0.0175/mm;         // r=77mm, f~10cm
  Aspher[0] = 72.0e-05/(mm);
  Aspher[1] = -5.0e-7/(mm3);
  Aspher[2] =  1.2e-13/(mm2*mm3);
}

  G4double TotAspher = 0.0*mm ;

  for(G4int k=1;k<9;k++){
    TotAspher += Aspher[k-1]*std::pow(radius,2*k) ;
  }

  G4double ArgSqrt = 1.0-(1.0+Conic)*std::pow(Curvature,2)*std::pow(radius,2) ; // note conic=-1, so ArgSqrt = 1.0

  if (ArgSqrt < 0.0){
     G4Exception("UltraFresnelLensParameterisation::Sagita: Square Root of <0 !");
  }
  G4double Sagita_value = Curvature*std::pow(radius,2)/(1.0+std::sqrt(ArgSqrt)) + TotAspher;

  return Sagita_value ;

                             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


