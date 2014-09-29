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
// $Id: ExN06DetectorConstruction.cc,v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "UltraFresnelLens.hh"
#include "UltraFresnelLensParameterisation.hh"
#include <cmath>
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::ExN06DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::~ExN06DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN06DetectorConstruction::Construct() {


//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

  G4Element* H = new G4Element("Hydrogen", "H",  z=1 , a=1.01*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N",  z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O",  z=8 , a=16.00*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.00*g/mole);



  // Water
 //
   G4Material* Water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
   Water->AddElement(H, 2);
   Water->AddElement(O, 1);

//
// ------------ Generate & Add Material Properties Table ------------
//
const G4int nEntries = 32;

G4double PhotonEnergy[nEntries] =
         { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
           2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
           2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
           2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
           2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
           3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
           3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
           3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
//xxx
// Water
//xxx

  // Air
  //
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  // Aerogel ------------------------------------- aerogel ------------------------------------------
  //
  G4Material* Aerogel = new G4Material("Aerogel", density= 0.02*g/cm3, nelements=2);
  Aerogel->AddElement(Si, 1);
  Aerogel->AddElement(O, 2);

  const G4int aeEntries = 32;
  G4double AePhotonEnergy[aeEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,     // 610, 600, 590, 580, (nm)
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     // 570, 560, 550, 540,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,     // 530, 520, 510, 500,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,     // 490, 480, 470, 460,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,     // 450, 440, 430, 420,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,     // 410, 400, 390, 380,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,     // 370, 360, 350, 340,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };   // 330, 320, 310, 300.

  G4double AeRefractiveIndex1[aeEntries] =
            { 1.02435, 1.0244,  1.02445, 1.0245,  1.02455,
              1.0246,  1.02465, 1.0247,  1.02475, 1.0248,
              1.02485, 1.02492, 1.025,   1.02505, 1.0251,
              1.02518, 1.02522, 1.02530, 1.02535, 1.0254,
              1.02545, 1.0255,  1.02555, 1.0256,  1.02568,
              1.02572, 1.0258,  1.02585, 1.0259,  1.02595,
              1.026,   1.02608};

  G4double AeAbsorption1[aeEntries] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

//  G4double AeScintilFast[aeEntries] =
//            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
//              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
//              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
//              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
//              1.00, 1.00, 1.00, 1.00 };
//  G4double AeScintilSlow[aeEntries] =
//            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
//              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
//              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
//              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
//              7.00, 6.00, 5.00, 4.00 };

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX",       AePhotonEnergy, AeRefractiveIndex1,aeEntries);
  myMPT3->AddProperty("ABSLENGTH",    AePhotonEnergy, AeAbsorption1,     aeEntries);
//  myMPT3->AddProperty("FASTCOMPONENT",AePhotonEnergy, AeScintilFast,     aeEntries);
//  myMPT3->AddProperty("SLOWCOMPONENT",AePhotonEnergy, AeScintilSlow,     aeEntries);

  G4double AeRScatLength[aeEntries];
//  static const G4double AerogelTypeAClarity = 0.00719*micrometer*micrometer*micrometer*micrometer/cm;
  static const G4double AerogelTypeAClarity = 0.0020*micrometer*micrometer*micrometer*micrometer/cm;
  G4double Cparam    =  AerogelTypeAClarity*cm/(micrometer*micrometer*micrometer*micrometer);
  G4double PhotMomWaveConv = 1239*eV*nm;
  if(Cparam != 0.0 ) {
    for(G4int ibinw=0; ibinw<aeEntries; ibinw++ ){
      G4double ephoton = AePhotonEnergy[ibinw];
      //In the following the 1000 is to convert form nm to micrometer
      G4double wphoton=(PhotMomWaveConv/ephoton)/(1000.0*nm);
      AeRScatLength[ibinw]=(std::pow(wphoton,4))/Cparam;
    }
  }
  myMPT3->AddProperty("RAYLEIGH",     AePhotonEnergy, AeRScatLength,     aeEntries);

  myMPT3->AddConstProperty("SCINTILLATIONYIELD",0./MeV);
  myMPT3->AddConstProperty("RESOLUTIONSCALE",1.0);
//  myMPT3->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
//  myMPT3->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
//  myMPT3->AddConstProperty("YIELDRATIO",0.8);

  Aerogel->SetMaterialPropertiesTable(myMPT3);


// --- end aerogel tables --------------------------------------- end aerogel -------------------------

  // PMMA C5H8O2 ( Acrylic )
   // -------------
      density = 1.19*g/cm3;
      G4Material* Acrylic = new G4Material("Acrylic", density, nelements=3);

      a = 12.01*g/mole;
      G4Element* C  = new G4Element("Carbon"  ,"C" , z= 6., a);
      Acrylic->AddElement(C, 5);
      Acrylic->AddElement(H, 8);     // molecular ratios
      Acrylic->AddElement(O, 2);

      G4double AcRefractiveIndex[aeEntries] =
                { 1.4902,  1.4907,  1.4913,  1.4918,  1.4924,     // 610, 600, 590, 580, 570,
                  1.4930,  1.4936,  1.4942,  1.4948,  1.4954,     // 560, 550, 540, 530, 520,  (this line is interpolated)
                  1.4960,  1.4965,  1.4971,  1.4977,  1.4983,     // 510, 500, 490, 480, 470,
                  1.4991,  1.5002,  1.5017,  1.5017,  1.5017,     // 460, 450, 440, 430, 420,
                  1.5017,  1.5017,  1.5017,  1.5017,  1.5017,     // 410,
                  1.5017,  1.5017,  1.5017,  1.5017,  1.5017,     // 360,     look up values below 435
                  1.5017,  1.5017, };                             // 310, 300.

           G4double AcAbsorption[aeEntries] =
                   { 00.448*cm,  00.082*cm,  00.329*cm,  00.174*cm,  00.346*cm,  00.889*cm,
                     00.152*cm,  00.241*cm,  00.868*cm,  00.000*cm,  00.316*cm,  00.714*cm,
                     00.455*cm,  00.619*cm,  00.632*cm,  00.632*cm,  00.556*cm,  00.632*cm,
                     00.632*cm,  00.619*cm,  00.455*cm,  00.667*cm,  00.037*cm,  00.333*cm,
                     00.001*cm,  00.001*cm,  00.001*cm,  00.001*cm,  00.001*cm,  00.001*cm,
                     00.001*cm,  00.001*cm };                                              // cutoff below 300 nm

      G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();
      myMPT4->AddProperty("RINDEX",       AePhotonEnergy, AcRefractiveIndex,aeEntries);
      myMPT4->AddProperty("ABSLENGTH",    AePhotonEnergy, AcAbsorption,     aeEntries);

      Acrylic->SetMaterialPropertiesTable(myMPT4);

      //
// Air ------------------------------------------------------------------------- Air:
//
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  
  Air->SetMaterialPropertiesTable(myMPT2);

  //
  // Aluminum --------------------------------------------------------------------------Aluminum:
  //
  density = 2.7*g/cm3;    // Aluminum
  G4Material* Aluminum = new G4Material("Aluminum", density, nelements=1);
  G4Element* Al = new G4Element("Aluminum"  ,"Al" , z= 13.0, a=26.98*g/mole);
  Aluminum->AddElement(Al, 1);

  //
///////////////////////////////////////////////   VOLUMES   /////////////////////////////////////////////////////////////////////

// The experimental Hall -------------------------------------------------------- hall ---------------
//
  expHall_x = expHall_y = expHall_z = 10.0*cm;
  G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  expHall_log -> SetVisAttributes (G4VisAttributes::GetInvisible());

////////////////////////////////////////////  Define key switches, sizes and locations   /////////////////////////////////////

  G4int make_agel    = 1;
  G4int make_lens    = 1;
  G4int make_mirrors = 1;

  G4double agel_halfx = 5.525*cm;
  G4double agel_halfy = agel_halfx;
  G4double agel_halfz = 1.0*cm;

  G4double BoxDelz = -2.0*mm;

  G4double lens_z = -2.5*cm + BoxDelz;
  G4double phodet_z     = 4.6*cm + BoxDelz;
  G4double phodet_halfz = 0.1*cm;
  //G4double readout_z = 10;
  G4double readout_halfz = 4.0;
  G4double readout_z     [2] = {phodet_z-phodet_halfz+3.0, phodet_z-phodet_halfz + 2.0*readout_halfz};
  G4double   LensDiameter        = 2.0 * agel_halfx * sqrt(2.0) ; // Size of the optical active area of the lens.


  ///////////////////////////////////////  HOLDER BOX   /////////////////////////////////////////////////////////

  G4double box_halfx = agel_halfx + 0.1*cm ;
  G4double box_halfy = box_halfx;
  G4double box_halfz = ( -lens_z + 2.0*agel_halfz + readout_z[1] + readout_halfz + 0.0*cm )/2.0;

  G4double box_x = 0.0*cm;
  G4double box_y = 0.0*cm;
  G4double box_z = 0.0*mm;

  G4Box* box_box = new G4Box("box", box_halfx, box_halfy, box_halfz);

  G4LogicalVolume* box_log = new G4LogicalVolume(box_box,Air,"box",0,0,0);

  G4VPhysicalVolume* box_phys
      = new G4PVPlacement(0,G4ThreeVector(box_x, box_y, box_z), box_log,"box",
                          expHall_log,false,0);

  //box_log -> SetVisAttributes (G4VisAttributes::GetInvisible());
  G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.0));
  SurfaceVisAtt->SetVisibility(true);
  SurfaceVisAtt->SetForceWireframe(true);
  box_log->SetVisAttributes(SurfaceVisAtt);


  ///////////////////////////////////////  TEST BOOLEAN OPERATION /////////////////////////////////////////////////////////


  ///////////////////////////////////////  AEROGEL BLOCK   /////////////////////////////////////////////////////////

  G4double agel_px;  G4double agel_py;  G4double agel_pz;
  if (make_agel == 1) {
	agel_halfx = agel_halfy;
	agel_px = 0.0*cm;
	agel_py = 0.0*cm;
	agel_pz = lens_z - agel_halfz - 1.5*mm + BoxDelz;
  }  // end make_agel
  else {
	    agel_halfx = 0.1*cm;   agel_halfy = 0.1*cm;   agel_halfz = 0.1*cm;
	    agel_px = 5.0*cm;	    agel_py = 5.0*cm;	    agel_pz =-8.0*cm;
  }

  G4Box* aerogel_box = new G4Box("agel",agel_halfx,agel_halfy,agel_halfz);

    G4LogicalVolume* aerogel_log
      = new G4LogicalVolume(aerogel_box,Aerogel,"agel",0,0,0);

    G4VPhysicalVolume* aerogel_phys
        = new G4PVPlacement(0,G4ThreeVector(agel_px,agel_py,agel_pz),aerogel_log,"agel",
                            box_log,false,0);

    SurfaceVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    SurfaceVisAtt->SetVisibility(true);
    SurfaceVisAtt->SetForceWireframe(true);
    aerogel_log->SetVisAttributes(SurfaceVisAtt);

////////////////////////////////////////////////////////// FRESNEL LENS ////////////////////////////////////////////////////////
//
// from /sw/share/geant4.9/examples/advanced/air_shower/src/UltraDetectorConstruction.cc

      G4int      LensNumOfGrooves    = 100;    // 400
      G4String name;

      G4Material   *LensMaterial        = G4Material::GetMaterial(name = "Acrylic") ;

      if (make_lens == 1) {
        G4ThreeVector LensPosition        = G4ThreeVector(0.0*mm,0.0*mm,lens_z+BoxDelz) ;
        new UltraFresnelLens(LensDiameter,LensNumOfGrooves,LensMaterial,box_phys,LensPosition) ;
      }   // end make_lens

////////////////////////////////////////////////////////// THE PHOTON DETECTOR    //////////////////////////////////////////////

        G4double factor = 0.8;
        G4double phodet_halfx  =  agel_halfx * factor;
        G4double phodet_halfy  =  agel_halfy * factor;
        G4Box* phoDet_box = new G4Box("phodet",phodet_halfx,phodet_halfy,phodet_halfz);

        G4LogicalVolume* phoDet_log
          = new G4LogicalVolume(phoDet_box,Aluminum,"phodet",0,0,0);

        G4VPhysicalVolume* phoDet_phys =
                new G4PVPlacement(0,G4ThreeVector(0,0,phodet_z),phoDet_log,"phodet",       // for lens_type=4
                                  box_log,false,0);

////////////////////////////////////////////////////////////   MIRRORS   ////////////////////////////////////////////////////////

        if (make_mirrors ==1 ) {

          G4double dx1 = agel_halfx;
          G4double dx2 = agel_halfx * factor;
          G4double dy1 = 0.01*cm;
          G4double dy2 = 0.01*cm;
          G4double dz = (phodet_z - lens_z - phodet_halfz - 0.3*cm)/2.0;
          G4double phi = atan2(agel_halfx-phodet_halfx,2.0*dz);
          G4double delxy = dz * sin(phi) + 0.1*cm;
          dz = sqrt(dz*dz + pow(agel_halfx-phodet_halfx,2) );

          G4Trd * trdMirror = new G4Trd ("trap",dx1, dx2, dy1, dy2, dz);
          G4LogicalVolume *logTrd ;
          logTrd = new G4LogicalVolume(trdMirror,Aluminum,"trap",0,0,0);
          G4VPhysicalVolume *physMirror5 ;

          G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
          SurfaceVisAtt->SetVisibility(true);
          SurfaceVisAtt->SetForceWireframe(true);
          logTrd->SetVisAttributes(SurfaceVisAtt);

          G4double mirror_halfx = agel_halfx;
          G4double mirror_halfy = 0.1*cm;
          G4double mirror_halfz = (phodet_z - lens_z - phodet_halfz - 0.1*cm)/2.0;
          //G4double lenshalfz = FresnelLens.LensThickness/2.0;
          G4double lens_halfz = 0.3*cm;

// 5 back mirror
          G4double Mirror_py = 0.0*cm;
          G4double Mirror_px = agel_halfy + mirror_halfy - delxy;
          G4double Mirror_pz = (lens_z + lens_halfz+(phodet_z-phodet_halfz))/2.0 + BoxDelz;
          G4ThreeVector SurfacePosition = G4ThreeVector(Mirror_px,Mirror_py,Mirror_pz) ;
          G4RotationMatrix *Surfrot5 = new G4RotationMatrix();
          Surfrot5->rotateZ(pi*0.5);
          Surfrot5->rotateX(-phi);
          physMirror5 = new G4PVPlacement(Surfrot5,SurfacePosition,"MirrorPV",logTrd,box_phys,false,0,false);
// 6 front
          Mirror_px = -Mirror_px;
          SurfacePosition = G4ThreeVector(Mirror_px,Mirror_py,Mirror_pz) ;
          G4VPhysicalVolume *physMirror6 ;
          G4RotationMatrix *Surfrot6 = new G4RotationMatrix();
          Surfrot6->rotateZ(pi*0.5);
          Surfrot6->rotateX(+phi);
          physMirror6 = new G4PVPlacement(Surfrot6,SurfacePosition,"MirrorPV",logTrd,box_phys,false,0,false);
// 7 top
          Mirror_py = -Mirror_px;
          Mirror_px = 0;
          SurfacePosition = G4ThreeVector(Mirror_px,Mirror_py,Mirror_pz) ;
          G4VPhysicalVolume *physMirror7 ;
          G4RotationMatrix *Surfrot7 = new G4RotationMatrix();
          Surfrot7->rotateX(-phi);
          physMirror7 = new G4PVPlacement(Surfrot7,SurfacePosition,"MirrorPV",logTrd,box_phys,false,0,false);
// 8 bottom
          Mirror_py = -Mirror_py;
          SurfacePosition = G4ThreeVector(Mirror_px,Mirror_py,Mirror_pz) ;
          G4VPhysicalVolume *physMirror8 ;
          G4RotationMatrix *Surfrot8 = new G4RotationMatrix();
          Surfrot8->rotateX(+phi);
          physMirror8 = new G4PVPlacement(Surfrot8,SurfacePosition,"MirrorPV",logTrd,box_phys,false,0,false);

          //////////////////////////////////////////////////////////////////////////////////////////
          //   Optical properties of the interface between the Air and Reflective Surface
          //   For Mirror, reflectivity is set at 95% and specular reflection is assumed.

          G4OpticalSurface *OpticalAirMirror = new G4OpticalSurface("AirMirrorSurface");
          OpticalAirMirror->SetModel(unified);
          OpticalAirMirror->SetType(dielectric_dielectric);
          OpticalAirMirror->SetFinish(polishedfrontpainted);

          const G4int NUM = 2;
          G4double lambda_min = 200*nm ;
          G4double lambda_max = 700*nm ;

          G4double XX[NUM] = {h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;
          G4double ICEREFLECTIVITY[NUM]      = { 0.95, 0.95 };

          G4MaterialPropertiesTable *AirMirrorMPT = new G4MaterialPropertiesTable();
          AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
          OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);

          new G4LogicalBorderSurface("Air/Mirror Surface",box_phys,physMirror5,OpticalAirMirror);
          new G4LogicalBorderSurface("Air/Mirror Surface",box_phys,physMirror6,OpticalAirMirror);
          new G4LogicalBorderSurface("Air/Mirror Surface",box_phys,physMirror7,OpticalAirMirror);
          new G4LogicalBorderSurface("Air/Mirror Surface",box_phys,physMirror8,OpticalAirMirror);

        }     // end if make mirrors

////////////////////////////////////////////  READOUT HARDWARE //////////////////////////////////////////////////

        G4double readout_rinner[2] = {phodet_halfx+1, phodet_halfx+1};
        G4double readout_router[2] = {agel_halfx, agel_halfx};
        G4Polyhedra* readout = new G4Polyhedra("readout", 0.25*pi, 2.0*pi, 4, 2, readout_z, readout_rinner, readout_router);

        G4LogicalVolume* readout_log = new G4LogicalVolume(readout,Aluminum,"readout",0,0,0);

        G4VPhysicalVolume* readout_phys =
        new G4PVPlacement(0,G4ThreeVector(0,0,0),readout_log,"readout",
                          box_log,false,0);

        G4VisAttributes* readoutVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));    // red
        readoutVisAtt->SetVisibility(true);
        readoutVisAtt->SetForceWireframe(true);
        readout_log->SetVisAttributes(readoutVisAtt);


/////////////////////////////////////////////////   END VOLUMES   ////////////////////////////////////////////////////////////////////////
//
////  ////  ////  ////  ////  ////  ////  ////  //// SURFACES ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
//
// Agel
//
  G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
  OpWaterSurface->SetType(dielectric_dielectric);
  OpWaterSurface->SetFinish(ground);
  OpWaterSurface->SetModel(unified);

  new G4LogicalBorderSurface("WaterSurface",
                                 aerogel_phys,expHall_phys,OpWaterSurface);


// Air plane
//
  G4OpticalSurface* OpAirSurface = new G4OpticalSurface("AirSurface");
  OpAirSurface->SetType(dielectric_dielectric);
  OpAirSurface->SetFinish(polished);
  OpAirSurface->SetModel(glisur);

//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double Ephoton[num] = {2.034*eV, 4.136*eV};

  //OpticalWaterSurface 
  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);

  //OpticalAirSurface
  G4double Reflectivity[num] = {0.3, 0.5};
  G4double Efficiency[num]   = {0.8, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  myST2->AddProperty("EFFICIENCY",   Ephoton, Efficiency,   num);

  OpAirSurface->SetMaterialPropertiesTable(myST2);

//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
