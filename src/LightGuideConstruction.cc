#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "LightGuideConstruction.hh"

// Class to define logical / virtual volume for the light guide

G4VSolid* LightGuideConstruction::ConstructPlate(){
  fSiliconPlate_h = 1.5*mm;
  fHolderWidth=90.00*mm;

  G4double TubsStartAngle =  0;
  G4double TubsSpanningAngle = 360 * deg;

  G4double rect_x = 20*mm;
  G4double rect_y = 14*mm;
  G4double chamfer_size = 2.25*mm;
  G4double slit_height = 0.125*mm;
  G4double slit_offset = 1.75*mm;
  G4double slit_depth = 0.5*mm;

  G4double pmt_diameter = 25.5*mm;
  G4double pmt_depth = 1*mm;

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(45.*deg);
  G4RotationMatrix* rm1 = new G4RotationMatrix();
  rm1->rotateZ(-35.*deg);
  G4RotationMatrix* rm2 = new G4RotationMatrix();
  rm2->rotateZ(-35.*deg);
  rm2->rotateX(45.*deg);
  G4RotationMatrix* rm3 = new G4RotationMatrix();
  rm3->rotateY(90.*deg);

  char solidname[100];

  G4VSolid* initial_block = new G4Box("lg_block", rect_x, rect_y, rect_y);
  G4VSolid* angle_block = new G4Box("angle_block", 35*mm, 20*mm, 15*mm);

  G4VSolid* chamfer_edge = new G4Box("chamfer", rect_x+1*mm, chamfer_size, chamfer_size);
  G4VSolid* chamfer_two = new G4Box("chamfer_two", 35*mm, chamfer_size+1*mm, chamfer_size+1*mm);

  G4VSolid* slit = new G4Box("slit", slit_depth, slit_height, rect_y);

  G4VSolid* pmt_hole = new G4Tubs("dip", 0, pmt_diameter/2, pmt_depth/2, 0, TubsSpanningAngle);

  G4SubtractionSolid* one_edge = new G4SubtractionSolid("one_edge", initial_block, chamfer_edge, rm, G4ThreeVector(0*mm,rect_y, rect_y));
  one_edge = new G4SubtractionSolid("one_edge", one_edge, chamfer_edge, rm, G4ThreeVector(0*mm,-rect_y, rect_y));
  one_edge = new G4SubtractionSolid("one_edge", one_edge, chamfer_edge, rm, G4ThreeVector(0*mm,-rect_y, -rect_y));
  one_edge = new G4SubtractionSolid("one_edge", one_edge, chamfer_edge, rm, G4ThreeVector(0*mm, rect_y, -rect_y));

  G4SubtractionSolid* two_edge = new G4SubtractionSolid("two_edge", one_edge,angle_block,rm1, G4ThreeVector(rect_x,-rect_y, 0));

  G4SubtractionSolid* three_edge = new G4SubtractionSolid("three_edge", two_edge, chamfer_two, rm2, G4ThreeVector(+chamfer_size, -chamfer_size, chamfer_size+rect_y));
  three_edge = new G4SubtractionSolid("three_edge", three_edge, chamfer_two, rm2, G4ThreeVector(+chamfer_size, -chamfer_size, -chamfer_size-rect_y));
  three_edge = new G4SubtractionSolid("three_edge", three_edge, pmt_hole, rm3, G4ThreeVector(-rect_x+0.5*mm, 0, 0));

  G4SubtractionSolid* guide = new G4SubtractionSolid("final_guide", three_edge, slit, 0, G4ThreeVector(rect_x, rect_y-slit_offset, 0));

  G4VSolid* final_plate = guide;
  return final_plate;

}

// G4LogicalVolume* LightGuideConstruction::ConstructGuideLog(){
//   G4VSolid* light_guide = ConstructPlate();
//
//   G4NistManager* man = G4NistManager::Instance();
//
//   G4Material* PMMA = man->FindOrBuildMaterial("G4_PLEXIGLASS");
//
//   G4double refractive_index[] = {1.49, 1.49, 1.49, 1.49, 1.49};
//   G4double abs[] = {0.5*m, 0.5*m, 0.5*m, 0.5*m, 0.5*m};
//   G4double refl[] = {0.9, 0.9, 0.9, 0.9, 0.9};
//   G4double energy[] = {2.48*eV, 2.58*eV, 2.68*eV, 2.78*eV, 3.1*eV};
//
//   G4MaterialPropertiesTable* pmmaMPT = new G4MaterialPropertiesTable();
//   pmmaMPT->AddProperty("RINDEX", energy, refractive_index, 5);
//   pmmaMPT->AddProperty("ABSLENGTH", energy, abs, 5);
//   pmmaMPT->AddProperty("REFLECTIVITY", energy, refl, 5);
//   PMMA->SetMaterialPropertiesTable(pmmaMPT);
//
//   G4LogicalVolume* guide_log = new G4LogicalVolume(light_guide, PMMA, "Ligh_guide_log");
//   return guide_log;
// }
