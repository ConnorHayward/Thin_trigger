#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "PMTConstruction.hh"

G4VSolid* PMTConstruction::ConstructPMT(){
  fSiliconPlate_h = 1.5*mm;
  fHolderWidth=90.00*mm;

  G4double TubsStartAngle =  0;
  G4double TubsSpanningAngle = 360 * deg;

  G4double casing_length = 30*mm/2;
  G4double casing_height = 32.5*mm/2;
  G4double casing_hole = 29*mm / 2;
  G4double casing_hole_height = 31.5*mm/2;

  G4double photo_cath_length = 26*mm/2;
  G4double photo_cath_eff_length = 23*mm/2;
  G4double photo_cath_height = 0.8*mm;

  G4Box* pmt_case = new G4Box("case",casing_length,casing_length,casing_height);
  G4Box* pHole = new G4Box("hole",casing_hole,casing_hole,casing_hole_height);
  G4SubtractionSolid* pmt_shell = new G4SubtractionSolid("pmt_shell",pmt_case,pHole);

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateZ(45.*deg);
  G4RotationMatrix* rm1 = new G4RotationMatrix()

  G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Green());
  pmt_cath_log->SetVisAttributes(detectorAttr);
  pmt_inactive_cath_log->SetVisAttributes(detectorAttr);
  return final_plate;

}
