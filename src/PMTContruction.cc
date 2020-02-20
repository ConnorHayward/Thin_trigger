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

  G4double casingLength = 30*mm/2;
  G4double casingHeight = 32.5*mm/2;
  G4double casingHole = 29*mm / 2;
  G4double casingHoleHeight = 31.5*mm/2;

  G4double photoCathLength = 26*mm/2;
  G4double photoCathEffLength = 23*mm/2;
  G4double photoCathHeight = 0.8*mm;

  G4Box* photoCase = new G4Box("case",casingLength,casingLength,casingHeight);
  G4Box* pHole = new G4Box("hole",casingHole,casingHole,casingHoleHeight);
  G4SubtractionSolid* pmtShell = new G4SubtractionSolid("pmtShell",photoCase,pHole);

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateZ(45.*deg);
  G4RotationMatrix* rm1 = new G4RotationMatrix()

  G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Green());
  pmtCathLog->SetVisAttributes(detectorAttr);
  pmtInactiveCathLog->SetVisAttributes(detectorAttr);
  return final_plate;

}
