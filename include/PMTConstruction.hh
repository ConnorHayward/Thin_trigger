#ifndef PMTCONSTUCTION_H
#define PMTCONSTUCTION_H

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"

class PMTConstruction
{
  public:
    G4VSolid* ConstructPMT();

  private:
    G4double fSiliconPlate_h;
    G4double fHolderWidth;
};
#endif
