#ifndef PDG_HXX
#define PDG_HXX

#ifndef INVLD
#define INVLD -9999
#endif

#include "Constants.hxx"

/*** Global variables ***/

TString gNumberingScheme;

Int_t gGammaID;
Int_t gPositronID;
Int_t gElectronID;
Int_t gNeutrinoID;
Int_t gMuonPlusID;
Int_t gMuonMinusID;
Int_t gPiZeroID;
Int_t gPiPlusID;
Int_t gPiMinusID;
Int_t gKaonZeroLongID;
Int_t gKaonPlusID;
Int_t gKaonMinusID;
Int_t gNeutronID;
Int_t gProtonID;
Int_t gKaonZeroShortID;
Int_t gEtaID;
Int_t gOmegaID;
Int_t gKaonZeroID;

void UpdateNumberingScheme() {
  if (gNumberingScheme == "PDG") {
    gGammaID = 22;
    gPositronID = -11;
    gElectronID = 11;
    gNeutrinoID = 12;
    gMuonPlusID = -13;
    gMuonMinusID = 13;
    gPiZeroID = 111;
    gPiPlusID = 211;
    gPiMinusID = -211;
    gKaonZeroLongID = 130;
    gKaonPlusID = 321;
    gKaonMinusID = -321;
    gNeutronID = 2112;
    gProtonID = 2212;
    gKaonZeroShortID = 310;
    gEtaID = 221;
    gOmegaID = 223;
    gKaonZeroID = 311;
  } else if (gNumberingScheme == "GEANT") {
    gGammaID = 1;
    gPositronID = 2;
    gElectronID = 3;
    gNeutrinoID = 4;
    gMuonPlusID = 5;
    gMuonMinusID = 6;
    gPiZeroID = 7;
    gPiPlusID = 8;
    gPiMinusID = 9;
    gKaonZeroLongID = 10;
    gKaonPlusID = 11;
    gKaonMinusID = 12;
    gNeutronID = 13;
    gProtonID = 14;
    gKaonZeroShortID = 16;
    gEtaID = 17;
    gOmegaID = 223;     // doesn't exist in GEANT nomenclature
    gKaonZeroID = 311;  // doesn't exist in GEANT nomenclature
  }
}

void SetNumberingScheme(TString fOption) {
  gNumberingScheme = fOption;
  UpdateNumberingScheme();
}

Double_t ToPDG(Double_t fID) {
  // ISSUE: one's gotta be careful with this function
  // in gNumberingScheme == "GEANT", if you put a fID in PDG format, it wil collapse!
  if (gNumberingScheme == "PDG") {
    return fID;  // if detected numbering scheme was PDG, then there's no problem
  } else if (gNumberingScheme == "GEANT") {
    if (fID == gGammaID)
      return 22;
    else if (fID == gPositronID)
      return -11;
    else if (fID == gElectronID)
      return 11;
    else if (fID == gNeutrinoID)
      return 12;
    else if (fID == gMuonPlusID)
      return -13;
    else if (fID == gMuonMinusID)
      return 13;
    else if (fID == gPiZeroID)
      return 111;
    else if (fID == gPiPlusID)
      return 211;
    else if (fID == gPiMinusID)
      return -211;
    else if (fID == gKaonZeroLongID)
      return 130;
    else if (fID == gKaonPlusID)
      return 321;
    else if (fID == gKaonMinusID)
      return -321;
    else if (fID == gNeutronID)
      return 2112;
    else if (fID == gProtonID)
      return 2212;
    else if (fID == gKaonZeroShortID)
      return 310;
    else if (fID == gEtaID)
      return 221;
    else if (fID == gOmegaID)
      return 223;
    else if (fID == gKaonZeroID)
      return 311;
  }
}

Int_t particleID(TString particleName) {
  if (particleName == "gamma")
    return 22;
  else if (particleName == "positron")
    return -11;
  else if (particleName == "electron")
    return 11;
  else if (particleName == "neutrino")
    return 12;
  else if (particleName == "muon+")
    return -13;
  else if (particleName == "muon")
    return 13;
  else if (particleName == "pi0")
    return 111;
  else if (particleName == "pi+")
    return 211;
  else if (particleName == "pi-")
    return -211;
  else if (particleName == "K0L")
    return 130;
  else if (particleName == "K+")
    return 321;
  else if (particleName == "K-")
    return -321;
  else if (particleName == "neutron")
    return 2112;
  else if (particleName == "proton")
    return 2212;
  else if (particleName == "K0S")
    return 310;
  else if (particleName == "eta")
    return 221;
  else if (particleName == "omega")
    return 223;
  else if (particleName == "K0")
    return 311;
  else
    return 0;
}

Double_t particleMass(Double_t fPid) {
  // only receives PDG
  if (fPid == INVLD)
    return INVLD;
  else if (fPid == 22)
    return kMassGamma;
  else if (fPid == -11)
    return kMassPositron;
  else if (fPid == 11)
    return kMassElectron;
  else if (fPid == 12)
    return kMassNeutrino;
  else if (fPid == -13)
    return kMassMuonPlus;
  else if (fPid == 13)
    return kMassMuonMinus;
  else if (fPid == 111)
    return kMassPiZero;
  else if (fPid == 211)
    return kMassPiPlus;
  else if (fPid == -211)
    return kMassPiMinus;
  else if (fPid == 130)
    return kMassKaonZeroLong;
  else if (fPid == 321)
    return kMassKaonPlus;
  else if (fPid == -321)
    return kMassKaonMinus;
  else if (fPid == 2112)
    return kMassNeutron;
  else if (fPid == 2212)
    return kMassProton;
  else if (fPid == 310)
    return kMassKaonZeroShort;
  else if (fPid == 221)
    return kMassEta;
  else if (fPid == 223)
    return kMassOmega;
  else if (fPid == 311)
    return kMassKaonZero;
  else
    return INVLD;  // closure
}

#endif
