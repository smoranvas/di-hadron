#include "GSTtree.h"
#include "TIdentificatorV2.h"
#include "Constants.hxx"
#include "Headers.hxx"
#include "Math.hxx"
#include "PDG.hxx"

#ifndef UI_HXX
#include "UI.hxx"
#endif

#ifndef TREEOPS_HXX
#define TREEOPS_HXX

void SetOutputBranches_REC(TTree* tree, rec_p& rec) {
  // reconstructed electron (46 variables)
  
  tree->Branch("Q2", &rec.Q2);
  tree->Branch("Nu", &rec.Nu);
  tree->Branch("Xb", &rec.Xb);
  tree->Branch("vxe", &rec.vxe);
  tree->Branch("vye", &rec.vye);
  tree->Branch("vze", &rec.vze);
  tree->Branch("SectorEl", &rec.SectorEl);
  tree->Branch("TargType", &rec.TargType);
  tree->Branch("Pex", &rec.Pex);
  tree->Branch("Pey", &rec.Pey);
  tree->Branch("Pez", &rec.Pez);
  tree->Branch("Pe", &rec.Pe);
  tree->Branch("vyec", &rec.vyec);
  tree->Branch("vzec", &rec.vzec);  
  tree->Branch("NpheEl", &rec.NpheEl);
  tree->Branch("FidCheckCutEl", &rec.FidCheckCutEl); // added by smoran
  tree->Branch("TargTypeSM", &rec.TargTypeSM);
  tree->Branch("TargTypeHH", &rec.TargTypeHH);  // added by smoran
  tree->Branch("CCMirrorMatchingEl", &rec.CCMirrorMatchingEl);  // added by smoran
  tree->Branch("SampFractionEl25", &rec.SampFractionEl25);  // added by smoran
  tree->Branch("SampFractionEl20", &rec.SampFractionEl20);  // added by smoran
  tree->Branch("W", &rec.W);

 
  //tree->Branch("W", &rec.W);
  //tree->Branch("Yb", &rec.Yb);
  //tree->Branch("BettaEl", &rec.BettaEl);
  //tree->Branch("Etote", &rec.Etote);
  //tree->Branch("Eine", &rec.Eine);
  //tree->Branch("Eoute", &rec.Eoute);
  //tree->Branch("vxec", &rec.vxec);  
  //tree->Branch("XECe", &rec.XECe);
  //tree->Branch("YECe", &rec.YECe);
  //tree->Branch("ZECe", &rec.ZECe);
  //tree->Branch("PhiLabEl", &rec.PhiLabEl);
  //tree->Branch("ThetaLabEl", &rec.ThetaLabEl);
  //tree->Branch("StatDCEl", &rec.StatDCEl);
  //tree->Branch("DCStatusEl", &rec.DCStatusEl);
  //tree->Branch("StatECEl", &rec.StatECEl);
  //tree->Branch("ECStatusEl", &rec.ECStatusEl);
  //tree->Branch("TimeECEl", &rec.TimeECEl);
  //tree->Branch("PathECEl", &rec.PathECEl);
  //tree->Branch("Chi2ECEl", &rec.Chi2ECEl);
  //tree->Branch("StatSCEl", &rec.StatSCEl);
  //tree->Branch("SCStatusEl", &rec.SCStatusEl);
  //tree->Branch("TimeSCEl", &rec.TimeSCEl);
  //tree->Branch("PathSCEl", &rec.PathSCEl);
  //tree->Branch("StatCCEl", &rec.StatCCEl);
  //tree->Branch("CCStatusEl", &rec.CCStatusEl);
  //tree->Branch("Chi2CCEl", &rec.Chi2CCEl);
  //tree->Branch("StatusEl", &rec.StatusEl);
  //tree->Branch("NRowsDCEl", &rec.NRowsDCEl);
  //tree->Branch("NRowsECEl", &rec.NRowsECEl);
  //tree->Branch("NRowsSCEl", &rec.NRowsSCEl);
  //tree->Branch("NRowsCCEl", &rec.NRowsCCEl);  

  // reconstructed particle (49 variables)
  
  tree->Branch("Zh", &rec.Zh);
  tree->Branch("Pt2", &rec.Pt2);
  tree->Branch("PhiPQ", &rec.PhiPQ);
  tree->Branch("ThetaLab", &rec.ThetaLab);
  tree->Branch("vxh", &rec.vxh);
  tree->Branch("vyh", &rec.vyh);
  tree->Branch("vzh", &rec.vzh);
  tree->Branch("Sector", &rec.Sector);
  tree->Branch("Px", &rec.Px);
  tree->Branch("Py", &rec.Py);
  tree->Branch("Pz", &rec.Pz);
  tree->Branch("P", &rec.P);
  tree->Branch("Betta", &rec.Betta);
  tree->Branch("pid", &rec.pid);
  tree->Branch("T4", &rec.T4);
  tree->Branch("T4_pi", &rec.T4_pi);
  tree->Branch("deltaZ", &rec.deltaZ);
  tree->Branch("Nphe", &rec.Nphe);
  tree->Branch("evnt", &rec.evnt); 
  tree->Branch("FidCheckCutPiPlus", &rec.FidCheckCutPiPlus); // added by smoran
  tree->Branch("StatCC", &rec.StatCC);
  tree->Branch("Chi2CC", &rec.Chi2CC);
  tree->Branch("NRowsCC", &rec.NRowsCC);   
  tree->Branch("FidCheckCut", &rec.FidCheckCut);            // added by smoran
  tree->Branch("SampFraction25", &rec.SampFraction25);      // added by smoran
  tree->Branch("SampFraction20", &rec.SampFraction20);      // added by smoran
  tree->Branch("CCMirrorMatching", &rec.CCMirrorMatching);  // added by smoran
  tree->Branch("Eh", &rec.Eh);
  tree->Branch("ThetaPQ", &rec.ThetaPQ);  
  tree->Branch("PhiLab", &rec.PhiLab);

  //event-related (1 variable)
  //tree->Branch("Pl2", &rec.Pl2);
  //tree->Branch("Mx2", &rec.Mx2);
  //tree->Branch("T", &rec.T);
  //tree->Branch("Mass2", &rec.Mass2);
  //tree->Branch("Etot", &rec.Etot);
  //tree->Branch("Ein", &rec.Ein);
  //tree->Branch("Eout", &rec.Eout);
  //tree->Branch("XEC", &rec.XEC);
  //tree->Branch("YEC", &rec.YEC);
  //tree->Branch("ZEC", &rec.ZEC);   
  //tree->Branch("StatDC", &rec.StatDC);
  //tree->Branch("DCStatus", &rec.DCStatus);
  //tree->Branch("StatEC", &rec.StatEC);
  //tree->Branch("ECStatus", &rec.ECStatus);
  //tree->Branch("TimeEC", &rec.TimeEC);
  //tree->Branch("PathEC", &rec.PathEC);
  //tree->Branch("Chi2EC", &rec.Chi2EC);
  //tree->Branch("StatSC", &rec.StatSC);
  //tree->Branch("SCStatus", &rec.SCStatus);
  //tree->Branch("TimeSC", &rec.TimeSC);
  //tree->Branch("PathSC", &rec.PathSC);
  //tree->Branch("CCStatus", &rec.CCStatus);
  //tree->Branch("Status", &rec.Status);
  //tree->Branch("NRowsDC", &rec.NRowsDC);
  //tree->Branch("NRowsEC", &rec.NRowsEC);
  //tree->Branch("NRowsSC", &rec.NRowsSC);
}

void SetOutputBranches_GEN(TTree* tree, gen_p& mc) {
  // generated electron (17 variables)
  tree->Branch("mc_Q2", &mc.Q2);
  //tree->Branch("mc_W", &mc.W);
  tree->Branch("mc_Nu", &mc.Nu);
  tree->Branch("mc_Xb", &mc.Xb);
  //tree->Branch("mc_Yb", &mc.Yb);
  tree->Branch("mc_vxe", &mc.vxe);
  tree->Branch("mc_vye", &mc.vye);
  tree->Branch("mc_vze", &mc.vze);
  tree->Branch("mc_SectorEl", &mc.SectorEl);
  tree->Branch("mc_TargType", &mc.TargType);
  tree->Branch("mc_TargTypeSM", &mc.TargTypeSM);
  tree->Branch("mc_TargTypeHH", &mc.TargTypeHH);
  tree->Branch("mc_Pex", &mc.Pex);
  tree->Branch("mc_Pey", &mc.Pey);
  tree->Branch("mc_Pez", &mc.Pez);
  tree->Branch("mc_Pe", &mc.Pe);
  tree->Branch("mc_BettaEl", &mc.BettaEl);
  tree->Branch("mc_ThetaLabEl", &mc.ThetaLabEl);
  //tree->Branch("mc_PhiLabEl", &mc.PhiLabEl);
  // generated particle (22 variables)
  //tree->Branch("mc_Eh", &mc.Eh);
  tree->Branch("mc_Zh", &mc.Zh);
  tree->Branch("mc_ThetaPQ", &mc.ThetaPQ);
  tree->Branch("mc_Pt2", &mc.Pt2);
  //tree->Branch("mc_Pl2", &mc.Pl2);
  tree->Branch("mc_PhiPQ", &mc.PhiPQ);
  //tree->Branch("mc_Mx2", &mc.Mx2);
  //tree->Branch("mc_T", &mc.T);
  tree->Branch("mc_ThetaLab", &mc.ThetaLab);
  tree->Branch("mc_PhiLab", &mc.PhiLab);
  tree->Branch("mc_vxh", &mc.vxh);
  tree->Branch("mc_vyh", &mc.vyh);
  tree->Branch("mc_vzh", &mc.vzh);
  tree->Branch("mc_Sector", &mc.Sector);
  tree->Branch("mc_Px", &mc.Px);
  tree->Branch("mc_Py", &mc.Py);
  tree->Branch("mc_Pz", &mc.Pz);
  tree->Branch("mc_P", &mc.P);
  tree->Branch("mc_Betta", &mc.Betta);
  //tree->Branch("mc_Mass2", &mc.Mass2);
  tree->Branch("mc_pid", &mc.pid);
  tree->Branch("mc_deltaZ", &mc.deltaZ);
}

void AssignElectronVar_REC(TIdentificatorV2* t, rec_p& rec, Int_t evnt) {
  // reconstructed electron (46 variables)
  rec.Q2 = t->Q2(0);
  rec.W  = t->W(0);
  rec.Nu = t->Nu(0);
  rec.Xb = t->Xb(0);
  rec.Yb = t->Yb(0);
  rec.vxe = t->X(0);
  rec.vye = t->Y(0);
  rec.vze = t->Z(0);
  rec.SectorEl = t->Sector(0);
  rec.Pex = t->Px(0);
  rec.Pey = t->Py(0);
  rec.Pez = t->Pz(0);
  rec.Pe = t->Momentum(0);
  TVector3* fVert = t->GetCorrectedVert();
  rec.vyec = fVert->Y();
  rec.vzec = fVert->Z();
  rec.TargType   = t->TargType(gDataKind, gTargetOption);
  rec.TargTypeSM = t->TargTypeSM(gDataKind, gTargetOption, 0);
  rec.TargTypeHH = t->TargTypeHH(0);
  rec.FidCheckCutEl = t->FidCheckCut(0);            // added by smoran
  rec.NpheEl = t->Nphe(0);
  rec.CCMirrorMatchingEl = t->CCMirrorMatching(0);
  rec.SampFractionEl25   = t->SampFracCheck(0, gDataKind, gTargetOption, 2.5);
  rec.SampFractionEl20   = t->SampFracCheck(0, gDataKind, gTargetOption, 2.0);
  rec.PhiLabEl = t->PhiLab(0);

  //rec.BettaEl = t->Betta(0);
  //rec.Etote = t->Etot(0);
  //rec.Eine = t->Ein(0);
  //rec.Eoute = t->Eout(0);
  //rec.vxec = fVert->X();
  //rec.XECe = t->XEC(0);
  //rec.YECe = t->YEC(0);
  //rec.ZECe = t->ZEC(0);
  //rec.ThetaLabEl = t->ThetaLab(0);
  //rec.StatDCEl = t->StatDC(0);
  //rec.DCStatusEl = t->DCStatus(0);
  //rec.StatECEl = t->StatEC(0);
  //rec.ECStatusEl = t->ECStatus(0);
  //rec.TimeECEl = t->TimeEC(0);
  //rec.PathECEl = t->PathEC(0);
  //rec.Chi2ECEl = t->Chi2EC(0);
  //rec.StatSCEl = t->StatSC(0);
  //rec.SCStatusEl = t->SCStatus(0);
  //rec.TimeSCEl = t->TimeSC(0);
  //rec.PathSCEl = t->PathSC(0);
  //rec.StatCCEl = t->StatCC(0);
  //rec.CCStatusEl = t->CCStatus(0);
  //rec.Chi2CCEl = t->Chi2CC(0);
  //rec.StatusEl = t->Status(0);
  //rec.NRowsDCEl = t->NRowsDC();
  //rec.NRowsECEl = t->NRowsEC();
  //rec.NRowsSCEl = t->NRowsSC();
  //rec.NRowsCCEl = t->NRowsCC();
  //event-related (1 variable)
  
  rec.evnt = evnt;
}

void AssignParticleVar_REC(TIdentificatorV2* t, rec_p& rec, Int_t row) {
  // reconstructed particles (49 variables)
  Int_t fPid = particleID(t->GetCategorization(row, gDataKind, gTargetOption));
  rec.pid.push_back(fPid);
  rec.vxh.push_back(t->X(row));
  rec.vyh.push_back(t->Y(row));
  rec.vzh.push_back(t->Z(row));
  rec.Sector.push_back(t->Sector(row));
  TVector3* fVert = t->GetCorrectedVert();
  rec.deltaZ.push_back(t->Z(row) - fVert->Z());
  rec.XEC.push_back(t->XEC(row));
  rec.YEC.push_back(t->YEC(row));
  rec.ZEC.push_back(t->ZEC(row));
  rec.Etot.push_back(t->Etot(row));
  rec.Ein.push_back(t->Ein(row));
  rec.Eout.push_back(t->Eout(row));
  Double_t fMass = particleMass(fPid);
  //mass of pion (so that T4 assuming pion, ie T4_pi, can be added to the tuple)
  Double_t fMassPi = particleMass(211);
  TLorentzVector* fGamma = t->GetCorrPhotonMomentum(row);
  Double_t fPx = (fPid == 22) * fGamma->Px() + (fPid != 22) * t->Px(row);
  Double_t fPy = (fPid == 22) * fGamma->Py() + (fPid != 22) * t->Py(row);
  Double_t fPz = (fPid == 22) * fGamma->Pz() + (fPid != 22) * t->Pz(row);
  Double_t fP = TMath::Sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
  Double_t fE = (fPid == 22) * fGamma->E() + (fPid != 22) * TMath::Sqrt(fMass * fMass + fP * fP);
  Double_t fZ = fE / t->Nu();
  rec.Px.push_back(fPx);
  rec.Py.push_back(fPy);
  rec.Pz.push_back(fPz);
  rec.P.push_back(fP);
  rec.PhiLab.push_back(PhiLab(fPx, fPy, fPz));
  rec.ThetaLab.push_back(ThetaLab(fPx, fPy, fPz));
  rec.Zh.push_back(fZ);
  rec.ThetaPQ.push_back(ThetaPQ(t->Px(0), t->Py(0), t->Pz(0), fPx, fPy, fPz));
  rec.PhiPQ.push_back(PhiPQ(t->Px(0), t->Py(0), t->Pz(0), fPx, fPy, fPz));
  Double_t fCosThetaPQ = (fPz * (kEbeam - t->Pz(0)) - fPx * t->Px(0) - fPy * t->Py(0)) / (TMath::Sqrt(t->Nu() * t->Nu() + t->Q2()) * fP);
  Double_t fPt2 = fP * fP * (1 - fCosThetaPQ * fCosThetaPQ);
  Double_t fPl2 = fP * fP * fCosThetaPQ * fCosThetaPQ;
  rec.Pt2.push_back(fPt2);
  rec.Pl2.push_back(fPl2);
  rec.Mx2.push_back(t->W() * t->W() + fMass * fMass - 2 * fZ * t->Nu() * t->Nu() + 2 * TMath::Sqrt(fPl2 * (t->Nu() * t->Nu() + t->Q2())) - 2 * kMassProton * fZ * t->Nu());
  rec.T.push_back(fMass * fMass - 2 * fZ * t->Nu() * t->Nu() + 2 * TMath::Sqrt(fPl2 * (t->Nu() * t->Nu() + t->Q2())) - t->Q2());
  rec.Betta.push_back(t->Betta(row));  // BettaMeasured
  rec.Mass2.push_back(fP * fP * (TMath::Power(t->Betta(row), -2) - 1));
  rec.T4.push_back(t->PathSC(0) / 30. - t->TimeSC(0) + t->TimeSC(row) - (t->PathSC(row) / 30.) * TMath::Sqrt(TMath::Power(fMass / fP, 2) + 1));
  rec.T4_pi.push_back(t->PathSC(0) / 30. - t->TimeSC(0) + t->TimeSC(row) - (t->PathSC(row) / 30.) * TMath::Sqrt(TMath::Power(fMassPi / fP, 2) + 1));
  rec.FidCheckCutPiPlus.push_back(t->FidCheckCutPiPlus(row)); // added by smoran
  rec.Nphe.push_back(t->Nphe(row));
  rec.Chi2CC.push_back(t->Chi2CC(row));
  rec.NRowsCC.push_back(t->NRowsCC());
  rec.StatCC.push_back(t->StatCC(row));
  rec.FidCheckCut.push_back(t->FidCheckCut(row)); // added by smoran
  rec.SampFraction25.push_back(t->SampFracCheck(row, gDataKind, gTargetOption, 2.5));
  rec.SampFraction20.push_back(t->SampFracCheck(row, gDataKind, gTargetOption, 2.0));
  rec.CCMirrorMatching.push_back(t->CCMirrorMatching(row));
  rec.Eh.push_back(fE);

  
  //rec.StatDC.push_back(t->StatDC(row));
  //rec.DCStatus.push_back(t->DCStatus(row));
  //rec.StatEC.push_back(t->StatEC(row));
  //rec.ECStatus.push_back(t->ECStatus(row));
  //rec.TimeEC.push_back(t->TimeEC(row));
  //rec.PathEC.push_back(t->PathEC(row));
  //rec.Chi2EC.push_back(t->Chi2EC(row));
  //rec.StatSC.push_back(t->StatSC(row));
  //rec.SCStatus.push_back(t->SCStatus(row));
  //rec.TimeSC.push_back(t->TimeSC(row));
  //rec.PathSC.push_back(t->PathSC(row));
  //rec.CCStatus.push_back(t->CCStatus(row));
  //rec.Status.push_back(t->Status(row));
  //rec.NRowsDC.push_back(t->NRowsDC());
  //rec.NRowsEC.push_back(t->NRowsEC());
  //rec.NRowsSC.push_back(t->NRowsSC());
  
}

void ClearParticleVar_REC(rec_p& rec) {
  // reconstructed particles (49 variables)
  rec.pid.clear();
  rec.vxh.clear();
  rec.vyh.clear();
  rec.vzh.clear();
  rec.Sector.clear();
  rec.deltaZ.clear();
  rec.XEC.clear();
  rec.YEC.clear();
  rec.ZEC.clear();
  rec.Etot.clear();
  rec.Ein.clear();
  rec.Eout.clear();
  rec.Px.clear();
  rec.Py.clear();
  rec.Pz.clear();
  rec.P.clear();
  rec.PhiLab.clear();
  rec.ThetaLab.clear();
  rec.Zh.clear();
  rec.ThetaPQ.clear();
  rec.PhiPQ.clear();
  rec.Pt2.clear();
  rec.Pl2.clear();
  rec.Mx2.clear();
  rec.T.clear();
  rec.Betta.clear();
  rec.Mass2.clear();
  rec.T4.clear();
  rec.T4_pi.clear();
  rec.StatCC.clear();
  rec.FidCheckCutPiPlus.clear();
  rec.Nphe.clear();
  rec.Chi2CC.clear();
  rec.NRowsCC.clear();
  rec.FidCheckCut.clear();
  rec.SampFraction25.clear();
  rec.SampFraction20.clear();
  rec.CCMirrorMatching.clear();
  rec.Eh.clear();
  
  //rec.StatDC.clear();
  //rec.DCStatus.clear();
  //rec.StatEC.clear();
  //rec.ECStatus.clear();
  //rec.TimeEC.clear();
  //rec.PathEC.clear();
  //rec.Chi2EC.clear();
  //rec.StatSC.clear();
  //rec.SCStatus.clear();
  //rec.TimeSC.clear();
  //rec.PathSC.clear();
  //rec.CCStatus.clear();
  //rec.Status.clear();
  //rec.NRowsDC.clear();
  //rec.NRowsEC.clear();
  //rec.NRowsSC.clear();
  
}

void AssignElectronVar_GEN(TIdentificatorV2* t, gen_p& mc) {
  // generated electron (17 variables)
  mc.Q2 = t->Q2(1);
  mc.W = t->W(1);
  mc.Nu = t->Nu(1);
  mc.Xb = t->Xb(1);
  //mc.Yb = t->Yb(1);
  mc.vxe = t->X(0, 1);
  mc.vye = t->Y(0, 1);
  mc.vze = t->Z(0, 1);
  mc.SectorEl = t->Sector(0, 1);
  mc.Pex = t->Px(0, 1);
  mc.Pey = t->Py(0, 1);
  mc.Pez = t->Pz(0, 1);
  mc.Pe = t->Momentum(0, 1);
  mc.BettaEl = t->Betta(0, 1);
  mc.ThetaLabEl = t->ThetaLab(0, 1);
  mc.PhiLabEl = t->PhiLab(0, 1);
  mc.TargType = t->TargTypeSM(gDataKind, gTargetOption, 1);

}

void AssignParticleVar_GEN(TIdentificatorV2* t, gen_p& mc, Int_t row) {
  // generated particles (22 variables)
  Int_t fPid = ToPDG(t->Id(row, 1));  // from GEANT to PDG
  mc.pid.push_back(fPid);
  mc.ThetaPQ.push_back(t->ThetaPQ(row, 1));
  mc.PhiPQ.push_back(t->PhiPQ(row, 1));
  mc.Pt2.push_back(t->Pt2(row, 1));
  mc.Pl2.push_back(t->Pl2(row, 1));
  Double_t fMass = particleMass(fPid);
  Double_t fE = t->Zh(row, 1, fMass) * t->Nu(1);
  mc.Zh.push_back(t->Zh(row, 1, fMass));
  //mc.Eh.push_back(fE);
  mc.Mx2.push_back(t->Mx2(row, 1, fMass));
  //mc.T.push_back(t->T(row, 1, fMass));
  mc.ThetaLab.push_back(t->ThetaLab(row, 1));
  mc.PhiLab.push_back(t->PhiLab(row, 1));
  mc.vxh.push_back(t->X(row, 1));
  mc.vyh.push_back(t->Y(row, 1));
  mc.vzh.push_back(t->Z(row, 1));
  mc.Sector.push_back(t->Sector(row, 1));
  mc.Px.push_back(t->Px(row, 1));
  mc.Py.push_back(t->Py(row, 1));
  mc.Pz.push_back(t->Pz(row, 1));
  mc.P.push_back(t->Momentum(row, 1));
  mc.Betta.push_back(t->Betta(row, 1));
  //mc.Mass2.push_back(t->Mass2(row, 1));
  mc.deltaZ.push_back(t->Z(row, 1) - t->Z(0, 1));
}

void ClearParticleVar_GEN(gen_p& mc) {
  // generated particle (22 variables)
  mc.pid.clear();
  mc.ThetaPQ.clear();
  mc.PhiPQ.clear();
  mc.Pt2.clear();
  mc.Pl2.clear();
  mc.Zh.clear();
  //mc.Eh.clear();
  mc.Mx2.clear();
  mc.T.clear();
  mc.ThetaLab.clear();
  mc.PhiLab.clear();
  mc.vxh.clear();
  mc.vyh.clear();
  mc.vzh.clear();
  mc.Sector.clear();
  mc.Px.clear();
  mc.Py.clear();
  mc.Pz.clear();
  mc.P.clear();
  mc.Betta.clear();
  mc.Mass2.clear();
  mc.deltaZ.clear();
}

void NullElectronVar_REC(rec_p& rec) {
  // reconstructed electron (46 variables)
  rec.Q2 = INVLD;
  rec.W = INVLD;
  rec.Nu = INVLD;
  rec.Xb = INVLD;
  rec.Yb = INVLD;
  rec.vxe = INVLD;
  rec.vye = INVLD;
  rec.vze = INVLD;
  rec.SectorEl = INVLD;
  rec.Pex = INVLD;
  rec.Pey = INVLD;
  rec.Pez = INVLD;
  rec.Pe = INVLD;
  rec.BettaEl = INVLD;
  rec.Etote = INVLD;
  rec.Eine = INVLD;
  rec.Eoute = INVLD;
  rec.vxec = INVLD;
  rec.vyec = INVLD;
  rec.vzec = INVLD;
  rec.TargType = INVLD;
  rec.XECe = INVLD;
  rec.YECe = INVLD;
  rec.ZECe = INVLD;
  rec.ThetaLabEl = INVLD;
  rec.PhiLabEl = INVLD;
  rec.StatDCEl = INVLD;
  rec.DCStatusEl = INVLD;
  rec.StatECEl = INVLD;
  rec.ECStatusEl = INVLD;
  rec.TimeECEl = INVLD;
  rec.PathECEl = INVLD;
  rec.Chi2ECEl = INVLD;
  rec.StatSCEl = INVLD;
  rec.SCStatusEl = INVLD;
  rec.TimeSCEl = INVLD;
  rec.PathSCEl = INVLD;
  rec.StatCCEl = INVLD;
  rec.CCStatusEl = INVLD;
  rec.NpheEl = INVLD;
  rec.CCMirrorMatchingEl = INVLD;
  rec.SampFractionEl25= INVLD;
  rec.SampFractionEl20= INVLD;
  rec.Chi2CCEl = INVLD;
  rec.StatusEl = INVLD;
  rec.NRowsDCEl = INVLD;
  rec.NRowsECEl = INVLD;
  rec.NRowsSCEl = INVLD;
  rec.NRowsCCEl = INVLD;
}

void NullParticleVar_REC(rec_p& rec) {
  // reconstructed particles (49 variables)
  rec.pid.push_back(INVLD);
  rec.Eh.push_back(INVLD);
  rec.Zh.push_back(INVLD);
  rec.ThetaPQ.push_back(INVLD);
  rec.PhiPQ.push_back(INVLD);
  rec.Pt2.push_back(INVLD);
  rec.Pl2.push_back(INVLD);
  rec.Mx2.push_back(INVLD);
  rec.T.push_back(INVLD);
  rec.ThetaLab.push_back(INVLD);
  rec.PhiLab.push_back(INVLD);
  rec.T4.push_back(INVLD);
  rec.T4_pi.push_back(INVLD);
  rec.vxh.push_back(INVLD);
  rec.vyh.push_back(INVLD);
  rec.vzh.push_back(INVLD);
  rec.Sector.push_back(INVLD);
  rec.Px.push_back(INVLD);
  rec.Py.push_back(INVLD);
  rec.Pz.push_back(INVLD);
  rec.P.push_back(INVLD);
  rec.Betta.push_back(INVLD);
  rec.Mass2.push_back(INVLD);
  rec.Etot.push_back(INVLD);
  rec.Ein.push_back(INVLD);
  rec.Eout.push_back(INVLD);
  rec.XEC.push_back(INVLD);
  rec.YEC.push_back(INVLD);
  rec.ZEC.push_back(INVLD);
  rec.deltaZ.push_back(INVLD);
  rec.StatCC.push_back(INVLD);
  rec.NRowsCC.push_back(INVLD);
  rec.Nphe.push_back(INVLD);
  rec.Chi2CC.push_back(INVLD);
  rec.FidCheckCutPiPlus.push_back(INVLD);
  rec.FidCheckCut.push_back(INVLD);
  
  
  //rec.StatDC.push_back(INVLD);
  //rec.DCStatus.push_back(INVLD);
  //rec.StatEC.push_back(INVLD);
  //rec.ECStatus.push_back(INVLD);
  //rec.TimeEC.push_back(INVLD);
  //rec.PathEC.push_back(INVLD);
  //rec.Chi2EC.push_back(INVLD);
  //rec.StatSC.push_back(INVLD);
  //rec.SCStatus.push_back(INVLD);
  //rec.TimeSC.push_back(INVLD);
  //rec.PathSC.push_back(INVLD);
  //rec.CCStatus.push_back(INVLD);
  //rec.Status.push_back(INVLD);
  //rec.NRowsDC.push_back(INVLD);
  //rec.NRowsEC.push_back(INVLD);
  //rec.NRowsSC.push_back(INVLD);

}

void NullElectronVar_GEN(gen_p& mc) {
  // gsim electron (17 variables)
  mc.Q2 = INVLD;
  mc.W = INVLD;
  mc.Nu = INVLD;
  mc.Xb = INVLD;
  mc.Yb = INVLD;
  mc.vxe = INVLD;
  mc.vye = INVLD;
  mc.vze = INVLD;
  mc.SectorEl = INVLD;
  mc.TargType = INVLD;
  mc.Pex = INVLD;
  mc.Pey = INVLD;
  mc.Pez = INVLD;
  mc.Pe = INVLD;
  mc.BettaEl = INVLD;
  mc.ThetaLabEl = INVLD;
  mc.PhiLabEl = INVLD;
}

void NullParticleVar_GEN(gen_p& mc) {
  // gsim particle (22 variables)
  //mc.Eh.push_back(INVLD);
  mc.Zh.push_back(INVLD);
  mc.ThetaPQ.push_back(INVLD);
  mc.Pt2.push_back(INVLD);
  mc.Pl2.push_back(INVLD);
  mc.PhiPQ.push_back(INVLD);
  mc.ThetaLab.push_back(INVLD);
  mc.PhiLab.push_back(INVLD);
  mc.Mx2.push_back(INVLD);
  mc.T.push_back(INVLD);
  mc.vxh.push_back(INVLD);
  mc.vyh.push_back(INVLD);
  mc.vzh.push_back(INVLD);
  mc.Sector.push_back(INVLD);
  mc.Px.push_back(INVLD);
  mc.Py.push_back(INVLD);
  mc.Pz.push_back(INVLD);
  mc.P.push_back(INVLD);
  mc.Betta.push_back(INVLD);
  mc.Mass2.push_back(INVLD);
  mc.deltaZ.push_back(INVLD);
  mc.pid.push_back(INVLD);
}

#endif
