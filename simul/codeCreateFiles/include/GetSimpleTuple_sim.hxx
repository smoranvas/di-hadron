#include "Headers.hxx"

#include "GSTtree.h"
#include "TreeOps.hxx"
#include "UI.hxx"

using namespace ROOT::VecOps;

void FindParticles(TClasTool* input, TIdentificatorV2* t, RVec<Int_t>& gsim_row, RVec<Int_t>& simrec_row) {
  if (input->GetNRows("GSIM") > 0) {  // prevent seg-fault
    // first, check numbering scheme
    if (t->Id(0, 1) == 11)
      SetNumberingScheme("PDG");
    else if (t->Id(0, 1) == 3)
      SetNumberingScheme("GEANT");

    if (t->Id(0, 1) == gElectronID) {
      for (Int_t q = 1; q < input->GetNRows("GSIM"); q++) {
        if (t->Id(q, 1) == gPiPlusID || 
	    t->Id(q, 1) == gPiMinusID || 
	    t->Id(q, 1) == gProtonID ||                                                                                        
	    t->Id(q, 1) == gKaonPlusID ||                                                                                       
	    t->Id(q, 1) == gKaonMinusID ||                                                                                      
	    t->Id(q, 1) == gElectronID ||                                                                                  
	    t->Id(q, 1) == gPositronID 
	 /*
	   || t->Id(q, 1) == gGammaID || 
	      t->Id(q, 1) == gElectronID || 
	      t->Id(q, 1) == gPositronID ||
              t->Id(q, 1) == gProtonID || 
	      t->Id(q, 1) == gNeutronID || 
	      t->Id(q, 1) == gKaonPlusID || 
	      t->Id(q, 1) == gKaonMinusID || 
	      t->Id(q, 1) == gKaonZeroLongID ||
              t->Id(q, 1) == gKaonZeroShortID || 
	      t->Id(q, 1) == gKaonZeroID
	 */   ) {
          gsim_row.push_back(q);
        }
      }  // end of loop in gsim-particles

      if (input->GetNRows("EVNT") > 0) {  // prevent seg-fault
        if (t->GetCategorization(0, gDataKind, gTargetOption) == "electron") {
          for (Int_t p = 1; p < input->GetNRows("EVNT"); p++) {
            if (t->GetCategorization(p, gDataKind, gTargetOption) == "pi+" || 
	        t->GetCategorization(p, gDataKind, gTargetOption) == "pi-" ||
                t->GetCategorization(p, gDataKind, gTargetOption) == "proton") {
              simrec_row.push_back(p);
            }
          }  // end of loop in simrec-particles
        }    // end of electron-in-simrec condition
      }      // end of smth-in-EVNT-bank

    }  // end of electron-in-gsim condition
  }    // end of smth-in-GSIM-bank
}

RVec<Int_t> SortByMomentum(TIdentificatorV2* t, RVec<Int_t> vector_row, Int_t kind) {
  // Returns a new "vector_row2" where particle indices are sorted by their respective momentum (from lower to higher)
  // first, fill the momentum vector
  RVec<Double_t> momentum;
  for (Size_t m = 0; m < vector_row.size(); m++) {
    momentum.push_back(t->Momentum(vector_row[m], kind));
  }
  // Argsort() creates an indices-vector with the indices sorted by the input-vector values
  RVec<Int_t> indices = Argsort(momentum);
  // Take(input-vector, indices-vector) creates a sorted-vector by moving all input-vector indices to match the order assigned by the indices-vector
  RVec<Int_t> vector_row2 = Take(vector_row, indices);
  return vector_row2;
}

void AngularMatching(TIdentificatorV2* t, RVec<Int_t>& simrec_row, RVec<Int_t>& gsim_row) {
  // Matches the "simrec_row" vector with the "gsim_row" vector under angular matching
  // if particles don't match, the counterpart is filled with null

  // worst from electron (CLAS paper @ P = 0.75 GeV):
  //       fDeltaThetaLab = 0.33, fDeltaPhiLab = 0.96
  // arbitrary values:
  //       fDeltaThetaLab = 0.5, fDeltaPhiLab = 2.0
  // measured values: (worst from pi+ with 0 < P < 0.35 GeV)
  //       fDeltaThetaLab = 0.78, fDeltaPhiLab = 1.41
  // worst from electron (CLAS paper @ P = 0.1 GeV):
  //const Double_t fDeltaThetaLab = 2.40;  // Delta_Theta = 3*sigma_Theta
  //const Double_t fDeltaPhiLab = 5.43;    // Delta_Phi = 3*sigma_Phi
  const Double_t fDeltaThetaLab = 2.40*5/3.;  // Delta_Theta = 5*sigma_Theta                                                    
  const Double_t fDeltaPhiLab = 5.43*5/3.;    // Delta_Phi = 5*sigma_Phi 


  // define output vectors - initially empty
  RVec<Int_t> simrec_new;
  RVec<Int_t> gsim_new;

  // declare PhiLab and ThetaLab
  Double_t simrec_phi;
  Double_t simrec_theta;
  Double_t gsim_phi;
  Double_t gsim_theta;

  // m, n are the vectors' indices
  for (Size_t m = 0; m < gsim_row.size(); m++) {
    for (Size_t n = 0; n < simrec_row.size(); n++) {

      // gsim angles don't need correction
      gsim_phi = PhiLab(t->Px(gsim_row[m], 1), t->Py(gsim_row[m], 1), t->Pz(gsim_row[m], 1));
      gsim_theta = ThetaLab(t->Px(gsim_row[m], 1), t->Py(gsim_row[m], 1), t->Pz(gsim_row[m], 1));

      // for reconstructed photons
      TLorentzVector* fGamma = t->GetCorrPhotonMomentum(simrec_row[n]);
      TString fParticleName = t->GetCategorization(simrec_row[n], gDataKind, gTargetOption);
      Double_t fPx = (fParticleName == "gamma") * fGamma->Px() + (fParticleName != "gamma") * t->Px(simrec_row[n], 0);
      Double_t fPy = (fParticleName == "gamma") * fGamma->Py() + (fParticleName != "gamma") * t->Py(simrec_row[n], 0);
      Double_t fPz = (fParticleName == "gamma") * fGamma->Pz() + (fParticleName != "gamma") * t->Pz(simrec_row[n], 0);

      simrec_phi = PhiLab(fPx, fPy, fPz);
      simrec_theta = ThetaLab(fPx, fPy, fPz);

      /*** MATCHING CONDITION ***/
      double dphi =simrec_phi - gsim_phi;
      if(dphi > 180)
        dphi -= 360;
      if(dphi < -180)
	dphi += 360;
      Bool_t fAngularMatching = TMath::Abs(dphi) < fDeltaPhiLab && TMath::Abs(simrec_theta - gsim_theta) < fDeltaThetaLab;
      if (fAngularMatching && std::find(gsim_new.begin(), gsim_new.end(), gsim_row[m]) == gsim_new.end() &&
          std::find(simrec_new.begin(), simrec_new.end(), simrec_row[n]) == simrec_new.end()) {
        // std::find function returns an iterator to the first element in the range ["begin","end"[ that compares equal to "row"
        // if no such element is found, the function returns "end", this is done to prevent particles that already matched
        gsim_new.push_back(gsim_row[m]);      // add gsim_row to gsim_new
        simrec_new.push_back(simrec_row[n]);  // add simrec_row to simrec_new
      }

    }  // end of simrec loop
  }    // end of gsim loop

  /*** FILL REMAINING ***/

  for (Size_t m = 0; m < gsim_row.size(); m++) {
    // if it's not in gsim_new
    if (std::find(gsim_new.begin(), gsim_new.end(), gsim_row[m]) == gsim_new.end()) {
      gsim_new.push_back(gsim_row[m]);  // add gsim_row to gsim_new
      simrec_new.push_back(-9999);         // add null simrec
    }
  }

  for (Size_t n = 0; n < simrec_row.size(); n++) {
    // if it's not in simrec_new
    if (std::find(simrec_new.begin(), simrec_new.end(), simrec_row[n]) == simrec_new.end()) {
      simrec_new.push_back(simrec_row[n]);  // add simrec_row to simrec_new
      gsim_new.push_back(-9999);               // add null gsim
    }
  }

  // assign results
  gsim_row = gsim_new;
  simrec_row = simrec_new;
}

void FillElectronNParticles(TIdentificatorV2* t, TClasTool* input, TTree* tParticles, RVec<Int_t>& gsim_row, RVec<Int_t>& simrec_row, gen_p& mc, rec_p& rec, Int_t i) {

  /*** Electron ***/

  if (input->GetNRows("GSIM") > 0) {  // prevent seg-fault
    if (t->Id(0, 1) == gElectronID) {
      AssignElectronVar_GEN(t, mc);       // (TIdentificatorV2, gen_e)
      rec.evnt = i;
      if (input->GetNRows("EVNT") > 0) {  // prevent seg-fault
        if (t->GetCategorization(0, gDataKind, gTargetOption) != "electron")
          NullElectronVar_REC(rec);
        else
          AssignElectronVar_REC(t, rec, i);
      }  // end of smth-in-EVNT-bank condition
    }    // end of electron-in-GSIM condition
  }      // end of smth-in-GSIM-bank condition

  /*** Particles ***/

  for (Size_t index = 0; index < simrec_row.size(); index++) {  // simrec_row.size() == gsim_row.size()
    // gsim
    if (gsim_row[index] == -9999)
      NullParticleVar_GEN(mc);
    else
      AssignParticleVar_GEN(t, mc, gsim_row[index]);  // (TIdentificatorV2, gen_p, row)
    // simrec
    if (simrec_row[index] == -9999)
      NullParticleVar_REC(rec);
    else
      AssignParticleVar_REC(t, rec, simrec_row[index]);  // (TIdentificatorV2, rec_p, row)
  }

  /*** Fill! ***/

  tParticles->Fill();
  ClearParticleVar_GEN(mc);
  ClearParticleVar_REC(rec);  
}
