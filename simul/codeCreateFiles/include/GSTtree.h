#ifndef GSTTREE_H
#define GSTTREE_H

#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

/*
   NAMING CONVENTIONS

   rec_p : output variables for reconstructed particles (data or simulations)
   gen_p : output variables for generated particles (simulations)
*/

struct rec_p {
  // electron (46)
  Float_t Q2;          // 
  Float_t W;           //
  Float_t Nu;          // virtual photon's energy, Nu = (electron's initial energy) - (electron's final energy), in (GeV)
  Float_t Xb;          //
  Float_t Yb;          //
  Float_t vxe;         // x-vertex of electron, in (cm)
  Float_t vye;         // y-vertex of electron, in (cm)
  Float_t vze;         // z-vertex of electron, in (cm)
  Int_t SectorEl;      //
  Float_t Pex;         // x-component of electron's momentum, in (GeV)
  Float_t Pey;         // y-component of electron's momentum, in (GeV)
  Float_t Pez;         // z-component of electron's momentum, in (GeV)
  Float_t Pe;          // electron's momentum, in (GeV)
  Float_t BettaEl;     //
  Float_t Etote;       //
  Float_t Eine;        //
  Float_t Eoute;       //
  Float_t vxec;        // corrected x-vertex of electron, in (cm)
  Float_t vyec;        // corrected y-vertex of electron, in (cm)
  Float_t vzec;        // corrected z-vertex of electron, in (cm)
  Int_t TargType;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargType()
  Float_t XECe;        //
  Float_t YECe;        //
  Float_t ZECe;        //
  Float_t ThetaLabEl;  //
  Float_t PhiLabEl;    //
  Float_t StatDCEl;    //
  Float_t DCStatusEl;  //
  Float_t StatECEl;    //
  Float_t ECStatusEl;  //
  Float_t TimeECEl;    //
  Float_t PathECEl;    //
  Float_t Chi2ECEl;    //
  Float_t StatSCEl;    //
  Float_t SCStatusEl;  //
  Float_t TimeSCEl;    //
  Float_t PathSCEl;    //
  Float_t StatCCEl;    //
  Float_t CCStatusEl;  //
  Float_t NpheEl;      //
  Float_t Chi2CCEl;    //
  Float_t StatusEl;    //
  Float_t NRowsDCEl;   //
  Float_t NRowsECEl;   //
  Float_t NRowsSCEl;   //
  Float_t NRowsCCEl;   //
  Bool_t FidCheckCutEl;     // added by smoran
  Int_t TargTypeSM;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargTypeSM()
  Int_t TargTypeHH;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargTypeHH()
  Bool_t CCMirrorMatchingEl;
  Bool_t SampFractionEl25;
  Bool_t SampFractionEl20;

  // particle (49)
  std::vector<Float_t> Eh;        // energy of the outgoing particle, in (GeV)
  std::vector<Float_t> Zh;        // fraction energy of the outgoing particle, Zh = Eh/Nu
  std::vector<Float_t> ThetaPQ;   //
  std::vector<Float_t> Pt2;       //
  std::vector<Float_t> Pl2;       //
  std::vector<Float_t> PhiPQ;     //
  std::vector<Float_t> Mx2;       //
  std::vector<Float_t> T;         //
  std::vector<Float_t> ThetaLab;  //
  std::vector<Float_t> PhiLab;    //
  std::vector<Float_t> vxh;       // x-vertex of particle, in (cm)
  std::vector<Float_t> vyh;       // y-vertex of particle, in (cm)
  std::vector<Float_t> vzh;       // z-vertex of particle, in (cm)
  std::vector<Int_t> Sector;      //
  std::vector<Float_t> Px;        // x-component of particle's momentum, in (GeV)
  std::vector<Float_t> Py;        // y-component of particle's momentum, in (GeV)
  std::vector<Float_t> Pz;        // z-component of particle's momentum, in (GeV)
  std::vector<Float_t> P;         // particle's momentum, in (GeV)
  std::vector<Float_t> Betta;     //
  std::vector<Float_t> Mass2;     //
  std::vector<Float_t> Etot;      //
  std::vector<Float_t> Ein;       //
  std::vector<Float_t> Eout;      //
  std::vector<Float_t> XEC;       //
  std::vector<Float_t> YEC;       //
  std::vector<Float_t> ZEC;       //
  std::vector<Int_t> pid;         // particle ID number (according to PDG)
  std::vector<Float_t> T4;        //
  std::vector<Float_t> T4_pi;     // T4, assuming the mass of a charged pion
  std::vector<Float_t> deltaZ;    //
  std::vector<Float_t> StatDC;    //
  std::vector<Float_t> DCStatus;  //
  std::vector<Float_t> StatEC;    //
  std::vector<Float_t> ECStatus;  //
  std::vector<Float_t> TimeEC;    //
  std::vector<Float_t> PathEC;    //
  std::vector<Float_t> Chi2EC;    //
  std::vector<Float_t> StatSC;    //
  std::vector<Float_t> SCStatus;  //
  std::vector<Float_t> TimeSC;    //
  std::vector<Float_t> PathSC;    //
  std::vector<Float_t> StatCC;    //
  std::vector<Float_t> CCStatus;  //
  std::vector<Float_t> Nphe;      //
  std::vector<Float_t> Chi2CC;    //
  std::vector<Float_t> Status;    //
  std::vector<Float_t> NRowsDC;   //
  std::vector<Float_t> NRowsEC;   //
  std::vector<Float_t> NRowsSC;   //
  std::vector<Float_t> NRowsCC;   //
  std::vector<Bool_t>  FidCheckCutPiPlus;  // Added by smoran
  std::vector<Bool_t>  FidCheckCut;  // Added by smoran
  std::vector<Bool_t>  SampFraction25;
  std::vector<Bool_t>  SampFraction20;
  std::vector<Bool_t>  CCMirrorMatching;
  // event (1)
  Float_t evnt;  // event number  

};

struct gen_p {
  // electron (17)
  Float_t Q2;          // 
  Float_t W;           //
  Float_t Nu;          // virtual photon's energy, Nu = (electron's initial energy) - (electron's final energy), in (GeV)
  Float_t Xb;          //
  Float_t Yb;          //
  Float_t vxe;         // x-vertex of electron, in (cm)
  Float_t vye;         // y-vertex of electron, in (cm)
  Float_t vze;         // z-vertex of electron, in (cm)
  Int_t SectorEl;      //
  Int_t TargType;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargType()
  Float_t Pex;         // x-component of electron's momentum, in (GeV)
  Float_t Pey;         // y-component of electron's momentum, in (GeV)
  Float_t Pez;         // z-component of electron's momentum, in (GeV)
  Float_t Pe;          // electron's momentum, in (GeV)
  Float_t BettaEl;     //
  Float_t ThetaLabEl;  //
  Float_t PhiLabEl;    //
  // particle (22)
  std::vector<Float_t> Eh;        // energy of the outgoing particle, in (GeV)
  std::vector<Float_t> Zh;        // fraction energy of the outgoing particle, Zh = Eh/Nu
  std::vector<Float_t> ThetaPQ;   //
  std::vector<Float_t> Pt2;       //
  std::vector<Float_t> Pl2;       //
  std::vector<Float_t> PhiPQ;     //
  std::vector<Float_t> ThetaLab;  //
  std::vector<Float_t> PhiLab;    //
  std::vector<Float_t> Mx2;       //
  std::vector<Float_t> T;         //
  std::vector<Float_t> vxh;       // x-vertex of electron, in (cm)
  std::vector<Float_t> vyh;       // y-vertex of electron, in (cm)
  std::vector<Float_t> vzh;       // z-vertex of electron, in (cm)
  std::vector<Int_t> Sector;      //
  std::vector<Float_t> Px;        // x-component of particle's momentum, in (GeV)
  std::vector<Float_t> Py;        // y-component of particle's momentum, in (GeV)
  std::vector<Float_t> Pz;        // z-component of particle's momentum, in (GeV)
  std::vector<Float_t> P;         // particle's momentum, in (GeV)
  std::vector<Float_t> Betta;     //
  std::vector<Float_t> Mass2;     //
  std::vector<Float_t> deltaZ;    //
  std::vector<Int_t> pid;    
  Int_t TargTypeSM;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargType()
  Int_t TargTypeHH;      // 1 for liquid target, 2 for solid target, criteria written in TIdentificatorV2.h/TargType()

       //
};

#endif
