#ifndef MATH_HXX
#define MATH_HXX

/*** MATHEMATICAL FUNCTIONS ***/

Double_t PhiPQ(Double_t fPex, Double_t fPey, Double_t fPez, Double_t fPx, Double_t fPy, Double_t fPz) {
  Double_t fPhiPQ;
  // two 3-momentum vectors are defined
  // one for the hadron
  // one for the virtual photon
  TVector3 hadr(fPx, fPy, fPz);
  TVector3 virt(-fPex, -fPey, kEbeam - fPez);
  // the function Phi() and Theta() converts the vector into spherical coordinates
  // and the function RotateZ(amount) rotates the vector around z-axis a certain amount
  // now, this paragraph turns y-component of the virtual photon vector into 0 and rotates the hadron vector the same amount
  // to have the hadron vector in function of the virtual photon direction
  // if its not clear, try to draw these two vectors in the xy plane
  Double_t Phi_z = TMath::Pi() - virt.Phi();
  virt.RotateZ(Phi_z);
  hadr.RotateZ(Phi_z);
  // the function Angle() returns the angle between two vectors
  // now, take into account that the virtual photon vector is in the y=0 plane
  // again, this paragraph turns the x-component of the virtual photon vector into 0 and rotates the hadron vector the same amount
  // to have the hadron vector in function of the virtual photon direction
  // if its not clear, try to draw these two vectors in the xz plane
  TVector3 unit(0., 0., 1.);
  Double_t Phi_y = virt.Angle(unit);
  virt.RotateY(Phi_y);
  hadr.RotateY(Phi_y);
  // finally, it obtains the phi component (in spherical coordinates) of the hadron vector (now in the virtual photon frame of reference)
  fPhiPQ = hadr.Phi() * TMath::RadToDeg();

  return fPhiPQ;
}

Double_t ThetaPQ(Double_t fPex, Double_t fPey, Double_t fPez, Double_t fPx, Double_t fPy, Double_t fPz) {
  // angle between virtual photon and particle
  // (and this makes perfect sense if one draws it)
  Double_t fThetaPQ;
  TVector3 hadr(fPx, fPy, fPz);
  TVector3 virt(-fPex, -fPey, kEbeam - fPez);
  fThetaPQ = virt.Angle(hadr) * TMath::RadToDeg();

  return fThetaPQ;
}

Double_t ThetaLab(Double_t fPx, Double_t fPy, Double_t fPz) {
  // polar angle in lab frame
  TVector3 v3p(fPx, fPy, fPz);
  return v3p.Theta() * TMath::RadToDeg();
}

Double_t PhiLab(Double_t fPx, Double_t fPy, Double_t fPz) {
  // azimuthal angle in lab frame, range: (-30, 330)
  TVector3 v3p(fPx, fPy, fPz);
  Double_t PhiLab_value = v3p.Phi() * TMath::RadToDeg();
  if (PhiLab_value < -30.) {
    return PhiLab_value + 360.;
  } else if (PhiLab_value > 330.) {
    return PhiLab_value - 360.;
  }                     // closure
  return PhiLab_value;  // default
}

#endif
