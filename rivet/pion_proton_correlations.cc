// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class pion_proton_correlations : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(pion_proton_correlations);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections. Note that the definition
      // of the scattered lepton can be influenced by sepcifying
      // options as declared in the .info file.
      std::cout << options() << std::endl;
      DISLepton lepton;
      declare(lepton, "Lepton");
      declare(DISKinematics(lepton), "Kinematics");

      
      //find the leading pi+ and subleading pi-s in the event

      IdentifiedFinalState pips_ifs(211);
      IdentifiedFinalState protons_ifs(2212);
      
      declare(pips_ifs, "pi+");
      declare(protons_ifs, "p");
      

      // Book histograms

      int nbins_dphi=8;
      vector< double > dphibins;
      for(int i = 0; i<nbins_dphi+1; i++)
	dphibins.push_back((M_PI*i)/nbins_dphi);

      book(_n_leading, "n_leading");
      book(_hist_dphi, "dphi", dphibins);
      
      vector<double> dybins={0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      book(_hist_dphi_dy, "dphi_dy", dphibins, dybins);

      vector<double> dystarbins={-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
      book(_hist_dphi_dystar, "dphi_dystar", dphibins, dystarbins);

      vector<double> dycmbins={-1.0,-0.75, -0.5,-0.25, 0.0, 0.25, 0.5,0.75, 1.0,1.25, 1.5,1.75, 2.0};
      book(_hist_dycm, "dycm", dycmbins);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) return;
      double x  = dk.x();
      double Ebeam=dk.beamLepton().E();
      //      cout <<"Ebeam=" << Ebeam <<endl;
      double nu  = dk.y()*Ebeam;
      double Q2 = dk.Q2();
      double W = sqrt(dk.W2());
      
      if (Q2 < 1*GeV2) vetoEvent;
      if (W< 2) vetoEvent;
      // use same y range as in the data
      if (nu<2.3 or nu>4.2) vetoEvent;
      

      // Momentum of the scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      if ( dl.failed() ) return;
      
      auto CM = (dk.beamLepton().momentum()+dk.beamHadron().momentum()-dk.scatteredLepton().momentum());
      double dycm= 0.5*log((CM.E()+CM.p())/(CM.E()-CM.p()));
   
      
       // Extract the leading pi+
      auto pips = apply<IdentifiedFinalState>(event, "pi+").particles();
      if(pips.size()==0)
	vetoEvent;
      //convert deg to rad
      double deg=M_PI/180;
      // Extract subleading pi-s
      auto protons = apply<IdentifiedFinalState>(event, "p").particles();
      double ptmin=0.070;
      double pmin=0.350;
      for(auto pip: pips) {
	//cout << "pip"<<endl;
	double z1=pip.E()/nu;
	//do it this way so that it is independent of how the electron direction is defined (+z or -z).  
	double theta = pip.momentum().angle(dk.beamLepton().momentum());
	double p1=pip.p();
	pip=pip.transformBy(dk.boostHCM());
	double pt1=pip.pt();
	if (z1<0.5 || theta<10*deg || pt1<ptmin)
	  continue;
	//cout << "pip passes cuts"<<endl;
	_n_leading->fill();
	for(auto proton:protons) {
	  double p=proton.p();
	  double theta=proton.momentum().angle(dk.beamLepton().momentum());
	  proton=proton.transformBy(dk.boostHCM());
	  double pt2=proton.pt();
	  if(p<pmin or p>2.8 or theta<10*deg or pt2<ptmin)
	    continue;
	  //cout << "found proton that passes cuts" << endl;
	  double dphi=pip.phi()-proton.phi();
	  if (dphi>M_PI)
	    dphi-=2*M_PI;
	  if (dphi<-M_PI)
            dphi+=2*M_PI;
	  dphi=abs(dphi);

	  _hist_dphi->fill(dphi);
	  double dy=pip.rap()-proton.rap();
	  //std::cout << dy << "  " <<  dphi << std::endl;
	  _hist_dphi_dy->fill(dphi,dy);
	  
	  //include this to compare with earlier version
	  double dystar=dy-dycm;
	  _hist_dphi_dystar->fill(dphi,dystar);
	  _hist_dycm->fill(dycm);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    
    }

    //@}


    /// The histograms.

    Histo2DPtr _hist_dphi_dy;
    CounterPtr _n_leading;

    //for debugging:
    Histo1DPtr _hist_dphi;
    Histo2DPtr _hist_dphi_dystar;
    Histo1DPtr _hist_dycm;

  };


  RIVET_DECLARE_PLUGIN(pion_proton_correlations);

}
