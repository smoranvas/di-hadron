
// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"



namespace Rivet {
  /// @brief Correlations between a leading pi+ and a subleading pi- in a nuclear DIS reaction
  class correlations : public Analysis {
  private:
    double deg_to_eta(double theta){
      return -log(tan(theta*M_PI/180/2.));
    }
  public:
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(correlations);
    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections. Note that the definition
      // of the scattered lepton can be influenced by sepcifying
      // options as declared in the .info file.
      DISLepton lepton(options());
      declare(lepton, "Lepton");
      declare(DISKinematics(lepton), "Kinematics");
      //declare(FinalState(), "FS");
      
      //find the leading pi+ and subleading pi-s in the event
      //const DISFinalState& disfs = declare(DISFinalState(DISFinalState::BoostFrame::LAB), "DISFS");

      IdentifiedFinalState pips_ifs(211);
      IdentifiedFinalState pims_ifs(-211);
      
      // TODO write cuts for the leading and subleading pions
      /*auto momentum= sqrt(Cut::E*Cut::E-Cut::mass*Cut::mass);

      //auto theta=180/M_PI*2*arctan(exp(-Cut::Quantity::eta));
      
      auto pip_cuts= Cuts::etaIn(deg_to_eta(120),deg_to_eta(10));
      auto pim_cuts= Cuts::etaIn(-999,deg_to_eta(20)) && momentum>0.5;
      
      FinalState pips_cut(pips_ifs, pip_cuts);
      FinalState pims_cut(pims_ifs, pim_cuts);

      DISFinalState pips_dis(DISFinalState::BoostFrame::LAB, pips_cut);
      DISFinalState pims_dis(DISFinalState::BoostFrame::LAB, pims_cut);*/

      declare(pips_ifs, "pi+");
      declare(pims_ifs, "pi-");
      

      // Book histograms

      int nbins_dphi=8;
      vector< double > dphibins;
      for(int i = 0; i<nbins_dphi+1; i++)
	dphibins.push_back((M_PI*i)/nbins_dphi);
      
      book(_hist_dphi, "dphi", dphibins);
      _n_leading=0.0;

      
      vector<double> dybins={-0.5, 0.5, 1.5,2.5};
      vector<double> pt1bins={0.25, 0.4, 0.6, 1.0};
      vector<double> pt2bins={0.25, 0.4, 0.6, 0.8};
      
      book(_hist_pt1, "pt1", pt1bins);
      book(_hist_dphi_dy, "dphi_dy", dphibins, dybins);
      book(_hist_dphi_pt1, "dphi_pt1", dphibins, pt1bins);
      book(_hist_dphi_pt2, "dphi_pt2", dphibins, pt2bins);
      
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
      
      // Weight of the event
      //_hist_Q2->fill(Q2);
      //_hist_y->fill(y);
      //_hist_x->fill(x);

      // Momentum of the scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      if ( dl.failed() ) return;
      //cout << "found dis electron"<<endl;
      //const FourMomentum leptonMom = dl.out();
      //const double ptel = leptonMom.pT();
      //const double enel = leptonMom.E();

      //if(enel<11*GeV) vetoEvent;
      
      //const double thel = leptonMom.angle(dk.beamHadron().mom())/degree;
      // _hist_ept->fill(ptel);


       // Extract the leading pi+
      auto pips = apply<IdentifiedFinalState>(event, "pi+").particles();
      if(pips.size()==0)
	vetoEvent;
      //convert deg to rad
      double deg=M_PI/180;
      // Extract subleading pi-s
      auto pims = apply<IdentifiedFinalState>(event, "pi-").particles();
      for(auto pip: pips) {
	//cout << "pip"<<endl;
	double z1=pip.E()/nu;
	//do it this way so that it is independent of how the electron direction is defined (+z or -z).  
	double theta = pip.momentum().angle(dk.beamLepton().momentum());
	if (z1<0.5 || theta<10*deg)
	  continue;
	//cout << "pip passes cuts"<<endl;
	pip=pip.transformBy(dk.boostBreit());
	_n_leading+=1;
	double pt1=pip.pt();
	_hist_pt1->fill(pt1);
	for(auto pim:pims) {
	  double z2=pim.E()/nu;
	  double theta=pim.momentum().angle(dk.beamLepton().momentum());
	  if(z2<0.05 or z2>0.45 or
	     not ((theta>25*deg and pim.p()>0.7) or (pim.p()>0.5 and theta>30*deg) or (theta>40*deg and pim.p()>0.35 )))
	    continue;
	  //cout << "found pim that passes cuts" << endl;
	  pim=pim.transformBy(dk.boostBreit());
	  double dphi=pip.phi()-pim.phi();
	  if (dphi>M_PI)
	    dphi-=2*M_PI;
	  if (dphi<-M_PI)
            dphi+=2*M_PI;
	  dphi=abs(dphi);
	  _hist_dphi->fill(dphi);
	  double dy=pip.rap()-pim.rap();
	  double pt2=pim.pt();
	  _hist_dphi_dy->fill(dphi,dy);
	  _hist_dphi_pt1->fill(dphi,pt1);
          _hist_dphi_pt2->fill(dphi,pt2);
	    
	}
	
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //_hist_dphi->Divide(_n_leading);
      //_hist_dphi->Divide(
      //scale(_hist_Q2, 1.0/(_inclusive_xs)); // normalize to unity
      ///normalize(_hist_y); // normalize to unity
      //normalize(_hist_ept);
      //scale(_hist_jetpt, 1.0/(_inclusive_xs));
      //normalize(_hist_jeteta);
      //scale(_hist_qt, 1.0/(_inclusive_xs));
    }

    //@}


    /// The histograms.
    Histo1DPtr _hist_dphi;
    Histo2DPtr _hist_dphi_dy;
    Histo2DPtr _hist_dphi_pt1;
    Histo2DPtr _hist_dphi_pt2;
    double _n_leading;
    //used for normalizing the pt1-sliced histogram
    Histo1DPtr _hist_pt1;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(correlations);


}
