// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// @brief Parton-level top-pair cross-sections at 7TeV and 8TeV
  class ATLAS_2014_I1301856 : public Analysis {
  public:

    // Default Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1301856);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Parton level tops which decay to electrons and muons
      declare(PartonicTops(PartonicTops::ELECTRON), "ElectronPartonTops") ;
      declare(PartonicTops(PartonicTops::MUON), "MuonPartonTops") ;


      // Book histograms
      //       _h_ttbar = bookHisto1D("checkHisto1D", "sigma_ttbar_hist01d");
      _h_ttbar_both = bookHisto1D("checkScatter2D_both");
      _h_ttbar_7 = bookHisto1D("checkScatter2D_7TeV");
      _h_ttbar_8 = bookHisto1D("checkScatter2D_8TeV");
      //      _c_ttbar = bookCounter("checkScatter2d");
      _s_ttbar = bookScatter1D("checkScatter1D" , "sigma_ttbar_scatter1d");   
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const Particles electronpartontops = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();
      const Particles muonpartontops     = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();


      // veto if not an e mu pair
      const bool eMuPair = ( electronpartontops.size() == 1 && muonpartontops.size() == 1 );      
      if ( !eMuPair ) vetoEvent;

      // veto if e and mu are not of opposite chare
      const bool oppositeCharge ( electronpartontops[0].charge() + muonpartontops[0].charge() == 0 );
      if ( !oppositeCharge ) vetoEvent;


      // Fill histos
	_h_ttbar_both->fill(8.0, event.weight());  
	_h_ttbar_7->fill(7.0, event.weight());  
	_h_ttbar_8->fill(8.0, event.weight());  
	//_c_ttbar->fill(event.weight());
    }

    

    /// Normalise histograms etc., after the run
    void finalize() {
      double scatterval =0.0;
      double BR = 0.032; // branching ratio
      double SF = crossSection()*picobarn/sumOfWeights()/BR;  // scale factor
      
      //      cout << "_c_->numEntries=" << _c_ttbar->numEntries() << endl;
       scatterval = SF*( _h_ttbar_both->numEntries() );  
      // cout<< "scattelrval=" << scatterval;
      //cout << "eventCtr=" << eventCtr << endl;
      //double selectionRatio1 = selectionCtr/sumOfWeights();
      //double selectionRatio2 = selectionCtr/eventCtr;
      //double selectionRatio3 = selectionCtr/ttbarCtr;
      //      cout << "percetage of selected events=" << selectionRatio*100 << endl;
      //double ttbarXsec = selectionRatio1*crossSection()*picobarn/BR;
      //      double ttbarXsec = selectionCtr*crossSection()*picobarn/sumOfWeights()/BR;

      //cout << "sumOfWeights()="   << sumOfWeights()  << endl;
      //cout << "ttbarCtr="         << ttbarCtr        << endl;
      //cout << "selection ratio1=" << selectionRatio1 << endl;
      //cout << "selection ratio2=" << selectionRatio2 << endl;
      //cout << "selection ratio3=" << selectionRatio3 << endl;
      //cout << "crossSection()="   << crossSection()  << endl;  
      //cout << "ttbarXsec="        << ttbarXsec       << endl;
      
      scale(_h_ttbar_both, SF);
      scale(_h_ttbar_7, SF);
      scale(_h_ttbar_8, SF);
      // scale(_c_ttbar, SF);
      
       _s_ttbar->point(1).setX(scatterval, 0.0);
      // _s_ttbar->point(2).setY(140, 0.0); //ARBITRARY VALUE JUST NOW
      
    }
    
    //@}
    //    double eventCtr;
    //double selectionCtr;
    //double ttbarCtr;

    /// @name Histograms

    //@{
    //    Profile1DPtr _h_XXXX;
    //    Histo1DPtr _h_YYYY;
    CounterPtr _c_ttbar;
    Histo1DPtr _h_ttbar_both, _h_ttbar_7, _h_ttbar_8;

    Scatter1DPtr _s_ttbar;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1301856);


}
