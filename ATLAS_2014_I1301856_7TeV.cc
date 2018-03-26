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
      _h_ttbar_7 = bookHisto1D("checkScatter2D_7TeV");
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
	_h_ttbar_7->fill(7.0, event.weight());  
    }

    

    /// Normalise histograms etc., after the run
    void finalize() {

      double BR = 0.032; // branching ratio
      double SF = crossSection()*picobarn/sumOfWeights()/BR;  // scale factor

      scale(_h_ttbar_7, SF);
    }
    
    //@}
   
    /// @name Histograms

    //@{
   
    Histo1DPtr _h_ttbar_7;

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1301856);


}
