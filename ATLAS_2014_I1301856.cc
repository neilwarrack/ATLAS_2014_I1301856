// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// @brief Parton-level top-pair cross-sections at 7TeV and 8TeV
  class ATLAS_2014_I1301856 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1301856);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(PartonicTops(PartonicTops::ELECTRON), "ElectronPartonTops") ;
      declare(PartonicTops(PartonicTops::MUON), "MuonPartonTops") ;

      // Book histograms
      //      _h_XXXX = bookProfile1D(1, 1, 1);
      //      _h_YYYY = bookHisto1D(2, 1, 1);
      _h_ZZZZ = bookCounter(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles electronpartontops = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt() ;
      const Particles muonpartontops = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt() ;
      const bool EMuPair = ( electronpartontops.size() == 1 && muonpartontops.size() == 1 ) ;

      if ( !EMuPair ) { 
	std::cout<< "1" << std::endl ; // VETO DUE TO: not being a electron-muon pair dileptonic decay 
	vetoEvent ;
      }


      _h_ZZZZ->fill(event.weight()) ;

      std::cout << "9" << std::endl ; // SUCCESS!
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //      normalize(_h_YYYY); // normalize to unity
      scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}


    /// @name Histograms
    //@{
    //    Profile1DPtr _h_XXXX;
    //    Histo1DPtr _h_YYYY;
    CounterPtr _h_ZZZZ;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1301856);


}
