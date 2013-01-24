// ------------------------------------------------------------
//  
//    Ntp1Finalizer_VHgluglu - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/interface/BTagSFUtil.h"



class Ntp1Finalizer_VHgluglu : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_VHgluglu( const std::string& dataset, const std::string& selectionType );
  virtual ~Ntp1Finalizer_VHgluglu() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  float get_btagThresh( const std::string& btag_OP_ );

  void set_leptSyst( int leptSyst ) { leptSyst_ = leptSyst; };
  void set_jes( int jes ) { jes_ = jes; };
  void set_btagSyst( int btagSyst ) { btagSyst_ = btagSyst; };

  void set_jer( bool jer ) { jer_ = jer; };

 private:

   std::string selectionType_;

   float  ptJet_thresh_;
   float  etaJet_thresh_;

   int njets_thresh_;

   float  ptLeptZ1_thresh_;
   float  ptLeptZ2_thresh_;
   float  etaLeptZ1_thresh_;
   float  etaLeptZ2_thresh_;


   float  mZll_threshLo_;
   float  mZll_threshHi_;

   float  ptZll_thresh_;

   float  absCosThetaStar_thresh_;

   float  qglJet1_thresh_;
   float  qglJet2_thresh_;


   // syst flags:
   int leptSyst_; // number of sigmas
   int jes_; // number of sigmas
   int btagSyst_; // number of sigmas
   bool jer_;

};

