#include "Ntp1Finalizer_VHgluglu.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TRandom3.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"





bool USE_MC_MASS=false;
bool choose_maxPtPair = false;

int DEBUG_EVENT = -1;


struct ValueAndError {

  float val;
  float err;

};


ValueAndError getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod );
ValueAndError getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod );
ValueAndError getEventHLTSF( ValueAndError effSingle1, ValueAndError effSingle2, ValueAndError effDouble1, ValueAndError effDouble2 );
ValueAndError getMuonRecoSF( float pt, float eta );
ValueAndError getMuonIsoSF( float pt, float eta );
ValueAndError getElectronRecoSF( float pt, float eta );
ValueAndError getElectronIsoSF( float pt, float eta );
float getJERSF( float eta );


// constructor:

Ntp1Finalizer_VHgluglu::Ntp1Finalizer_VHgluglu( const std::string& dataset, const std::string& selectionType ) : Ntp1Finalizer( "VHgluglu", dataset ) {


  leptSyst_ = 0;
  jes_ = 0;
  btagSyst_ = 0;
  jer_ = false;

  setSelectionType(selectionType);

}




void Ntp1Finalizer_VHgluglu::finalize( ) {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 160000);

  std::string fullFlags = selectionType_;
  if( jes_ == 1 ) fullFlags = fullFlags + "_JESUP";
  else if( jes_ == -1 ) fullFlags = fullFlags + "_JESDOWN";
  else if( jes_ != 0 ) {
    char fullFlags_char[100];
    if( jes_ < 0 )
      sprintf( fullFlags_char, "%s_JESDOWN%d", fullFlags.c_str(), jes_ );
    else
      sprintf( fullFlags_char, "%s_JESUP%d", fullFlags.c_str(), jes_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( btagSyst_ == 1 ) fullFlags = fullFlags + "_BTagUP";
  else if( btagSyst_ == -1 ) fullFlags = fullFlags + "_BTagDOWN";
  else if( btagSyst_ != 0 ) {
    char fullFlags_char[100];
    if( btagSyst_ < 0 )
      sprintf( fullFlags_char, "%s_BTagDOWN%d", fullFlags.c_str(), btagSyst_ );
    else
      sprintf( fullFlags_char, "%s_BTagUP%d", fullFlags.c_str(), btagSyst_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( leptSyst_ == 1 ) fullFlags = fullFlags + "_LeptUP";
  else if( leptSyst_ == -1 ) fullFlags = fullFlags + "_LeptDOWN";
  else if( leptSyst_ != 0 ) {
    char fullFlags_char[100];
    if( leptSyst_ < 0 )
      sprintf( fullFlags_char, "%s_LeptDOWN%d", fullFlags.c_str(), leptSyst_ );
    else
      sprintf( fullFlags_char, "%s_LeptUP%d", fullFlags.c_str(), leptSyst_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( jer_ ) fullFlags = fullFlags + "_JERUP";
  this->set_flags(fullFlags); //this is for the outfile name
  this->createOutputFile();



  TTree* tree_passedEvents = new TTree("tree_passedEvents", "Unbinned data for statistical treatment");


  TH1D* h1_nCounter = new TH1D("nCounter", "", 1, 0., 1.);
  h1_nCounter->Sumw2();
  TH1D* h1_nCounterW = new TH1D("nCounterW", "", 1, 0., 1.);
  h1_nCounterW->Sumw2();
  TH1D* h1_nCounterPU = new TH1D("nCounterPU", "", 1, 0., 1.);
  h1_nCounterPU->Sumw2();


  TH1D* h1_nvertex = new TH1D("nvertex", "", 36, -0.5, 35.5);
  h1_nvertex->Sumw2();
  TH1D* h1_nvertex_PUW = new TH1D("nvertex_PUW", "", 36, -0.5, 35.5);
  h1_nvertex_PUW->Sumw2();
  TH1D* h1_nvertex_PUW_ave = new TH1D("nvertex_PUW_ave", "", 36, -0.5, 35.5);
  h1_nvertex_PUW_ave->Sumw2();


  TH1D* h1_rhoPF_noPUW = new TH1D("rhoPF_noPUW", "", 50, 0., 30.);
  h1_rhoPF_noPUW->Sumw2();
  TH1D* h1_rhoPF_prepresel = new TH1D("rhoPF_prepresel", "", 50, 0., 30.);
  h1_rhoPF_prepresel->Sumw2();
  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 30.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 30.);
  h1_rhoPF->Sumw2();


  TH1D* h1_ptLeptZ1 = new TH1D("ptLeptZ1", "", 500, 20., 520.);
  h1_ptLeptZ1->Sumw2();
  TH1D* h1_ptLeptZ2 = new TH1D("ptLeptZ2", "", 200, 20., 220.);
  h1_ptLeptZ2->Sumw2();
  TH1D* h1_etaLeptZ1 = new TH1D("etaLeptZ1", "", 50, -2.5, 2.5);
  h1_etaLeptZ1->Sumw2();
  TH1D* h1_etaLeptZ2 = new TH1D("etaLeptZ2", "", 50, -2.5, 2.5);
  h1_etaLeptZ2->Sumw2();


  TH1D* h1_deltaRZll = new TH1D("deltaRZll", "", 500, 0., 5.);
  h1_deltaRZll->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 400., 0., 400.);
  h1_ptZll->Sumw2();
  TH1D* h1_etaZll = new TH1D("etaZll", "", 200, -5., 5.);
  h1_etaZll->Sumw2();
  TH1D* h1_mZll = new TH1D("mZll", "", 220, 50., 160.);
  h1_mZll->Sumw2();


  TH1D* h1_nJets = new TH1D("nJets", "", 11, -0.5, 10.5);
  h1_nJets->Sumw2();
  TH1D* h1_nBJets_loose = new TH1D("nBJets_loose", "", 11, -0.5, 10.5);
  h1_nBJets_loose->Sumw2();
  TH1D* h1_nBJets_medium = new TH1D("nBJets_medium", "", 11, -0.5, 10.5);
  h1_nBJets_medium->Sumw2();

  TH1D* h1_mjj = new TH1D("mjj", "", 100, 50., 250.);
  h1_mjj->Sumw2();
  TH1D* h1_ptjj = new TH1D("ptjj", "", 50, 0., 200.);
  h1_ptjj->Sumw2();
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 100, 0., 300.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 100, 0., 200.);
  h1_ptJet2->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 100, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 100, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_qglJet1 = new TH1D("qglJet1", "", 50, 0., 1.0001);
  h1_qglJet1->Sumw2();
  TH1D* h1_qglJet2 = new TH1D("qglJet2", "", 50, 0., 1.0001);
  h1_qglJet2->Sumw2();
  TH1D* h1_qglProd = new TH1D("qglProd", "", 50, 0., 1.0001);
  h1_qglProd->Sumw2();
  TH1D* h1_cosThetaStar = new TH1D("cosThetaStar", "", 50, -1., 1.0001);
  h1_cosThetaStar->Sumw2();






  Int_t nPU;
  tree_->SetBranchAddress("nPU", &nPU);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  unsigned int event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeightPU;
  tree_->SetBranchAddress("eventWeightPU", &eventWeightPU);
  Float_t eventWeightPU_ave;
  tree_->SetBranchAddress("eventWeightPU_ave", &eventWeightPU_ave);
  //Float_t eventWeight_Zee;
  //tree_->SetBranchAddress("eventWeight_Zee", &eventWeight_Zee);
  //Float_t eventWeight_Zmm;
  //tree_->SetBranchAddress("eventWeight_Zmm", &eventWeight_Zmm);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree_->SetBranchAddress("epfMet", &pfMet);
  Float_t metSignificance;
  tree_->SetBranchAddress("metSignificance", &metSignificance);
  Float_t mEtSig;
  tree_->SetBranchAddress("mEtSig", &mEtSig);
  Float_t phiMet;
  tree_->SetBranchAddress("phipfMet", &phiMet);


  int leptType;
  tree_->SetBranchAddress("leptType", &leptType);

  Float_t eLeptZ1;
  tree_->SetBranchAddress("eLeptZ1", &eLeptZ1);
  Float_t ptLeptZ1;
  tree_->SetBranchAddress("ptLeptZ1", &ptLeptZ1);
  Float_t etaLeptZ1;
  tree_->SetBranchAddress("etaLeptZ1", &etaLeptZ1);
  Float_t phiLeptZ1;
  tree_->SetBranchAddress("phiLeptZ1", &phiLeptZ1);
  Int_t chargeLeptZ1;
  tree_->SetBranchAddress("chargeLeptZ1", &chargeLeptZ1);
  Float_t combinedIsoRelLeptZ1;
  tree_->SetBranchAddress("combinedIsoRelLeptZ1", &combinedIsoRelLeptZ1);

  Float_t eLeptZ2;
  tree_->SetBranchAddress("eLeptZ2", &eLeptZ2);
  Float_t ptLeptZ2;
  tree_->SetBranchAddress("ptLeptZ2", &ptLeptZ2);
  Float_t etaLeptZ2;
  tree_->SetBranchAddress("etaLeptZ2", &etaLeptZ2);
  Float_t phiLeptZ2;
  tree_->SetBranchAddress("phiLeptZ2", &phiLeptZ2);
  Int_t chargeLeptZ2;
  tree_->SetBranchAddress("chargeLeptZ2", &chargeLeptZ2);
  Float_t combinedIsoRelLeptZ2;
  tree_->SetBranchAddress("combinedIsoRelLeptZ2", &combinedIsoRelLeptZ2);


  Int_t nLept;
  tree_->SetBranchAddress("nLept", &nLept);
  Int_t leptTypeLept[10];
  tree_->SetBranchAddress("leptTypeLept", leptTypeLept);
  Float_t eLept[10];
  tree_->SetBranchAddress("eLept", eLept);
  Float_t ptLept[10];
  tree_->SetBranchAddress("ptLept", ptLept);
  Float_t etaLept[10];
  tree_->SetBranchAddress("etaLept", etaLept);
  Float_t phiLept[10];
  tree_->SetBranchAddress("phiLept", phiLept);
  Int_t chargeLept[10];
  tree_->SetBranchAddress("chargeLept", chargeLept);
  Float_t combinedIsoRelLept[10];
  tree_->SetBranchAddress("combinedIsoRelLept", combinedIsoRelLept);


  Int_t nJets;
  tree_->SetBranchAddress("nJets", &nJets);

  Float_t eJet[50];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[50];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t ptUncertJet[50];
  tree_->SetBranchAddress("ptUncertJet", ptUncertJet);
  Float_t etaJet[50];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[50];
  tree_->SetBranchAddress("phiJet", phiJet);
  Float_t eGenJet[50];
  tree_->SetBranchAddress("eGenJet", eGenJet);
  Float_t ptGenJet[50];
  tree_->SetBranchAddress("ptGenJet", ptGenJet);
  Float_t etaGenJet[50];
  tree_->SetBranchAddress("etaGenJet", etaGenJet);
  Float_t phiGenJet[50];
  tree_->SetBranchAddress("phiGenJet", phiGenJet);
  Float_t eChargedHadronsJet[50];
  tree_->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t rmsCandJet[50];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[50];
  tree_->SetBranchAddress("ptDJet", ptDJet);
  Float_t ptD_QCJet[50];
  tree_->SetBranchAddress("ptD_QCJet", ptD_QCJet);
  Float_t axis2_QCJet[50];
  tree_->SetBranchAddress("axis2_QCJet", axis2_QCJet);
  Int_t nPFCand_QC_ptCutJet[50];
  tree_->SetBranchAddress("nPFCand_QC_ptCutJet", nPFCand_QC_ptCutJet);
  Int_t nChargedJet[50];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[50];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t eMuonsJet[50];
  tree_->SetBranchAddress("eMuonsJet", eMuonsJet);
  Float_t eElectronsJet[50];
  tree_->SetBranchAddress("eElectronsJet", eElectronsJet);
  Float_t trackCountingHighEffBJetTagJet[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet", trackCountingHighEffBJetTagJet);
  Float_t trackCountingHighPurBJetTagJet[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet", trackCountingHighPurBJetTagJet);
  Float_t simpleSecondaryVertexHighEffBJetTagJet[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet", simpleSecondaryVertexHighEffBJetTagJet);
  Float_t simpleSecondaryVertexHighPurBJetTagJet[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet", simpleSecondaryVertexHighPurBJetTagJet);
  Float_t combinedSecondaryVertexBJetTagJet[50];
  tree_->SetBranchAddress("combinedSecondaryVertexBJetTagJet", combinedSecondaryVertexBJetTagJet);
  Float_t jetBProbabilityBJetTagJet[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet", jetBProbabilityBJetTagJet);
  Float_t jetProbabilityBJetTagJet[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet", jetProbabilityBJetTagJet);
  Int_t pdgIdPartJet[50];
  tree_->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  Float_t etaPartJet[50];
  tree_->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[50];
  tree_->SetBranchAddress("phiPartJet", phiPartJet);



  //Int_t nPart;
  //tree_->SetBranchAddress("nPart", &nPart);
  //Float_t ePart[20];
  //tree_->SetBranchAddress("ePart", ePart);
  //Float_t ptPart[20];
  //tree_->SetBranchAddress("ptPart", ptPart);
  //Float_t etaPart[20];
  //tree_->SetBranchAddress("etaPart", etaPart);
  //Float_t phiPart[20];
  //tree_->SetBranchAddress("phiPart", phiPart);
  //Int_t pdgIdPart[20];
  //tree_->SetBranchAddress("pdgIdPart", pdgIdPart);


  // HLT:
  Bool_t passed_HLT_DoubleMu6;
  tree_->SetBranchAddress("passed_HLT_DoubleMu6", &passed_HLT_DoubleMu6);
  Bool_t passed_HLT_DoubleMu7;
  tree_->SetBranchAddress("passed_HLT_DoubleMu7", &passed_HLT_DoubleMu7);
  Bool_t passed_HLT_Mu13_Mu8;
  tree_->SetBranchAddress("passed_HLT_Mu13_Mu8", &passed_HLT_Mu13_Mu8);
  Bool_t passed_HLT_IsoMu17;
  tree_->SetBranchAddress("passed_HLT_IsoMu17", &passed_HLT_IsoMu17);
  Bool_t passed_HLT_IsoMu24;
  tree_->SetBranchAddress("passed_HLT_IsoMu24", &passed_HLT_IsoMu24);
  Bool_t passed_HLT_Mu8_Jet40;
  tree_->SetBranchAddress("passed_HLT_Mu8_Jet40", &passed_HLT_Mu8_Jet40);
  Bool_t passed_HLT_L2DoubleMu23_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu23_NoVertex", &passed_HLT_L2DoubleMu23_NoVertex);
  Bool_t passed_HLT_L2DoubleMu30_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu30_NoVertex", &passed_HLT_L2DoubleMu30_NoVertex);
  Bool_t passed_HLT_TripleMu5;
  tree_->SetBranchAddress("passed_HLT_TripleMu5", &passed_HLT_TripleMu5);

  Bool_t passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  Bool_t passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
  





  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  BTagSFUtil* btsfutil = new BTagSFUtil("CSV", 13);
  std::string meanminmax;
  if( btagSyst_==0 )       meanminmax = "mean";
  else if( btagSyst_==1 )  meanminmax = "max";
  else if( btagSyst_==-1 ) meanminmax = "min";
  else {
    std::cout << "Only allowed values for btagSyst are 0, 1, -1. Exiting." << std::endl;
    exit(55);
  }

  
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
 
  float Zmass = 91.1876;
  float tmass = 172.9;


  std::string puType = "Summer11_S4";
  TString dataset_tstr(dataset_);
  if( dataset_tstr.Contains("Fall11") ) {
    puType = "Fall11";
  }




  QGLikelihoodCalculator* qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();



  float mZll_t, ptZll_t;
  float ptLeptZ1_t, ptLeptZ2_t, etaLeptZ1_t, etaLeptZ2_t;
  float ptJetB1_t, ptJetB2_t, etaJetB1_t, etaJetB2_t;
  //float bTagJetB1_t, bTagJetB2_t;
  //float ptJet3_t, ptJet4_t, etaJet3_t, etaJet4_t;
  float HLTSF, leptonSF;
  bool isMZllSignalRegion;
  bool passed_btag;
  int njets;
  int nBjets_loose;
  int nBjets_medium;
  float ptJet1, ptJet2;
  float etaJet1, etaJet2;
  float qglJet1, qglJet2;
  float mjj;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "ptLeptZ1", &ptLeptZ1_t, "ptLeptZ1_t/F" );
  tree_passedEvents->Branch( "ptLeptZ2", &ptLeptZ2_t, "ptLeptZ2_t/F" );
  tree_passedEvents->Branch( "etaLeptZ1", &etaLeptZ1_t, "etaLeptZ1_t/F" );
  tree_passedEvents->Branch( "etaLeptZ2", &etaLeptZ2_t, "etaLeptZ2_t/F" );
  tree_passedEvents->Branch( "ptZll", &ptZll_t, "ptZll_t/F" );
  tree_passedEvents->Branch( "mZll", &mZll_t, "mZll_t/F" );
  tree_passedEvents->Branch( "njets", &njets, "njets/I" );
  tree_passedEvents->Branch( "nBjets_loose", &nBjets_loose, "nBjets_loose/I" );
  tree_passedEvents->Branch( "nBjets_medium", &nBjets_medium, "nBjets_medium/I" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1, "ptJet1/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2, "ptJet2/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1, "etaJet1/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2, "etaJet2/F" );
  tree_passedEvents->Branch( "qglJet1", &qglJet1, "qglJet1/F" );
  tree_passedEvents->Branch( "qglJet2", &qglJet2, "qglJet2/F" );
  tree_passedEvents->Branch( "mjj", &mjj, "mjj/F" );
      



ofstream ofs("run_event.txt");


  TRandom3* rand = new TRandom3(1313);


  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  if( jes_!= 0 )
    std::cout << "----> JES SYST: " << jes_ << " sigma" << std::endl;
  if( btagSyst_!= 0 )
    std::cout << "----> BTag SYST: " << btagSyst_ << " sigma" << std::endl;
  if( leptSyst_!= 0 )
    std::cout << "----> LEPT SYST: " << leptSyst_ << " sigma" << std::endl;
  if( jer_ )
    std::cout << "----> JER SYST ON" << std::endl;
  std::cout << std::endl << std::endl;



  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    // need this to be compatible with trigger!!! (just put it in analyzer, but dont want to rerun yet)
    if( ptLeptZ1<20. ) continue;
    if( ptLeptZ2<20. ) continue;


//std::cout << "new event" << std::endl;

    if( eventWeight <= 0. ) eventWeight = 1.;



    if( !isMC ) { 

      // remove duplicate events:

      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }


    } //if is not mc

   

    

    h1_nvertex->Fill(nvertex, eventWeight);
    h1_rhoPF_noPUW->Fill( rhoPF, eventWeight);


    // first leptons:
    
    TLorentzVector leptZ1, leptZ2;
    leptZ1.SetPtEtaPhiE( ptLeptZ1, etaLeptZ1, phiLeptZ1, eLeptZ1 );
    leptZ2.SetPtEtaPhiE( ptLeptZ2, etaLeptZ2, phiLeptZ2, eLeptZ2 );

    TLorentzVector diLepton = leptZ1+leptZ2;

    if( event==DEBUG_EVENT ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "** LOG FOR RUN: " << run << "   EVENT: " << DEBUG_EVENT << std::endl << std::endl;
      std::cout << "eventWeight: " << eventWeight << std::endl;
      std::cout << "leptType: " << leptType << std::endl; 
      std::cout << "leptZ1.Pt(): " << leptZ1.Pt() << " leptZ1.Eta(): " << leptZ1.Eta() << std::endl;
      std::cout << "leptZ2.Pt(): " << leptZ2.Pt() << " leptZ2.Eta(): " << leptZ2.Eta() << std::endl;
      std::cout << "diLepton.M(): " << diLepton.M() << std::endl;
    }


    if( diLepton.M()<50. ) continue; // gen cut in DY sample
    if( diLepton.M()<mZll_threshLo_ || diLepton.M()>mZll_threshHi_ ) continue; 
    if( diLepton.Pt()<ptZll_thresh_ ) continue; 



    // then jets:
    
    njets=0;
    nBjets_loose = 0;
    nBjets_medium  = 0;

  
    std::vector< TLorentzVector > selectedJets;
    std::vector< TLorentzVector > matchedJets;

    AnalysisJet jet1;
    AnalysisJet jet2;


    for( unsigned iJet=0; iJet<nJets; ++iJet) {

      
      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);


      if( thisJet.Pt() < ptJet_thresh_ ) continue;
      if( fabs(thisJet.Eta()) > etaJet_thresh_ ) continue;

      thisJet.rmsCand = rmsCandJet[iJet];
      thisJet.ptD = ptDJet[iJet];
      thisJet.nCharged = nChargedJet[iJet];
      thisJet.nNeutral = nNeutralJet[iJet];
      thisJet.eMuons = eMuonsJet[iJet]/thisJet.Energy();
      thisJet.eElectrons = eElectronsJet[iJet]/thisJet.Energy();

      thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet[iJet];
      thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet[iJet];
      thisJet.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet[iJet];
      thisJet.jetProbabilityBJetTag               = jetProbabilityBJetTagJet[iJet];
      thisJet.combinedSecondaryVertexBJetTag      = combinedSecondaryVertexBJetTagJet[iJet];

      thisJet.QGLikelihood = qglikeli->computeQGLikelihood2012( thisJet.Pt(), thisJet.Eta(), rhoPF, nPFCand_QC_ptCutJet[iJet], ptD_QCJet[iJet], axis2_QCJet[iJet] );

      thisJet.ptGen = ptGenJet[iJet];
      thisJet.etaGen = etaGenJet[iJet];
      thisJet.phiGen = phiGenJet[iJet];
      thisJet.eGen = eGenJet[iJet];


      float thisBtag;
      thisBtag = thisJet.combinedSecondaryVertexBJetTag;

      bool isBtagged_loose = thisBtag > 0.244;
      bool isBtagged_medium = thisBtag > 0.679;


      njets += 1;
      if( isBtagged_loose ) nBjets_loose += 1;
      if( isBtagged_medium ) nBjets_medium += 1;

      selectedJets.push_back(thisJet);

      if( thisJet.Pt() > jet1.Pt() ) {
        jet2 = jet1;
        jet1 = thisJet;
      } else if( thisJet.Pt() > jet2.Pt() ) {
        jet2 = thisJet;
      }
   


    } // for jets


    if( njets < 2 ) continue;
    if( njets < njets_thresh_ ) continue;

    // bjet veto:
    if( nBjets_medium > 0 ) continue;



    if( choose_maxPtPair ) {

      float maxPt = 0.;

      for( unsigned i=0; i<selectedJets.size(); ++i ) {
        for( unsigned j=i+1; j<selectedJets.size(); ++j ) {

          TLorentzVector dijet_t = selectedJets[i] + selectedJets[j];
          if( dijet_t.Pt() > maxPt ) {
            maxPt = dijet_t.Pt();
            jet1 = selectedJets[i];
            jet2 = selectedJets[j];
          }
        } // for j
      } // for i

    }


    if( jet1.QGLikelihood > qglJet1_thresh_ ) continue;
    if( jet2.QGLikelihood > qglJet2_thresh_ ) continue;


    TLorentzVector diJet = jet1+jet2;

    


    HelicityLikelihoodDiscriminant::HelicityAngles hangles;
    if( chargeLeptZ1<0 ) hangles = LD->computeHelicityAngles(leptZ1, leptZ2, jet1, jet2);
    else                 hangles = LD->computeHelicityAngles(leptZ2, leptZ1, jet1, jet2);
    
    float cosThetaStar = hangles.helCosThetaStar;

    if( fabs(cosThetaStar) > absCosThetaStar_thresh_ ) continue;



    // fill some histos before requiring third lepton:
    h1_nvertex_PUW->Fill(nvertex, eventWeight);

    h1_rhoPF_prepresel->Fill( rhoPF, eventWeight);

    h1_nJets->Fill( njets, eventWeight );
    h1_nBJets_loose->Fill( nBjets_loose, eventWeight );
    h1_nBJets_medium->Fill( nBjets_medium, eventWeight );

    h1_mZll->Fill( diLepton.M(), eventWeight );
    h1_ptZll->Fill( diLepton.Pt(), eventWeight );

    h1_mjj->Fill( diJet.M(), eventWeight );
    h1_ptjj->Fill( diJet.Pt(), eventWeight );

    h1_ptJet1->Fill( jet1.Pt(), eventWeight );
    h1_ptJet2->Fill( jet2.Pt(), eventWeight );

    h1_etaJet1->Fill( jet1.Eta(), eventWeight );
    h1_etaJet2->Fill( jet2.Eta(), eventWeight );

    h1_qglJet1->Fill( jet1.QGLikelihood, eventWeight );
    h1_qglJet2->Fill( jet2.QGLikelihood, eventWeight );
    h1_qglProd->Fill( jet1.QGLikelihood*jet2.QGLikelihood, eventWeight );

    h1_cosThetaStar->Fill( cosThetaStar, eventWeight );


    ptLeptZ1_t = leptZ1.Pt();
    ptLeptZ2_t = leptZ2.Pt();
    etaLeptZ1_t = leptZ1.Eta();
    etaLeptZ2_t = leptZ2.Eta();

    ptZll_t = diLepton.Pt();
    mZll_t = diLepton.M();


    ptJet1 = jet1.Pt();
    ptJet2 = jet2.Pt();

    etaJet1 = jet1.Eta();
    etaJet2 = jet2.Eta();

    qglJet1 = jet1.QGLikelihood;
    qglJet2 = jet2.QGLikelihood;

    mjj = diJet.M();



    // and fill tree (remember this includes the mZll sidebands):
    tree_passedEvents->Fill();

  

  
  } //for entries



  h1_nCounter->SetBinContent(1, nCounter_);
  h1_nCounterW->SetBinContent(1, nCounterW_);
  h1_nCounterPU->SetBinContent(1, nCounterPU_);




  // write all stuff in files:

  outFile_->cd();

  tree_passedEvents->Write();

  h1_nCounter->Write();
  h1_nCounterW->Write();
  h1_nCounterPU->Write();

  h1_nvertex->Write();
  h1_nvertex_PUW->Write();
  h1_nvertex_PUW_ave->Write();


  h1_rhoPF_noPUW->Write();
  h1_rhoPF_prepresel->Write();
  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();


  h1_ptLeptZ1->Write();
  h1_ptLeptZ2->Write();
  h1_etaLeptZ1->Write();
  h1_etaLeptZ2->Write();

  h1_ptZll->Write();
  h1_mZll->Write();


  h1_nJets->Write();
  h1_nBJets_loose->Write();
  h1_nBJets_medium->Write();

  h1_ptJet1->Write();
  h1_ptJet2->Write();

  h1_etaJet1->Write();
  h1_etaJet2->Write();

  h1_mjj->Write();
  h1_ptjj->Write();

  h1_qglJet1->Write();
  h1_qglJet2->Write();
  h1_qglProd->Write();

  h1_cosThetaStar->Write();



  outFile_->Close();


} // finalize()



void Ntp1Finalizer_VHgluglu::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  ptLeptZ1_thresh_ = 20.;
  ptLeptZ2_thresh_ = 20.;
  etaLeptZ1_thresh_ = 3.;
  etaLeptZ2_thresh_ = 3.;

  ptJet_thresh_ = 20.;
  etaJet_thresh_ = 5;

  njets_thresh_ = 2;

  mZll_threshLo_ = 0.;
  mZll_threshHi_ = 10000.;

  ptZll_thresh_ = 0.;

  absCosThetaStar_thresh_ = 1.;

  qglJet1_thresh_ = 1.;
  qglJet2_thresh_ = 1.;


  if( selectionType_=="presel" ) {


  } else if( selectionType_=="sel0" ) {

    etaJet_thresh_ = 2.5;

    njets_thresh_ = 2;

    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;

    ptZll_thresh_ = 70.;

    absCosThetaStar_thresh_ = 0.6;


  } else if( selectionType_=="sel1" ) {

    etaJet_thresh_ = 2.5;

    njets_thresh_ = 2;

    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;

    ptZll_thresh_ = 70.;

    absCosThetaStar_thresh_ = 0.6;

    qglJet1_thresh_ = 0.1;
    qglJet2_thresh_ = 0.075;


  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType








ValueAndError getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod ) {

  float hltsf = 0.;
  float hltsf_err = 0.;

  // these numbers taken from AN2011-399-v4
  if( runPeriod=="Run2011A" ) {

    if( fabs(eta)<0.8 ) {
      hltsf = 0.975;
      hltsf_err = 0.004;
    } else if( fabs(eta)<2.1 ) {
      hltsf = 0.950;
      hltsf_err = 0.005;
    } else  {
      hltsf = 0.910;
      hltsf_err = 0.01;
    }

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) {
      if( pt<40. ) {
        hltsf = 0.977;
        hltsf_err = 0.001;
      } else {
        hltsf = 0.975;
        hltsf_err = 0.001;
      }
    } else if( fabs(eta)<2.1 ){
      if( pt<40. ) {
        hltsf = 0.955;
        hltsf_err = 0.002;
      } else {
        hltsf = 0.955;
        hltsf_err = 0.001;
      }
    } else { // eta 2.1 -> 2.4
      if( pt<40. ) {
        hltsf = 0.89;
        hltsf_err = 0.007;
      } else {
        hltsf = 0.90;
        hltsf_err = 0.006;
      }
    }

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  ValueAndError ve_hlt;
  ve_hlt.val = hltsf;
  ve_hlt.err = hltsf_err;

  return ve_hlt;

}




ValueAndError getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod ) {

  if( pt<25. ) { 
    ValueAndError ve_hlt;
    ve_hlt.val = 0.;
    ve_hlt.err = 0.;
    return ve_hlt;
  }

  float hltsf = 0.;
  float hltsf_err = 0.;

  if( runPeriod=="Run2011A1" ) { //up to may10 technical stop

    if( fabs(eta)<0.8 ) {
      hltsf = 0.896;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.807;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.608;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011A2" ) { //from may10 to EPS

    if( fabs(eta)<0.8 ) {
      hltsf = 0.895;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.838;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.738;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011A3" ) { //from EPS to end of Run2011A

    if( fabs(eta)<0.8 ) {
      hltsf = 0.890;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.809;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.493;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) {
      hltsf = 0.87;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.79;
      hltsf_err = 0.001;
    } else  { //using HLT_IsoMu24_eta2p1
      hltsf = 0.;
      hltsf_err = 0.;
    }

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  ValueAndError ve_hlt;
  ve_hlt.val = hltsf;
  ve_hlt.err = hltsf_err;

  return ve_hlt;


}


ValueAndError getEventHLTSF( ValueAndError effSingle1, ValueAndError effSingle2, ValueAndError effDouble1, ValueAndError effDouble2 ) {

  float HLTSF = effDouble1.val * effDouble2.val +
                effSingle2.val * (1. - effDouble2.val ) +
                effSingle1.val * (1. - effDouble1.val );

  float HLTSF_err = effDouble1.err * effDouble1.err* effDouble2.val * effDouble2.val +
                    effDouble1.val * effDouble1.val* effDouble2.err * effDouble2.err +
                    effSingle2.err * effSingle2.err * (1. - effDouble2.val ) * (1. - effDouble2.val ) +
                    effSingle2.val * effSingle2.val * effDouble2.err * effDouble2.err +
                    effSingle1.err * effSingle1.err * (1. - effDouble1.val ) * (1. - effDouble1.val ) +
                    effSingle1.val * effSingle1.val * effDouble1.err * effDouble1.err;

  HLTSF_err = sqrt( HLTSF_err );

  ValueAndError ve_hlt;
  ve_hlt.val = HLTSF;
  ve_hlt.err = HLTSF_err;

  return ve_hlt;

}



ValueAndError getMuonRecoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( fabs(eta)<1.2 ) {
    recoSF = 0.996;
    recoSF_err = 0.001;
  } else {
    recoSF = 0.986;
    recoSF_err = 0.001;
  }

  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}


ValueAndError getElectronRecoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( fabs(eta)<0.8 ) {
    recoSF = 0.999;
    recoSF_err = 0.005;
  } else if( fabs(eta)<1.44 ) {
    recoSF = 0.964;
    recoSF_err = 0.003;
  } else if( fabs(eta)<1.57 ) {
    recoSF = 0.99;
    recoSF_err = 0.04;
  } else if( fabs(eta)<2.0 ) {
    recoSF = 0.992;
    recoSF_err = 0.006;
  } else {
    recoSF = 1.001;
    recoSF_err = 0.006;
  }

  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}



ValueAndError getMuonIsoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( pt<40. ) {

    if( fabs(eta)<0.9 ) {
      recoSF = 0.987;
      recoSF_err = 0.006;
    } else {
      recoSF = 0.995;
      recoSF_err = 0.005;
    }
   
  } else {

    if( fabs(eta)<0.9 ) {
      recoSF = 0.994;
      recoSF_err = 0.002;
    } else {
      recoSF = 0.996;
      recoSF_err = 0.002;
    }
   
  }


  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}


ValueAndError getElectronIsoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( pt<40. ) {

    if( fabs(eta)<1.5 ) {
      recoSF = 0.988;
      recoSF_err = 0.006;
    } else {
      recoSF = 0.998;
      recoSF_err = 0.011;
    }
   
  } else {

    if( fabs(eta)<1.5 ) {
      recoSF = 0.988;
      recoSF_err = 0.003;
    } else {
      recoSF = 1.016;
      recoSF_err = 0.064;
    }
   
  }


  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}





float getJERSF( float eta ) {

  // scale factors taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
 
  float SF=-1.;

  if( fabs(eta)<0.5 ) {
    SF = 1.052;
  } else if( fabs(eta)<1.1 ) {
    SF = 1.057;
  } else if( fabs(eta)<1.7 ) {
    SF = 1.096;
  } else if( fabs(eta)<2.3 ) {
    SF = 1.134;
  } else {
    SF = 1.288;
  }

  if( SF<0. ) {
    std::cout << "JERSF is negative, this shouldn't be possible." << std::endl;
    exit(1123);
  }

  return SF;

}


