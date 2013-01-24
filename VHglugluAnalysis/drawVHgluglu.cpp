#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"





int main(int argc, char* argv[]) {

  if(  argc != 2 ) {
    std::cout << "USAGE: ./drawVHgluglu [(string)selType]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);



  DrawBase* db = new DrawBase("VHgluglu");

  db->set_lumiOnRightSide();

  std::string outputdir_str = "VHglugluPlots_MConly_" + selType;
  db->set_outputdir(outputdir_str);

  int signalFillColor = 46;

  std::string VHFileName = "VHgluglu_ZH_ZToLL_HToGluGlu_M_125_8TeV_powheg_herwigpp";
  VHFileName += "_" + selType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  db->add_mcFile( VHFile, "VH", "VH (125)", signalFillColor, 0);


  std::string mcDYFileName = "VHgluglu_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2";
  mcDYFileName += "_" + selType;
  mcDYFileName += ".root";
  TFile* mcDYFile = TFile::Open(mcDYFileName.c_str());
  db->add_mcFile( mcDYFile, "DY", "DY", 29, 3005);



  db->set_shapeNormalization();




  bool log = true;


  //db->set_legendTitle("Trilepton channel");

  db->set_yAxisMaxScale(1.1);

  db->drawHisto("nJets", "Jet Multiplicity", "", "Events", log);
  db->set_xAxisMax(6.5);
  db->drawHisto("nBJets_loose", "b-Jet Multiplicity (loose)", "", "Events", log);
  db->drawHisto("nBJets_medium", "b-Jet Multiplicity (medium)", "", "Events", log);

  db->set_xAxisMax();
  db->set_rebin(4);
  db->drawHisto("mZll", "Dilepton mass", "GeV", "Events");
  db->set_rebin(10);
  db->drawHisto("ptZll", "Dilepton p_{T}", "GeV", "Events");

  db->set_rebin(4);
  db->drawHisto("ptJet1", "Lead Jet p_{T}", "GeV", "Events");
  db->drawHisto("ptJet2", "Sublead Jet p_{T}", "GeV", "Events");

  db->set_rebin(2);
  db->drawHisto("qglJet1", "Lead Jet Quark/Gluon LD", "", "Events");
  db->drawHisto("qglJet2", "Sublead Jet Quark/Gluon LD", "", "Events");
  db->drawHisto("qglProd", "Quark/Gluon LD Product", "", "Events");

  db->drawHisto("cosThetaStar", "cos(#theta^{*})", "", "Events");


  db->set_lumiNormalization(20000.);
  db->set_noStack(false);
  db->drawHisto("mjj", "DiJet Mass", "GeV", "Events");

  TH1F::AddDirectory(kTRUE);

  // computeYields
  TTree* sigTree = (TTree*)VHFile->Get("tree_passedEvents");
  TTree* bgTree = (TTree*)mcDYFile->Get("tree_passedEvents");

  TH1D* h1_mjj_sig = new TH1D("mjj_sig", "", 100, 0., 250.);
  h1_mjj_sig->Sumw2();
  TH1D* h1_mjj_bg = new TH1D("mjj_bg", "", 100, 0., 250.);
  h1_mjj_bg->Sumw2();


  sigTree->Project( "mjj_sig", "mjj", "eventWeight*(mjj>105. && mjj<120.)" );
  bgTree->Project( "mjj_bg", "mjj",   "eventWeight*(mjj>105. && mjj<120.)" );

  float lumi = 20000.;

  float s = lumi*h1_mjj_sig->Integral();
  float b = lumi*h1_mjj_bg->Integral();

  float signal_xsec = 0.3943*8.57E-02*0.0337*3.;
  float total_signal = signal_xsec*lumi;
  float effS = s/total_signal;

  std::cout << "s: " << s << " (" << 100.*effS << ") bg: " << b << std::endl;

  float ul = CLA( lumi, 0., effS, 0., b, 0. );

  std::cout << "UL/SM: " << ul/signal_xsec << std::endl;

  return 0;

}  


