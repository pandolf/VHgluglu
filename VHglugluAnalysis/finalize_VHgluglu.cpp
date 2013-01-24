
#include "Ntp1Finalizer_VHgluglu.h"
#include "TMath.h"
#include <iostream>








int main( int argc, char* argv[] ) {

  if( argc!=3 ) {
    std::cout << "USAGE: ./finalize_VHgluglu [dataset] [selectionType]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);



  Ntp1Finalizer_VHgluglu* nf = new Ntp1Finalizer_VHgluglu( dataset, selectionType );


  if( dataset=="DATA_Run2011_FULL" ) {
   
    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v2");
    nf->addFile("MuEG_Run2011A_FULL");
    nf->addFile("MuEG_Run2011B");
    //nf->addFile("SingleMu_Run2011A_FULL");
    //nf->addFile("SingleMu_Run2011B_v2");

  } else if( dataset=="VV_Summer11" ) {

    nf->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");

  } else if( dataset=="DY_VV" ) {

    nf->addFile("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11");
    nf->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");

  } else if( dataset=="BG" ) {

    nf->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("TTJ_Fall11_highstat");
    nf->addFile("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11");
    nf->addFile("GVJets_7TeV-madgraph_Fall11");

  } else if( dataset=="DYToLL_M-20_CT10_TuneZ2_7TeV-powheg-pythia" ) {

    nf->addFile("DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia");
    nf->addFile("DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia");

  } else if( dataset=="DYJetsToLL_TuneZ2_7TeV-madgraph-tauola_Fall11" ) {

    nf->addFile("DYJetsToLL_M-10To50_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1");
    nf->addFile("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11");
  
  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}


