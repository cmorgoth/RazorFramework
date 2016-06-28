//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TROOT.h>
//LOCAL INCLUDES
#include "RunTwoFitMgg.hh"
#include "CommandLineInput.hh"
#include "PrintAIC.hh"

int main( int argc, char* argv[])
{
  gROOT->Reset();
  std::string inputFile = ParseCommandLine( argc, argv, "-inputFile=" );
  if (  inputFile == "" )
    {
      std::cerr << "[ERROR]: please provide an input file using --inputFile=<path_to_file>" << std::endl;
      return -1;
    }

  std::string inputFile2 = ParseCommandLine( argc, argv, "-inputFile2=" );
  if (  inputFile2 == "" )
    {
      std::cerr << "[INFO]: using a single root file, provide an second input file using --inputFile2=<path_to_file>" << std::endl;
    }
  
  std::string treeName = ParseCommandLine( argc, argv, "-treeName=" );
  if (  treeName == "" )
    {
      std::cout << "[WARNING]: user did not provide a tree name, if you wish to provide one use --treeName=<yourTreeName>" << std::endl;
      treeName = "HighRes";
    }

  std::string LowMRcut = ParseCommandLine( argc, argv, "-LowMRcut=" );
  std::string HighMRcut = ParseCommandLine( argc, argv, "-HighMRcut=" );
  std::string LowRSQcut = ParseCommandLine( argc, argv, "-LowRSQcut=" );
  std::string HighRSQcut = ParseCommandLine( argc, argv, "-HighRSQcut=" );
  
  std::string BinSelection = ParseCommandLine( argc, argv, "-BinSelection=" );
  if (  BinSelection == "" && ( LowMRcut == "" || LowRSQcut == "" ) )
    {
      std::cout << "[ERROR]: user must provide either a BinSelection or an MR cut and Rsq cut. use --BinSelection=<yourBin> or --LowMRcut=<yourLowMRCut> --LowRSQcut=<yourLowRSQCut>" << std::endl;
      return -1;
    }

  if (  BinSelection == "" && ( HighMRcut == "" || HighRSQcut == "" ) )
    {
      std::cout << "[WARNING]: user does not provide either a High-MR cut or High-Rsq cut. use --HighMRcut=<yourHighMRCut> --HighRSQcut=<HighRSQCut>" << std::endl;
    }


  std::string dataMode = ParseCommandLine( argc, argv, "-dataMode=" );
  if (  dataMode == "" )
    {
      std::cerr << "[ERROR]: please provide a data type, it could be either mc or data. Use --dataMode=<data/mc>" << std::endl;
      return -1;
    }
  
  std::string categoryMode = ParseCommandLine( argc, argv, "-category=" );
  if (  categoryMode == "" )
    {
      std::cerr << "[ERROR]: please provide the category. Use --category=<highpt,highres,lowres>" << std::endl;
      return -1;
    }

  std::string AIC = ParseCommandLine( argc, argv, "-AIC=" );
  if ( AIC == "" )
    {
      std::cerr << "[INFO]: AIC not provided. Use --category=<yes/no>. Using <no> in this execution" << std::endl;
    }
  
  std::string fitMode = ParseCommandLine( argc, argv, "-fitMode=" );
  if (  fitMode == "" )
    {
      std::cerr << "[ERROR]: please provide a fit mode using --fitMode=<mode> , options are: \nsideband\nsb (signal+bkg)\nAIC\nAIC2\nbias\nsignalFit\nbiasSignal\nchooseBinning" << std::endl;
      return -1;
    }

  std::string runPeriod = ParseCommandLine( argc, argv, "-runPeriod=" );
  if (  runPeriod == "" )
    {
      std::cerr << "[ERROR]: please provide a fit mode using --runPeriod=<run1/run2>" << std::endl; 
      return -1;
    }

  std::string f1 = ParseCommandLine( argc, argv, "-f1=" );
  if (  f1 == "" && fitMode == "bias" )
    {
      std::cerr << "[WARNING]: f1 name not provided, using singleExp. Options are: singleExp\ndoubleExp\nsinglePow\ndoublePow";
      std::cerr << "\nmodExp\npoly2\npoly3" << std::endl;
      f1 = "singleExp";
    }
  
  std::string f2 = ParseCommandLine( argc, argv, "-f2=" );
  if (  f2 == "" && fitMode == "bias" )
    {
      std::cerr << "[WARNING]: f2 name not provided, using doubleExp. Options are: singleExp\ndoubleExp\nsinglePow\ndoublePow";
      std::cerr << "\nmodExp\npoly2\npoly3" << std::endl;
      f2 = "doubleExp";
    }

  std::string invertIso = ParseCommandLine( argc, argv, "-invertedIso=" );
  if ( invertIso != "yes" && invertIso != "no" )
    {
      std::cerr << "[WARNING]: invertedIso not provided, provide using: --invertedIso=<yes/no>" << std::endl;
      std::cerr << "[WARNING]: deafault value is <no>" << std::endl;
      invertIso = "no";
    }
  
  std::string _fraction = ParseCommandLine( argc, argv, "-signalFraction=" );
  float _signalFraction = 0.0;
  if ( _fraction != "" )
    {
      _signalFraction = atof( _fraction.c_str() );
    }
  
  if ( fitMode == "biasSignal" && _fraction == "" )
    {
      std::cerr << "[WARNING] signal fraction for bias study is set to default: " << _signalFraction << std::endl;
      std::cerr << "[WARNING] to set a signal fraction use --signalFraction=<sfrac>: " << std::endl;
    }
  
  std::string _ntoys = ParseCommandLine( argc, argv, "-nToys=" );
  int ntoys = 0;
  if ( _ntoys != "" )
    {
      ntoys = atoi( _ntoys.c_str() );
    }
  else
    {
      std::cerr << "[INFO] nToys has not been set, please use --nToys=" << std::endl;
    }

  std::string _ndataset = ParseCommandLine( argc, argv, "-nDataset=" );
  int ndataset = 10;
  if ( _ndataset != "" )
    {
      ndataset = atoi( _ndataset.c_str() );
    }
  else
    {
      std::cerr << "[INFO] nDataset has not been set, please use --nDataset=" << std::endl;
    }
  std::string outputfilename = ParseCommandLine( argc, argv, "-outputFile=" );
  
  std::cout << "[INFO]: tree name is  :" << treeName << std::endl;
  std::cout << "[INFO]: data/mc mode  :" << dataMode << std::endl;
  std::cout << "[INFO]: category mode :" << categoryMode << std::endl;
  std::cout << "[INFO]: fit mode      :" << fitMode << std::endl;
  std::cout << "[INFO]: outputfile    :" << outputfilename << std::endl;
  std::cout << "[INFO]: inverted Iso  :" << invertIso << std::endl;
  std::cout << "[INFO]: nToys         :" << ntoys << std::endl;
  std::cout << "[INFO]: nDataset      :" << ndataset << std::endl;
  
  if ( fitMode == "biasSignal" ) std::cout << "[INFO]: signal Fraction  :" << _signalFraction << std::endl;
  

  //------------------------
  //combine datacard related
  //------------------------
  std::string inputFileSignal = ParseCommandLine( argc, argv, "-inputFileSignal=" );
  if (  inputFileSignal == "" && fitMode == "datacard" )
    {
      std::cerr << "[WARNING]: please provide an input file using --inputFileSignal=<path_to_file>" << std::endl;
    }

  std::string inputFileSMH = ParseCommandLine( argc, argv, "-inputFileSMH=" );
  if (  inputFileSMH == "" && fitMode == "datacard")
    {
      std::cerr << "[WARNING]: please provide an input file using --inputFileSMH=<path_to_file>" << std::endl;
    }

  //SMH nominal yield
  std::string SMH_Yield = ParseCommandLine( argc, argv, "-SMH_Yield=" );
  float _SMH_Yield = 1.e-2;
  if (  SMH_Yield == "" && fitMode == "datacard" )
    {
      std::cerr << "[WARNING]: please provide an input SMH_Yield, --SMH_Yield=<Yield>" << std::endl;
    }
  else
    {
      _SMH_Yield = atof( SMH_Yield.c_str() );
    }

  //SMH Config Line
  std::string SMH_CL = ParseCommandLine( argc, argv, "-SMH_CL=" );
  std::stringstream _SMH_CL;
  if (  SMH_CL == "" && fitMode == "datacard" )
    {
      std::cerr << "[ERROR]: please provide an input SMH_CL, --SMH_CL=<SMH config line>" << std::endl;
      return -1;
    }
  _SMH_CL << SMH_CL;

  //Signal Config Line
  std::string Signal_CL = ParseCommandLine( argc, argv, "-Signal_CL=" );
  std::stringstream _Signal_CL;
  if (  Signal_CL == "" && fitMode == "datacard" )
    {
      std::cerr << "[ERROR]: please provide an input Signal_CL, --Signal_CL=<Signal config line>" << std::endl;
      return -1;
    }
  _Signal_CL << Signal_CL;

  //signal nominal yield
  std::string Signal_Yield = ParseCommandLine( argc, argv, "-Signal_Yield=" );
  float _Signal_Yield = 1.;
  if (  Signal_Yield == "" && fitMode == "datacard" )
    {
      std::cerr << "[WARNING]: please provide an input Signal_Yield, --Signal_Yield=<Yield>" << std::endl;
    }
  else
    {
      _Signal_Yield = atof( Signal_Yield.c_str() );
    }

  std::string binNumber = ParseCommandLine( argc, argv, "-binNumber=" );
  TString _binNumber = "-666";
  if (  binNumber == "" && fitMode == "datacard" )
    {
      std::cerr << "[WARNING]: please provide a binNumber, --binNumber=<binNumber>" << std::endl;
    }
  else
    {
      _binNumber = binNumber;
    }
  

  if (  f1 != "" ) std::cout << "[INFO]: f1    :" << f1 << std::endl;
  if (  f2 != "" ) std::cout << "[INFO]: f2    :" << f2 << std::endl;

  //-------------------------------
  //O u t p u t   R O O T   f i l e
  //-------------------------------
  TFile* f;
  TFile* fs;
  TFile* fsmh;
  TTree* tree;
  TTree* treeSignal;
  TTree* treeSMH;
  bool _getSignal = false;
  if ( dataMode == "data" )
    {
      f = new TFile( inputFile.c_str() , "READ");
      tree = (TTree*)f->Get( treeName.c_str() );
    }
  else if ( dataMode == "mc" )
    {
      f = new TFile( inputFile.c_str() , "READ");
      tree = (TTree*)f->Get( treeName.c_str() );
    }
  else if ( dataMode == "data+signal" )
    {
      f = new TFile( inputFile.c_str() , "READ");
      tree = (TTree*)f->Get( treeName.c_str() );
      fs = new TFile( inputFileSignal.c_str() , "READ");
      treeSignal = (TTree*)fs->Get( treeName.c_str() );
      fsmh = new TFile( inputFileSMH.c_str() , "READ");
      treeSMH = (TTree*)fsmh->Get( treeName.c_str() );
      _getSignal = true;
    }
 
  TTree* mc_tree;
  if ( inputFile2 != "" )
    {
      f = new TFile( inputFile2.c_str(), "READ" );
      mc_tree = (TTree*)f->Get( treeName.c_str() );
    }
  
  //****************************************************
  //MR-Rsq Bin Cut String
  //****************************************************  
  TString BinCutString = "";
  if (  LowMRcut != "" || LowRSQcut != "" )
    {
      BinCutString = " && MR > " + LowMRcut + " && t1Rsq > " + LowRSQcut;
    }
  
  if (  HighMRcut != "" )
    {
      BinCutString += " && MR < " + HighMRcut;
    }
  
  if ( HighRSQcut != "" )
    {
      BinCutString += " && t1Rsq < " + HighRSQcut;
    }
  
  if (categoryMode == "highpt") {
    if (BinSelection == "0") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq < 0.05 ";
    if (BinSelection == "1") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "2") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.10 && t1Rsq < 0.15 ";
    if (BinSelection == "3") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.15 && t1Rsq < 0.20 ";
    if (BinSelection == "4") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.20 ";
    if (BinSelection == "5") BinCutString = " && MR > 200 && MR <= 300 && t1Rsq < 0.05";
    if (BinSelection == "6") BinCutString = " && MR > 200 && MR <= 300 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "7") BinCutString = " && MR > 200 && MR <= 300 && t1Rsq >= 0.10 && t1Rsq < 0.15 ";
    if (BinSelection == "8") BinCutString = " && MR > 200 && MR <= 300 && t1Rsq >= 0.15 ";
    if (BinSelection == "9") BinCutString = " && MR > 300 && MR <= 500 && t1Rsq < 0.05";
    if (BinSelection == "10") BinCutString = " && MR > 300 && MR <= 500 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "11") BinCutString = " && MR > 300 && MR <= 500 && t1Rsq >= 0.10 ";
    if (BinSelection == "12") BinCutString = " && MR > 500 && MR <= 1600 && t1Rsq < 0.05";
    if (BinSelection == "13") BinCutString = " && MR > 500 && MR <= 1600 && t1Rsq >= 0.05 ";
    if (BinSelection == "14") BinCutString = " && MR > 1600  ";
  }
  if (categoryMode == "hbb") {
    if (BinSelection == "0") BinCutString = " && MR > 150 && MR <= 300 && t1Rsq < 0.05 ";
    if (BinSelection == "1") BinCutString = " && MR > 150 && MR <= 300 && t1Rsq >= 0.05 ";
    if (BinSelection == "2") BinCutString = " && MR > 300 && MR <= 3000 ";
  }
  if (categoryMode == "zbb") {
    if (BinSelection == "0") BinCutString = " && MR > 150 && MR <= 450 && t1Rsq < 0.05 ";
    if (BinSelection == "1") BinCutString = " && MR > 150 && MR <= 450 && t1Rsq >= 0.05 ";
    if (BinSelection == "2") BinCutString = " && MR > 450 && MR <= 3000 ";
  }
  if (categoryMode == "highres") {
    if (BinSelection == "0") BinCutString = " && MR > 150 && MR <= 250 && t1Rsq < 0.05 ";
    if (BinSelection == "1") BinCutString = " && MR > 150 && MR <= 250 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "2") BinCutString = " && MR > 150 && MR <= 250 && t1Rsq >= 0.10 && t1Rsq < 0.15 ";
    if (BinSelection == "3") BinCutString = " && MR > 150 && MR <= 250 && t1Rsq >= 0.15 ";
    if (BinSelection == "4") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq < 0.05";
    if (BinSelection == "5") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "6") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq >= 0.10 ";
    if (BinSelection == "7") BinCutString = " && MR > 400 && MR <= 1400 && t1Rsq < 0.05";
    if (BinSelection == "8") BinCutString = " && MR > 400 && MR <= 1400 && t1Rsq >= 0.05  ";
    if (BinSelection == "9") BinCutString = " && MR > 1400 ";
  }
  if (categoryMode == "lowres") {
    if (BinSelection == "0") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq < 0.05 ";
    if (BinSelection == "1") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "2") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.10 && t1Rsq < 0.15 ";
    if (BinSelection == "3") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.15 && t1Rsq < 0.20 ";
    if (BinSelection == "4") BinCutString = " && MR > 150 && MR <= 200 && t1Rsq >= 0.20 ";
    if (BinSelection == "5") BinCutString = " && MR > 200 && MR <= 250 && t1Rsq < 0.05";
    if (BinSelection == "6") BinCutString = " && MR > 200 && MR <= 250 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "7") BinCutString = " && MR > 200 && MR <= 250 && t1Rsq >= 0.10 && t1Rsq < 0.15 ";
    if (BinSelection == "8") BinCutString = " && MR > 200 && MR <= 250 && t1Rsq >= 0.15 ";
    if (BinSelection == "9") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq < 0.05";
    if (BinSelection == "10") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq >= 0.05 && t1Rsq < 0.10 ";
    if (BinSelection == "11") BinCutString = " && MR > 250 && MR <= 400 && t1Rsq >= 0.10 ";
    if (BinSelection == "12") BinCutString = " && MR > 400 && MR <= 1200 && t1Rsq < 0.05";
    if (BinSelection == "13") BinCutString = " && MR > 400 && MR <= 1200 && t1Rsq >= 0.05 ";
    if (BinSelection == "14") BinCutString = " && MR > 1200  ";
  }
 
  //----------------
  // M a i n   C u t
  //----------------
  /*CP's Tree Format is default*/
  
  TString cut = "mGammaGamma >103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 25. && pho2Pt>25.";
  //assymetric cut on photon PT
  //TString cut = "mGammaGamma >103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 40. && pho2Pt>40.";
  //TString cutMETfilters = "&& (Flag_HBHENoiseFilter == 1 && Flag_CSCTightHaloFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1)";
  TString cutMETfilters = "";
  TString cutTrigger = "";
  if(fitMode == "AIC2")
     {
	cutMETfilters = "&& (Flag_HBHENoiseFilter == 1 && Flag_CSCTightHaloFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1)";
	cutTrigger = "&&  HLTDecision[65] == 1";
     }
  //****************************************************
  //Category Cut String
  //****************************************************  
  TString categoryCutString = "";
  if ( treeName == "SusyHggTree" )
    {
      if (categoryMode == "highpt") categoryCutString = " && ptgg >= 110 ";
      else if (categoryMode == "hbb") categoryCutString = " && ptgg < 110 && abs(mbb-125)<25";
      else if (categoryMode == "zbb") categoryCutString = " && ptgg < 110 && abs(mbb-91.2)<25 ";
      else if (categoryMode == "highres") categoryCutString = " && ptgg < 110 && abs(mbb-125)>=25 && abs(mbb-91.2)>=25 && pho1_sigEoE < 0.015 && pho2_sigEoE < 0.015 ";
      else if (categoryMode == "lowres") categoryCutString = " && ptgg < 110  && abs(mbb-125)>=25 && abs(mbb-91.2)>=25 && !(pho1_sigEoE < 0.015 && pho2_sigEoE < 0.015) ";
      else if (categoryMode == "inclusive") categoryCutString = "";
    }
  else if ( runPeriod == "run1" )
    {
      if (categoryMode == "highpt") categoryCutString = " && pTGammaGamma >= 110 ";
      else if (categoryMode == "hbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.) < 15";
      else if (categoryMode == "zbb") categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbZ-91.) < 15 && abs(mbbH-125.) >= 15 )";
      else if (categoryMode == "hzbb") categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbH-125.) < 15 || ( abs(mbbZ-91.) < 15 && abs(mbbH-125.) >= 15 ) )";
      else if (categoryMode == "highres") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.)>=15 && abs(mbbZ-91.)>=15 && (pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015)";
      else if (categoryMode == "lowres") categoryCutString = " && pTGammaGamma < 110  && abs(mbbH-125.)>=15 && abs(mbbZ-91.)>=15 && !(pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015)";
    }
  else
    {
      if (categoryMode == "highpt") categoryCutString = " && pTGammaGamma >= 110 ";
      else if (categoryMode == "hbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH_L-125.) < 15";
      else if (categoryMode == "zbb") categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbZ_L-91.) < 15 && abs(mbbH_L-125.) >= 15 )";
      else if (categoryMode == "hzbb") categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbH_L-125.) < 15 || ( abs(mbbZ_L-91.) < 15 && abs(mbbH_L-125.) >= 15 ) )";
      else if (categoryMode == "highres") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH_L-125.)>=15 && abs(mbbZ_L-91.)>=15 && sigmaMoverM < 0.0085";
      else if (categoryMode == "lowres") categoryCutString = " && pTGammaGamma < 110  && abs(mbbH_L-125.)>=15 && abs(mbbZ_L-91.)>=15 && sigmaMoverM >= 0.0085";
      else if (categoryMode == "inclusive") categoryCutString = "";
    }
  //---------------------------------------------
  // A l e x ' s   T r e e   F o r m at   C a s e
  //---------------------------------------------
  if ( treeName == "SusyHggTree" )
    {
      if ( invertIso != "yes" )
	{
	  cut = "ptgg > 20 && abs(pho1_eta) < 1.48 && abs(pho2_eta) < 1.48 && (pho1_pt>40 || pho2_pt>40) && pho1_pt > 25 && pho2_pt > 25 && pho1_pass_id == 1 && pho1_pass_iso == 1 && pho2_pass_id == 1 && pho2_pass_iso == 1 && mgg > 103 && mgg < 160" + categoryCutString + BinCutString;
	}
      else
	{
	  cut = "ptgg > 20 && abs(pho1_eta) < 1.48 && abs(pho2_eta) < 1.48 && (pho1_pt>40 || pho2_pt>40) && pho1_pt > 25 && pho2_pt > 25 && pho1_pass_id == 1 && pho2_pass_id == 1 && mgg > 103 && mgg < 160 && ( (pho1_pass_iso == 1 && pho2_pass_iso == 0) || (pho1_pass_iso == 0 && pho2_pass_iso == 1) )" + categoryCutString + BinCutString;
	}
    }
  
  float forceSigma = 1.5;
  bool constrainMu = false;
  bool sameMu = false;
  float forceMu    = -1;
  //TString mggName  = "mgg";
  TString mggName  = "mGammaGamma";
  RooWorkspace* w;
  RooWorkspace* w2;

  //---
  //AIC
  //---
  RooWorkspace* w_aic[8];
  double aic[8];
  double aic_2[8];
  double aic_3[8];
  double fitStatus_1[8];
  double fitStatus_2[8];
  double fitStatus_3[8];
  double fitStatus_4[8];
  double Nbkg = 0;
  std::map< std::string, double > aic_map;
  std::map< std::string, double > aic_map_2;
  std::map< std::string, double > aic_map_3;

  std::map< std::string, double > fitStatus_1_map;
  std::map< std::string, double > fitStatus_2_map;
  std::map< std::string, double > fitStatus_3_map;
  std::map< std::string, double > fitStatus_4_map;
  //-----------------------------------
  //C o n c a t e n a t i n g   C u t s
  //-----------------------------------
  cut = cut + categoryCutString + BinCutString + cutMETfilters + cutTrigger;
  std::cout << "[INFO]: cut -> " << cut << std::endl;
  //return -1;
  
  TString outName = outputfilename.substr( 0, outputfilename.find(".root") );
  TFile* fout = 0;
  if (outputfilename == "")
    {
      fout = new TFile( "test_out.root", "recreate" );
    }
  else
    {
      fout = new TFile( outputfilename.c_str(), "recreate" );
    }
  
  if ( fitMode == "sideband" )
    {
      w = MakeSideBandFit( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName );
      MakePlot( tree->CopyTree( cut ),  *w, "sideband_fitpdf_dExp_N1N2", mggName );
      GetIntegral( *w, "sideband_fitpdf_dExp_N1N2", mggName );
    }
  else if ( fitMode == "sb" )
    {
      RooWorkspace* w_sb;
      if ( !_getSignal ) 
	{
	  std::cout << "[INFO]: running S+B fit no SMH" << std::endl;
	  w_sb = MakeSignalBkgFit( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName );
	}
      else
	{
	  std::cout << "[INFO]: running S+B fit with SMH" << std::endl;
	  w_sb = MakeSignalBkgFit( tree->CopyTree( cut ), treeSignal->CopyTree( cut ), treeSMH->CopyTree( cut ), mggName );
	}
      w_sb->Write("w_sb");
    }
  else if ( fitMode == "datacard" )
    {
      RooWorkspace* w_sb;
      std::cout << "calling MakeDataCard" << std::endl;
      w_sb = MakeDataCard( tree->CopyTree( cut ), treeSignal->CopyTree( cut ), treeSMH->CopyTree( cut ), mggName, _SMH_Yield, SMH_CL,
      			   _Signal_Yield, Signal_CL, binNumber, categoryMode );
      std::cout << "finish MakeDataCard" << std::endl;
      w_sb->Write("w_sb");
    }
  else if ( fitMode == "AIC" )
    {
      w_aic[0] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[0], "doubleExp" );
      if( aic_map.find("doubleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doubleExp",aic[0]));
      w_aic[1] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[1], "singleExp" );
      if( aic_map.find("singleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singleExp",aic[1]));
      w_aic[2] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[2], "singlePow" );
      if( aic_map.find("singlePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singlePow",aic[2]));
      w_aic[3] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[3], "doublePow" );
      if( aic_map.find("doublePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doublePow",aic[3]));
      w_aic[4] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[4], "poly2" );
      if( aic_map.find("poly2") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly2",aic[4]));
      w_aic[5] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[5], "poly3" );
      if( aic_map.find("poly3") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly3",aic[5]));
      w_aic[6] = MakeSideBandFitAIC( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[6], "modExp" );
      if( aic_map.find("modExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("modExp",aic[6]));
    }
  else if(fitMode == "AIC2")
    {
      //w_aic[0] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[0], aic_2[0], aic_3[0], "singleExp" );
      //if( aic_map.find("singleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singleExp",aic[0]));

    
      w_aic[0] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[0], aic_2[0], aic_3[0], fitStatus_1[0],fitStatus_2[0], fitStatus_3[0],fitStatus_4[0] , "doubleExp");
      w_aic[1] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[1], aic_2[1], aic_3[1], fitStatus_1[1],fitStatus_2[1], fitStatus_3[1],fitStatus_4[1] , "singleExp");
      w_aic[2] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[2], aic_2[2], aic_3[2], fitStatus_1[2],fitStatus_2[2], fitStatus_3[2],fitStatus_4[2] , "singlePow");
      w_aic[3] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[3], aic_2[3], aic_3[3], fitStatus_1[3],fitStatus_2[3], fitStatus_3[3],fitStatus_4[3] , "doublePow");
      w_aic[4] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[4], aic_2[4], aic_3[4], fitStatus_1[4],fitStatus_2[4], fitStatus_3[4],fitStatus_4[4] , "poly2");
      w_aic[5] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[5], aic_2[5], aic_3[5], fitStatus_1[5],fitStatus_2[5], fitStatus_3[5],fitStatus_4[5] , "poly3");
      w_aic[6] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[6], aic_2[6], aic_3[6], fitStatus_1[6],fitStatus_2[6], fitStatus_3[6],fitStatus_4[6] , "modExp");
      w_aic[7] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName,&Nbkg, aic[7], aic_2[7], aic_3[7], fitStatus_1[7],fitStatus_2[7], fitStatus_3[7],fitStatus_4[7] , "poly4");
      
	double min_aic_tmp = 999.0;

       for(int i=0;i<8;i++)
        {
        if(aic[i]<min_aic_tmp)
                {
                min_aic_tmp=aic[i];
                }
        }

       for(int i=0;i<8;i++)
        {
        if(aic[i]>990.0)
                {
                aic[i]=min_aic_tmp+999.99;
                }
        }
      
     
 
      if( aic_map.find("doubleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doubleExp",aic[0]));
      if( aic_map.find("singleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singleExp",aic[1]));
      if( aic_map.find("singlePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singlePow",aic[2]));
      if( aic_map.find("doublePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doublePow",aic[3]));
      if( aic_map.find("poly2") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly2",aic[4]));
      if( aic_map.find("poly3") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly3",aic[5]));
      if( aic_map.find("modExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("modExp",aic[6]));
      if( aic_map.find("poly4") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly4",aic[7]));
      
      if( aic_map_2.find("doubleExp") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("doubleExp",aic_2[0]));
      if( aic_map_2.find("singleExp") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("singleExp",aic_2[1]));
      if( aic_map_2.find("singlePow") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("singlePow",aic_2[2]));
      if( aic_map_2.find("doublePow") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("doublePow",aic_2[3]));
      if( aic_map_2.find("poly2") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("poly2",aic_2[4]));
      if( aic_map_2.find("poly3") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("poly3",aic_2[5]));
      if( aic_map_2.find("modExp") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("modExp",aic_2[6]));
      if( aic_map_2.find("poly4") == aic_map_2.end() ) aic_map_2.insert( std::pair<std::string, double>("poly4",aic_2[7]));
      
      if( aic_map_3.find("doubleExp") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("doubleExp",aic_3[0]));
      if( aic_map_3.find("singleExp") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("singleExp",aic_3[1]));
      if( aic_map_3.find("singlePow") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("singlePow",aic_3[2]));
      if( aic_map_3.find("doublePow") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("doublePow",aic_3[3]));
      if( aic_map_3.find("poly2") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("poly2",aic_3[4]));
      if( aic_map_3.find("poly3") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("poly3",aic_3[5]));
      if( aic_map_3.find("modExp") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("modExp",aic_3[6]));
      if( aic_map_3.find("poly4") == aic_map_3.end() ) aic_map_3.insert( std::pair<std::string, double>("poly4",aic_3[7]));
 
      if( fitStatus_1_map.find("doubleExp") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("doubleExp",fitStatus_1[0]));
      if( fitStatus_1_map.find("singleExp") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("singleExp",fitStatus_1[1]));
      if( fitStatus_1_map.find("singlePow") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("singlePow",fitStatus_1[2]));
      if( fitStatus_1_map.find("doublePow") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("doublePow",fitStatus_1[3]));
      if( fitStatus_1_map.find("poly2") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("poly2",fitStatus_1[4]));
      if( fitStatus_1_map.find("poly3") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("poly3",fitStatus_1[5]));
      if( fitStatus_1_map.find("modExp") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("modExp",fitStatus_1[6]));
      if( fitStatus_1_map.find("poly4") == fitStatus_1_map.end() ) fitStatus_1_map.insert( std::pair<std::string, double>("poly4",fitStatus_1[7]));
   
      if( fitStatus_2_map.find("doubleExp") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("doubleExp",fitStatus_2[0]));
      if( fitStatus_2_map.find("singleExp") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("singleExp",fitStatus_2[1]));
      if( fitStatus_2_map.find("singlePow") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("singlePow",fitStatus_2[2]));
      if( fitStatus_2_map.find("doublePow") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("doublePow",fitStatus_2[3]));
      if( fitStatus_2_map.find("poly2") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("poly2",fitStatus_2[4]));
      if( fitStatus_2_map.find("poly3") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("poly3",fitStatus_2[5]));
      if( fitStatus_2_map.find("modExp") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("modExp",fitStatus_2[6]));
      if( fitStatus_2_map.find("poly4") == fitStatus_2_map.end() ) fitStatus_2_map.insert( std::pair<std::string, double>("poly4",fitStatus_2[7]));
 
      if( fitStatus_3_map.find("doubleExp") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("doubleExp",fitStatus_3[0]));
      if( fitStatus_3_map.find("singleExp") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("singleExp",fitStatus_3[1]));
      if( fitStatus_3_map.find("singlePow") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("singlePow",fitStatus_3[2]));
      if( fitStatus_3_map.find("doublePow") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("doublePow",fitStatus_3[3]));
      if( fitStatus_3_map.find("poly2") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("poly2",fitStatus_3[4]));
      if( fitStatus_3_map.find("poly3") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("poly3",fitStatus_3[5]));
      if( fitStatus_3_map.find("modExp") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("modExp",fitStatus_3[6]));
      if( fitStatus_3_map.find("poly4") == fitStatus_3_map.end() ) fitStatus_3_map.insert( std::pair<std::string, double>("poly4",fitStatus_3[7]));
      
      if( fitStatus_4_map.find("doubleExp") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("doubleExp",fitStatus_4[0]));
      if( fitStatus_4_map.find("singleExp") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("singleExp",fitStatus_4[1]));
      if( fitStatus_4_map.find("singlePow") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("singlePow",fitStatus_4[2]));
      if( fitStatus_4_map.find("doublePow") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("doublePow",fitStatus_4[3]));
      if( fitStatus_4_map.find("poly2") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("poly2",fitStatus_4[4]));
      if( fitStatus_4_map.find("poly3") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("poly3",fitStatus_4[5]));
      if( fitStatus_4_map.find("modExp") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("modExp",fitStatus_4[6]));
      if( fitStatus_4_map.find("poly4") == fitStatus_4_map.end() ) fitStatus_4_map.insert( std::pair<std::string, double>("poly4",fitStatus_4[7]));
     
    }
  else if ( fitMode == "bias" )
    {
      RooWorkspace* w_bias = DoBiasTest( tree->CopyTree( cut ), mggName, f1, f2, 1e4, ntoys);
      w_bias->Write("w_bias");
    }
  else if ( fitMode == "biasSignal")
    {
      RooWorkspace* w_biasSignal = DoBiasTestSignal( tree->CopyTree( cut ), mggName, f1, f2, ntoys, _signalFraction, outName );
      fout->cd();
      w_biasSignal->Write("w_biasSignal");
    }
  else if ( fitMode == "signalFit" )
    {
      //RooWorkspace* w_sFit = DoubleGausFit( tree->CopyTree( cut ), forceSigma, sameMu, forceMu, mggName );
      RooWorkspace* w_sFit = DoubleCBFit( tree->CopyTree( cut ), mggName, 125., 2. );
      w_sFit->Write("w_sFit");
    }
  else if ( fitMode == "chooseBinning")
    {
      TTree* newTree = tree->CopyTree( cut );
      newTree->Draw("mGammaGamma>>tmp1(38,103.,160)", "weight*(1)", "goff");
      double lumi = 2300.0;
      double kf   = 1.37;
      TH1F* massHisto = (TH1F*)gDirectory->Get("tmp1");
      massHisto->Scale( kf*lumi );
      RooWorkspace* w_binning = SelectBinning( massHisto, mggName, "singleExp", f2, ntoys, ndataset, outName );
      w_binning->Write("w_binning");
    }
  else
    {
      std::cout << "[ERROR]: please provide a valid fitMode!!" << std::endl;
      return -1;
    }
  
  //w->Write("w1");
  //w2->Write("w2");
  if ( fitMode == "AIC" )
    {
      w_aic[0]->Write("w3");
      w_aic[1]->Write("w4");
      w_aic[2]->Write("w5");
      w_aic[3]->Write("w6");
      w_aic[4]->Write("w7");
      w_aic[5]->Write("w8");
      w_aic[6]->Write("w9");
      std::cout.precision(10);
      double min_aic = 99999999;
      std::string func_min = "";
      std::cout << "MIN AIC->" << min_aic << std::endl;
      for ( auto tmp : aic_map )
	{
	  if (  tmp.second < min_aic )
	    {
	      min_aic  = tmp.second;
	      func_min = tmp.first;
	    }
	}
      std::cout << "MIN AIC function is: " << func_min << " AICc is: " << aic_map[func_min] << std::endl;
      //------------------------------
      //c o m p u t e   d e l t a  AIC
      //------------------------------
      std::map<std::string, double> delta_aic_map;
      double sumweights = 0.0;
      for ( auto tmp : aic_map )
	{
	  delta_aic_map[tmp.first] = tmp.second - min_aic;
	  sumweights += TMath::Exp( -0.5*delta_aic_map[tmp.first] ); 
	  //std::cout << tmp.first << " deltaAICc: " << delta_aic_map[tmp.first] << std::endl;
	}
      //------------------------------
      //c o m p u t e   AIC  weights
      //------------------------------
      std::map<std::string, double> aic_weight_map;
      for ( auto tmp : delta_aic_map )
	{
	  aic_weight_map[tmp.first] = TMath::Exp( -0.5*delta_aic_map[tmp.first] )/sumweights;
	  std::cout << tmp.first << " AICc Weights: " << aic_weight_map[tmp.first] << std::endl;
	}
      
    }
  else if( fitMode == "AIC2")
    {
      w_aic[0]->Write("w3");
      w_aic[1]->Write("w4");
      w_aic[2]->Write("w5");
      w_aic[3]->Write("w6");
      w_aic[4]->Write("w7");
      w_aic[5]->Write("w8");
      w_aic[6]->Write("w9");
      w_aic[7]->Write("w10");
      std::cout.precision(10);
      double min_aic = 99999999;
      double min_aic_2 = 99999999;
      double min_aic_3 = 99999999;
      std::string func_min = "";
      std::string func_min_2 = "";
      std::string func_min_3 = "";
      std::cout << "MIN AIC->" << min_aic << std::endl;
      for ( auto tmp : aic_map )
	{
	  if (  tmp.second < min_aic )
	    {
	      min_aic  = tmp.second;
	      func_min = tmp.first;
	    }
	}
      for ( auto tmp : aic_map_2 )
	{
	  if (  tmp.second < min_aic_2 )
	    {
	      min_aic_2  = tmp.second;
	      func_min_2 = tmp.first;
	    }
	}
      for ( auto tmp : aic_map_3 )
	{
	  if (  tmp.second < min_aic_3 )
	    {
	      min_aic_3  = tmp.second;
	      func_min_3 = tmp.first;
	    }
	}
      std::cout << "MIN AIC function is: " << func_min << " AICc is: " << aic_map[func_min] << std::endl;
      //------------------------------
      //c o m p u t e   d e l t a  AIC
      //------------------------------
      std::map<std::string, double> delta_aic_map;
      std::map<std::string, double> delta_aic_map_2;
      std::map<std::string, double> delta_aic_map_3;
      double sumweights = 0.0;
      double sumweights_2 = 0.0;
      double sumweights_3 = 0.0;
      for ( auto tmp : aic_map )
	{
	  delta_aic_map[tmp.first] = tmp.second - min_aic;
	  sumweights += TMath::Exp( -0.5*delta_aic_map[tmp.first] ); 
	  //std::cout << tmp.first << " deltaAICc: " << delta_aic_map[tmp.first] << std::endl;
	}
      for ( auto tmp : aic_map_2 )
	{
	  delta_aic_map_2[tmp.first] = tmp.second - min_aic_2;
	  sumweights_2 += TMath::Exp( -0.5*delta_aic_map_2[tmp.first] ); 
	  std::cout << tmp.first << " deltaAICc2: " << delta_aic_map_2[tmp.first] << std::endl;
	}
      for ( auto tmp : aic_map_3 )
	{
	  delta_aic_map_3[tmp.first] = tmp.second - min_aic_3;
	  sumweights_3 += TMath::Exp( -0.5*delta_aic_map_3[tmp.first] ); 
	  //std::cout << tmp.first << " deltaAICc: " << delta_aic_map[tmp.first] << std::endl;
	}
      //------------------------------
      //c o m p u t e   AIC  weights
      //------------------------------
      std::map<std::string, double> aic_weight_map;
      std::map<std::string, double> aic_weight_map_2;
      std::map<std::string, double> aic_weight_map_3;
      for ( auto tmp : delta_aic_map )
	{
	  aic_weight_map[tmp.first] = TMath::Exp( -0.5*delta_aic_map[tmp.first] )/sumweights;
	  aic_weight_map_2[tmp.first] = TMath::Exp( -0.5*delta_aic_map_2[tmp.first] )/sumweights_2;
	  aic_weight_map_3[tmp.first] = TMath::Exp( -0.5*delta_aic_map_3[tmp.first] )/sumweights_3;
	  std::cout << tmp.first << " AICc Weights: " << aic_weight_map[tmp.first] << std::endl;
	}
      
      PrintAICTable(Nbkg, categoryMode, LowMRcut, HighMRcut, LowRSQcut, HighRSQcut, delta_aic_map,delta_aic_map_2,delta_aic_map_3,aic_weight_map,aic_weight_map_2,aic_weight_map_3,w_aic, fitStatus_1_map, fitStatus_2_map, fitStatus_3_map, fitStatus_4_map);
      //PlotSidebandFit(LowMRcut,LowRSQcut,w_aic);
      
    }

  std::cout << "[INFO]: CLOSING FILE" << std::endl;
  fout->Close();
  gDirectory->Clear();
  gDirectory->Close();
  std::cout << "[INFO]: FILE CLOSED" << std::endl;
  
  return 0;	
}
