//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
//LOCAL INCLUDES
#include "RunTwoFitMgg.hh"
#include "CommandLineInput.hh"
#include "PrintAIC.hh"

int main( int argc, char* argv[])
{
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

  std::string MRcut = ParseCommandLine( argc, argv, "-MRcut=" );
  std::string RSQcut = ParseCommandLine( argc, argv, "-RSQcut=" );
  std::string BinSelection = ParseCommandLine( argc, argv, "-BinSelection=" );
  if (  BinSelection == "" && ( MRcut == "" || RSQcut == "" ) )
    {
      std::cout << "[ERROR]: user must provide either a BinSelection or an MR cut and Rsq cut. use --BinSelection=<yourBin> or --MRcut=<yourMRCut> --RSQcut=<yourRSQCut>" << std::endl;
      return -1;
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
      std::cerr << "[ERROR]: please provide a fit mode using --fitMode=<mode> , options are: \nsideband\nsb (signal+bkg)\nAIC\nbias\nsignalFit\nbiasSignal" << std::endl;
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
  
  std::string outputfilename = ParseCommandLine( argc, argv, "-outputFile=" );
  
  std::cout << "[INFO]: tree name is  :" << treeName << std::endl;
  std::cout << "[INFO]: data/mc mode  :" << dataMode << std::endl;
  std::cout << "[INFO]: category mode :" << categoryMode << std::endl;
  std::cout << "[INFO]: fit mode      :" << fitMode << std::endl;
  std::cout << "[INFO]: outputfile    :" << outputfilename << std::endl;
  std::cout << "[INFO]: inverted Iso  :" << invertIso << std::endl;
  if ( fitMode == "biasSignal" ) std::cout << "[INFO]: signal Fraction  :" << _signalFraction << std::endl;
  
  if (  f1 != "" ) std::cout << "[INFO]: f1    :" << f1 << std::endl;
  if (  f2 != "" ) std::cout << "[INFO]: f2    :" << f2 << std::endl;
  TFile* f;
  TTree* tree;
  if ( dataMode == "data" )
    {
      f = new TFile( inputFile.c_str() , "READ");
      //f = new TFile( inputFile.c_str() , "UPDATE");
      tree = (TTree*)f->Get( treeName.c_str() );
    }
  else if ( dataMode == "mc" )
    {
      f = new TFile( inputFile.c_str() , "READ");
      tree = (TTree*)f->Get( treeName.c_str() );
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
  if (  MRcut != "" || RSQcut != "" )
    {
      BinCutString = " && MR > " + MRcut + " && t1Rsq > " + RSQcut;
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
  
  TString cut = "pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 25. && pho2Pt>25.";

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
  else
    {
      if (categoryMode == "highpt") categoryCutString = " && pTGammaGamma >= 110 ";
      else if (categoryMode == "hbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.)<25";
      else if (categoryMode == "zbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbZ-91.2)<25 ";
      else if (categoryMode == "highres") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.)>=25 && abs(mbbZ-91.2)>=25 && pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015 ";
      else if (categoryMode == "lowres") categoryCutString = " && pTGammaGamma < 110  && abs(mbbH-125.)>=25 && abs(mbbZ-91.2)>=25 && !(pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015) ";
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
  std::map< std::string, double > aic_map;
  std::map< std::string, double > aic_map_2;
  std::map< std::string, double > aic_map_3;

  //-----------------------------------
  //C o n c a t e n a t i n g   C u t s
  //-----------------------------------
  cut = cut + categoryCutString + BinCutString;
  std::cout << "[INFO]: cut -> " << cut << std::endl;
  //return -1;
  //-------------------------------
  //O u t p u t   R O O T   f i l e
  //-------------------------------
  TFile* fout = 0;
  if (outputfilename == "")
    {
      fout = new TFile( "test_out.root", "recreate" );
    }
  else
    {
      fout = new TFile( outputfilename.c_str(), "recreate" );
    }
  
  TString outName = outputfilename.substr( 0, outputfilename.find(".root") );
  
  if ( fitMode == "sideband" )
    {
      w = MakeSideBandFit( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName );
      MakePlot( tree->CopyTree( cut ),  *w, "sideband_fitpdf_dExp_N1N2", mggName );
      GetIntegral( *w, "sideband_fitpdf_dExp_N1N2", mggName );
    }
  else if ( fitMode == "sb" )
    {
      RooWorkspace* w_sb = MakeSignalBkgFit( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName );
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
      w_aic[0] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[0], aic_2[0], aic_3[0], "doubleExp" );
      if( aic_map.find("doubleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doubleExp",aic[0]));
      w_aic[1] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[1], aic_2[1], aic_3[1], "singleExp" );
      if( aic_map.find("singleExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singleExp",aic[1]));
      w_aic[2] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[2], aic_2[2], aic_3[2], "singlePow" );
      if( aic_map.find("singlePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("singlePow",aic[2]));
      w_aic[3] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[3], aic_2[3], aic_3[3], "doublePow" );
      if( aic_map.find("doublePow") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("doublePow",aic[3]));
      w_aic[4] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[4], aic_2[4], aic_3[4], "poly2" );
      if( aic_map.find("poly2") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly2",aic[4]));
      w_aic[5] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[5], aic_2[5], aic_3[5], "poly3" );
      if( aic_map.find("poly3") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("poly3",aic[5]));
      w_aic[6] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[6], aic_2[6], aic_3[6], "modExp" );
      if( aic_map.find("modExp") == aic_map.end() ) aic_map.insert( std::pair<std::string, double>("modExp",aic[6]));
      w_aic[7] = MakeSideBandFitAIC_2( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName, aic[7], aic_2[7], aic_3[7], "poly4" );
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
    }
  else if ( fitMode == "bias" )
    {
      RooWorkspace* w_bias = DoBiasTest( tree->CopyTree( cut ), mggName, f1, f2, 1e4, 1e4);
      w_bias->Write("w_bias");
    }
  else if ( fitMode == "biasSignal")
    {
      RooWorkspace* w_biasSignal = DoBiasTestSignal( tree->CopyTree( cut ), mggName, f1, f2, 1e4, _signalFraction, outName );
      fout->cd();
      w_biasSignal->Write("w_biasSignal");
    }
  else if ( fitMode == "signalFit" )
    {
      RooWorkspace* w_sFit = DoubleGausFit( tree->CopyTree( cut ), forceSigma, constrainMu, forceMu, mggName );
      w_sFit->Write("w_sFit");
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

     PrintAICTable(MRcut, RSQcut,delta_aic_map,delta_aic_map_2,delta_aic_map_3,aic_weight_map,aic_weight_map_2,aic_weight_map_3,w_aic);
//     PlotSidebandFit(MRcut,RSQcut,w_aic);

}
  
  fout->Close();
    
  return 0;	
}
