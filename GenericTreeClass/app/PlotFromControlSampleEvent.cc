//C++ INCLUDES
#include <iostream>
#include <map>
//ROOT INCLUDES
#include <TString.h>
#include <THStack.h>
//LOCAL INCLUDES
#include "ControlSampleEvent.hh"
#include "GenericTreeClass.hh"
#include "HggAux.hh"
#include "CommandLineInput.hh"
#include "TChainTools.hh"
#include "PlotCosmetics.hh"

//--------------------------
//HighPt
const int N_HighPt = 5;
float MR_HighPt[N_HighPt+1] = {150,200,300,500,1600,3000};
float Rsq_HighPt[N_HighPt+1] = {0.01,0.05,0.10,0.15,0.20,1.00};

//HighRes
const int N_HighRes = 4;
float MR_HighRes[N_HighRes+1] = {150,250,400,1400,3000};
float Rsq_HighRes[N_HighRes+1] = {0.01,0.05,0.10,0.15,1.00};

//Hbb
const int N_Hbb = 2;
float MR_Hbb[N_Hbb+1] = {150,300,3000};
float Rsq_Hbb[N_Hbb+1] = {0.01,0.05,1.00};

float hpt_k[2]  = { 1.0, 1.0};
float hres_k[2] = { 1.0, 1.0};
float lres_k[2] = { 1.0, 1.0};
//A p p l y   B a s e l i n e   C u t
//-----------------------------------
//DATA
//EBEB
TString cut = "pho1passEleVeto == 0 && pho2passEleVeto == 0 && pho1passIso == 1 && pho2passIso == 1 && abs( pho1Eta ) < 1.479 && abs( pho2Eta ) < 1.479 && ( HLTDecision[30] == 1 || HLTDecision[31] == 1 ) && mGammaGamma>60. && mGammaGamma<120. && pho1Pt>30 && pho2Pt>20";

//MC
//EBEB
TString cut_mc = "pho1passEleVeto == 0 && pho2passEleVeto == 0 && pho1passIso == 1 && pho2passIso == 1 && abs( pho1Eta ) < 1.479 && abs( pho2Eta ) < 1.479  && ( HLTDecision[30] == 1 || HLTDecision[31] == 1 ) && mGammaGamma>60. && mGammaGamma<120. && pho1Pt>30 && pho2Pt>20";

TString mggCut = "1";

#define _debug 1


struct binning
{
  int nbins;
  float x_low;
  float x_high;
  binning(){};
  binning( int nbins, float x_low, float x_high ): nbins(nbins), x_low(x_low), x_high(x_high) {};
};

int main ( int argc, char* argv[] )
{
  std::cout << "[INFO]: Initializing program" << std::endl;
    
  //Map Containing the lists for different processes 
  std::map< std::string, std::string > mapList;
  std::string inputFile = ParseCommandLine( argc, argv, "-inputFile=" );
  if (  inputFile == "" )
    {
      std::cerr << "[ERROR]: please provide an input file using --inputFile=<path_to_file>" << std::endl;
      return -1;
    }
  
  std::string run = ParseCommandLine( argc, argv, "-run=" );
  if (  run == "" )
    {
      std::cout << "[WARNING]: please provide a valid run, use --run=<run1/run2>" << std::endl;
      run = "run1";
    }

  std::string treeType = ParseCommandLine( argc, argv, "-treeType=" );
  if (  treeType == "" )
    {
      std::cout << "[WARNING]: please provide a valid treeType, use --treeType=<inclusive/category>" << std::endl;
      treeType = "inclusive";
    }
  
  std::cout << "=================================" << std::endl;
  std::cout << "===========set parameters========" << std::endl;
  std::cout << "[INFO]: input file: " << inputFile << std::endl;
  std::cout << "[INFO]: run: " << run << std::endl;
  std::cout << "=================================" << std::endl;
  
  FillMapList( mapList, inputFile );
  if ( _debug ) std::cout << "[DEBUG]: map size: " << mapList.size() << std::endl;
  
  for( auto& myMap : mapList )
    {
      if ( _debug ) std::cout << "[DEBUG]: first: " << myMap.first << " second: " << myMap.second << std::endl;
    }
  
  TFile* f;
  TTree* tree;
  TChain* chain;
  TTree* cutTree;
  GenericTreeClass* myClass = new GenericTreeClass();
  THStack* stack;
  TLegend* leg;
  TH1F* data;
  TH1F* mc;
  TH1F* mc2 = new TH1F();

  const int nprocesses = 2;
  const int nplots = 4;
  double k_f = 1.0;
  const double lumi_frac = 1.0; // (5./19.8)
  const int mod = 0; 
  

  TString plots[3] = {"MR", "Rsq", "NJets40"};
  for ( int i = 0; i < 3; i++)
    {
      myClass->CreateGenericHisto( plots[i], plots[i], 25, 0, 2500 );
    }
  
  myClass->PrintStoredHistos();
  
  myClass->Loop();
  for ( int i = 0; i < 3; i++)
    {
      //Loop the processes in here (not implemented)
      std::cout << "deb 0" << std::endl;
      stack = new THStack( "hs" , "Hgg Stack " );
      leg = new TLegend( 0.7, 0.58, 0.93, 0.89, NULL, "brNDC" );
      std::cout << "deb 2" << std::endl;
      TH1F* tmp_h = new TH1F( *(myClass->map_1D_Histos[ plots[i] ]) );
      std::cout << "deb 2.1" << std::endl;
      TH1F* h_s = GetStyledHisto( tmp_h, Process::gammaJet );
      std::cout << "deb 2" << std::endl;
      stack->Add( h_s, "histo" );
      TH1F* h_data = GetStyledHisto( tmp_h, Process::data );
      data = new TH1F ( *h_data );
      mc = new TH1F( *tmp_h );
      std::cout << "deb 0" << std::endl;
      AddLegend( h_s, leg, Process::gammaJet );
      MakeStackPlot( stack, data, mc, plots[i], "plots/" + plots[i] + "_" + "INCLUSIVE", leg );
    }
  return 0;
}
