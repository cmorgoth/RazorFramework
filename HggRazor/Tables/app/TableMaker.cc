//C++ INCLUDES
#include <iostream>
#include <map>
#include <math.h>
#include <fstream>
#include <iomanip> 
//ROOT INCLUDES
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>
//LOCAL INCLUDES
#include "HggTree.hh"
#include "HggRazorClass.hh"
#include "CommandLineInput.hh"
#include "HggAux.hh"
#include "CommandLineInput.hh"
#include "TChainTools.hh"
//#include "TableMakerAux.hh"

#define _debug 1
#define TYPES 4

// D e f i n e  B i n n i n g
//---------------------------
int HggRazorClass::n_mgg = 57;
float HggRazorClass::mgg_l = 103.;
float HggRazorClass::mgg_h = 160.;
/*
int HggRazorClass::n_mgg = 5;
float HggRazorClass::mgg_l = 117.5;
float HggRazorClass::mgg_h = 132.5;
*/
int HggRazorClass::n_ptgg = 70;
float HggRazorClass::ptgg_l = 20.;
float HggRazorClass::ptgg_h = 720.;

int HggRazorClass::n_mr = 80;
float HggRazorClass::mr_l = 130.;
float HggRazorClass::mr_h = 8130.;

int HggRazorClass::n_rsq = 125;
float HggRazorClass::rsq_l = .0;
float HggRazorClass::rsq_h = 5.0;

std::string plots[TYPES] = { "mgg", "ptgg", "mr", "rsq" };

TString cut = "MR > 0.0 && Rsq > 0.0 && abs( pho1Eta ) < 1.44 && abs( pho2Eta ) < 1.44 && ( pho1Pt > 40. || pho2Pt > 40. ) && pho1Pt > 25. && pho2Pt> 25. && pTGammaGamma>20.";
//TString cut = "MR > 150.0 && t1Rsq > 0.01 && abs( pho1Eta ) < 1.44 && abs( pho2Eta ) < 1.44 && ( pho1Pt > 40. || pho2Pt > 40. ) && pho1Pt > 25. && pho2Pt> 25. && trigger == 1";
//TString mggCut = "mGammaGamma > 117. 5 && mGammaGamma < 132.5";
TString mggCut = "1";

struct FullBkg
{
  TH1F* ggH;
  TH1F* vH;
  TH1F* vbfH;
  TH1F* ttH;
  TH1F* gammaJet;
};

float hpt_k[2]  = { 0.89, 1.03};
float hres_k[2] = { 0.57, 0.67};
float lres_k[2] = { 0.58, 0.77};

int main ( int argc, char* argv[] )
{
  //std::cout.unsetf ( std::ios::floatfield );
  //std::cout.precision(1);
  std::cout << std::fixed;
  std::cout << std::setprecision(2);

  std::cout << "[INFO]: Initializing program" << std::endl;
  std::cout << "[INFO]: Hgg Branching Fraction = " << HggRazorClass::GetHggBF() << std::endl;
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
  HggRazorClass* hggclass;
  
  Yields incYield;
  Yields signalYield;
  
  double k_f = 1.0;
  const double lumi_frac = 0.253; // (5./19.8)
  const int mod = 1; 
  for( const auto& box : Boxes() )
    {
      std::string boxName = GetBoxString( box );
      for( const auto& process : Process() )
	{
	  std::string processName = GetProcessString( process );
	  if ( process == Process::data && run == "run2" )
	    {
	      std::cout << "[INFO]: run2 will skip data category" << std::endl;
	      signalYield.data[0] = 0.0;
	      signalYield.data[1] = 0.0;
	      continue;
	    }
	  if ( process == Process::qcd ) continue;
	  //if ( process == Process::qcd || process == Process::gammaJet || process == Process::diphoton ) continue;
	  std::cout << "====================================" << std::endl;
	  std::cout << "[INFO]: process name: " << processName << std::endl;
	  std::cout << "====================================" << std::endl;
	  k_f = 1.0;
	  if ( run == "run2" && (process == Process::gammaJet || process == Process::diphoton || process == Process::qcd) )
	    {
	      //13/8 TeV
	      k_f = hres_k[mod]*lumi_frac;
	      if ( box == Boxes::HighPt ) k_f = hpt_k[mod]*lumi_frac;
	      if ( box == Boxes::LowRes ) k_f = lres_k[mod]*lumi_frac;
	    }
	  
	  std::cout << "[INFO]: process name: " << processName << std::endl;
	  // R e t r i e v i n g  T r e e
	  //-----------------------------
	  chain   = new TChain( boxName.c_str() );
	  AddTChain( chain, mapList[processName] );
	  
	  //need to create temporary root file to store cutTree
	  TFile* tmp = new TFile("tmp","recreate");
	  //A p p l y i n g  C u t s
	  //------------------------
	  cutTree = (TTree*)chain->CopyTree( cut + " && " + mggCut );
	  if ( cutTree == NULL )
	    {
	      std::cout << "[WARNING]: Empty selected tree: " << boxName << std::endl;
	      continue;
	    }
	  //C r e a t in g   S e l e c t i o n   O b j e c t
	  //------------------------------------------------
	  hggclass = new HggRazorClass( cutTree, processName, boxName, false, false );
	  double err;
	  float n_events = hggclass->GetYields( 250.0, 0.05, 121., 129., err );
	  n_events *= k_f;
	  
	  //[0] -> yields, [1] -> error
	  if ( process == Process::gammaJet )
	    {
	      signalYield.gammaJet[0] = n_events;
	      signalYield.gammaJet[1] = err;
	    }
	  if ( process == Process::diphoton )
	    {
	      signalYield.diphoton[0] = n_events;
	      signalYield.diphoton[1] = err;
	    }
	  if ( process == Process::ttH )
	    {
	      signalYield.ttH[0] = n_events;
	      signalYield.ttH[1] = err;
	    }
	  if ( process == Process::ggH )
	    {
	      signalYield.ggH[0] = n_events;
	      signalYield.ggH[1] = err;
	    }
	  if ( process == Process::vbfH )
	    {
	      signalYield.vbfH[0] = n_events;
	      signalYield.vbfH[1] = err;
	    }
	  if ( process == Process::vH )
	    {
	      signalYield.vH[0] = n_events;
	      signalYield.vH[1] = err;
	    }
	  if ( process == Process::data )
	    {
	      signalYield.data[0] = n_events;
	      signalYield.data[1] = err;
	    }
	  //hggclass->Loop();
	  //hggclass->WriteOutput( boxName );
	}
      double nonRes[2];
      double Res[2];
      double totalMC[2];
      nonRes[0] = signalYield.gammaJet[0] + signalYield.diphoton[0];
      nonRes[1] = sqrt( signalYield.gammaJet[1]*signalYield.gammaJet[1] + signalYield.diphoton[1]*signalYield.diphoton[1] );
      Res[0] = signalYield.ggH[0] + signalYield.vbfH[0] + signalYield.vH[0] + signalYield.ttH[0];
      Res[1] = sqrt( signalYield.ggH[1]*signalYield.ggH[1] + signalYield.vbfH[1]*signalYield.vbfH[1]
		     + signalYield.vH[1]*signalYield.vH[1] + signalYield.ttH[1]*signalYield.ttH[1] );
      totalMC[0] = nonRes[0] + Res[0];
      totalMC[1] = sqrt( nonRes[1]*nonRes[1] + Res[1]*Res[1] );
      std::string tmp_s = boxName+"_yields.txt";
      std::ofstream ofs ( tmp_s.c_str(), std::fstream::out );
      ofs.precision(2);
      ofs << "\\begin{center}\n";
      ofs << "\\captionof{table}{" << boxName << "}\n";
      ofs << "\\begin{tabular}{|c|c|}\n";
      ofs << "\\hline\n";
      ofs << "process & yield \\\\ \\hline\n";
      ofs << "gammaJet & " << signalYield.gammaJet[0] << " \\pm " << signalYield.gammaJet[1] << " \\\\ \\hline\n";
      ofs << "diphoton & " << signalYield.diphoton[0] << " \\pm " << signalYield.diphoton[1] << "\\\\ \\hline\n";
      ofs << "non-resonant& " << nonRes[0] << " \\pm " << nonRes[1] << "\\\\ \\hline\n";
      ofs << "\\hline\n";
      ofs << "ttH & " << signalYield.ttH[0] << " \\pm "<< signalYield.ttH[1] << "\\\\ \\hline\n";
      ofs << "ggH & " << signalYield.ggH[0] << " \\pm " << signalYield.ggH[1] << "\\\\ \\hline\n";
      ofs << "vH & " << signalYield.vH[0] << " \\pm " << signalYield.vH[1] << "\\\\ \\hline\n";
      ofs << "vbfH &" << signalYield.vbfH[0] << " \\pm " << signalYield.vbfH[1] << "\\\\ \\hline\n";
      ofs << "resonant& " << Res[0] << " \\pm " << Res[1] << "\\\\ \\hline\n";
      ofs << "\\hline\n";
      ofs << "total MC & " << totalMC[0] << " \\pm " << totalMC[1] << "\\\\ \\hline\n";
      ofs << "observed &" << signalYield.data[0] << "\\\\ \\hline\n";
      ofs << "\\end{tabular}\n";
      ofs << "\\end{center}\n";
      
      ofs.close();
      std::cout << "[INFO]: Box name: " << boxName << std::endl;
      std::cout << "gammaJet: " << signalYield.gammaJet[0] << " +- " << signalYield.gammaJet[1] << std::endl;
      std::cout << "diphoton: " << signalYield.diphoton[0] << " +- " << signalYield.diphoton[1] << std::endl;
      std::cout << "ttH: " << signalYield.ttH[0] << " +- " << signalYield.ttH[1] << std::endl;
      std::cout << "ggH: " << signalYield.ggH[0] << " +- " << signalYield.ggH[1] << std::endl;
      std::cout << "vH: " << signalYield.vH[0] << " +- " << signalYield.vH[1] << std::endl;
      std::cout << "vbfH: " << signalYield.vbfH[0] << " +-" << signalYield.vbfH[1] << std::endl;
      std::cout << "data: " << signalYield.data[0] << std::endl;
    }

  return 0;
}
