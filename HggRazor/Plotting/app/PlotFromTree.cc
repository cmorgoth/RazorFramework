//C++ INCLUDES
#include <iostream>
#include <map>
//ROOT INCLUDES
#include <TString.h>
#include <THStack.h>
//LOCAL INCLUDES
#include "HggTree.hh"
#include "HggRazorClass.hh"
#include "HggAux.hh"
#include "CommandLineInput.hh"
#include "TChainTools.hh"
#include "PlotCosmetics.hh"

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

//C u s t o m MR_RSQ_BINNING
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


//A p p l y   B a s e l i n e   C u t
//-----------------------------------
TString cut = "MR > 250.0 && Rsq > 0.04 && abs( pho1Eta ) < 1.44 && abs( pho2Eta ) < 1.44 && ( pho1Pt > 40. || pho2Pt > 40. ) && pho1Pt > 25. && pho2Pt> 25.";
//TString mggCut = "mGammaGamma > 117. 5 && mGammaGamma < 132.5";
TString mggCut = "1";

#define _debug 1

int main ( int argc, char* argv[] )
{
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
  THStack* stack;
  TLegend* leg;
  TH1F* data;
  TH1F* mc;

  const int nprocesses = 8;
  const int nplots = 4;
  for( const auto& box : Boxes() )
    {
      std::string boxName = GetBoxString( box );
      int ctr = 0;
      Histos histos[nprocesses];
      for( const auto& process : Process() )
	{
	  std::string processName = GetProcessString( process );
	  if ( process == Process::data && run == "run2" )
	    {
	      std::cout << "[INFO]: run2 will skip data category" << std::endl;
	      continue;
	    }
	  std::cout << "====================================" << std::endl;
	  std::cout << "[INFO]: process name: " << processName << std::endl;
	  std::cout << "====================================" << std::endl;
     	  
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
	  if ( box == Boxes::HighPt )
	    {
	      hggclass->InitMrRsqCustomHisto( N_HighPt, MR_HighPt, N_HighPt, Rsq_HighPt );
	    }
	  else if ( box == Boxes::HighRes || box == Boxes::LowRes )
	    {
	      hggclass->InitMrRsqCustomHisto( N_HighRes, MR_HighRes, N_HighRes, Rsq_HighRes );
	    }
	  else if( box == Boxes::Hbb || box == Boxes::Zbb )
	    {
	      hggclass->InitMrRsqCustomHisto( N_Hbb, MR_Hbb, N_Hbb, Rsq_Hbb );
	    }
	  else
	    {
	      std::cout << "[WARNING]: Undefined box->" << std::endl;
	    }
	  //Creating Histos
	  hggclass->Loop();
	  //Assigning Histos
	  histos[ctr].process = process;
	  for ( const auto& htmp : HistoTypes() ) histos[ctr].AssignHisto( htmp, hggclass->GetHisto( htmp ) );
	  ctr++;
	  hggclass->WriteOutput( boxName );
	  //Cosmetics
	}
      //Plotting
      for ( const auto& htmp : HistoTypes() )
	{
	  stack = new THStack( "hs" , "Hgg Stack " );
	  leg = new TLegend( 0.7, 0.58, 0.95, 0.89, NULL, "brNDC" );
	  for (  int i  = 0; i < nprocesses; i++ )
	    {
	      TH1F* tmp_h = new TH1F( histos[i].GetHisto( htmp ) );
	      TH1F* h_s = GetStyledHisto( tmp_h, histos[i].process );
	      if ( histos[i].process == Process::data )
		{
		  data = new TH1F ( *h_s ); 
		}
	      else
		{
		  stack->Add( h_s, "histo" );
		  if ( mc == NULL )
		    {
		      mc = new TH1F( *h_s );
		    }
		  else
		    {
		      mc->Add( h_s );
		    }
		}
	      std::cout << i << " " << GetProcessString( histos[i].process ) << std::endl;
	      AddLegend( h_s, leg, histos[i].process );
	    }
	  std::string histoName = GetHistoTypesString( htmp );
	  std::cout << "making " <<  histoName << std::endl;
	  MakeStackPlot( stack, data, mc, histoName, histoName + "_" + boxName, leg );
	  std::cout << "leaving " <<  histoName << std::endl;
	  delete leg;
	  delete stack;
	  delete mc;
	  delete data;
	}
    }
  return 0;
}
