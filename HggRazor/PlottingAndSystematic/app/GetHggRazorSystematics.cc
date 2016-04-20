#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <assert.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TROOT.h>
//LOCAL INCLUDES
#include "HggRazorSystematics.hh"
#include "CommandLineInput.hh"

const bool _debug = false;

//---------------
//Binning
//---------------
const int nMR = 4;

//highpt
int highpt_nRsq[nMR] = {7,5,3,2};
float  highptMRedges[] = {150.0, 312.5, 625, 1250, 10000};
float  highptRSQedges0[7] = {0,0.027,0.052,0.077,0.102,0.127,1.0};
float  highptRSQedges1[5] = {0,0.022,0.047,0.072,1.0};
float  highptRSQedges2[3] = {0,0.021,1.0};
float  highptRSQedges3[2] = {0,1.0};

//highres
int highres_nRsq[nMR] = {9,5,3,2};
float  highresMRedges[] = {150,237.5,475,950,10000};
float  highresRSQedges0[9] = {0,0.028,0.053,0.078,0.103,0.128,0.153,0.178,1.0};
float  highresRSQedges1[5] = {0,0.035,0.06,0.085,1.0};
float  highresRSQedges2[3] = {0,0.018,1.0};
float  highresRSQedges3[2] = {0,1.0};

//lowres
int lowres_nRsq[nMR] = {7,6,3,2};
float  lowresMRedges[] = {150,200,400,800,10000};
float  lowresRSQedges0[7] = {0,0.049,0.074,0.099,0.124,0.149,1.0};
float  lowresRSQedges1[6] = {0,0.023,0.048,0.073,0.098,1.0};
float  lowresRSQedges2[3] = {0,0.02,1.0};
float  lowresRSQedges3[2] = {0,1.0};



//-----------------------------------------
//New Binning From Significance Calculation
//-----------------------------------------
//HighPt
float bin_highpt0[4] = {150,0.13,10000,1};
float bin_highpt1[4] = {750,0,10000,0.01};
float bin_highpt2[4] = {750,0.01,10000,0.13};
float bin_highpt3[4] = {500,0,750,0.03};
float bin_highpt4[4] = {500,0.03,750,0.13};
float bin_highpt5[4] = {150,0,400,0.13};
float bin_highpt6[4] = {400,0,500,0.045};
float bin_highpt7[4] = {400,0.045,500,0.13};
std::vector<float*> SetBinning_highpt()
{
  std::vector<float*> myVec;
  myVec.push_back(bin_highpt0);
  myVec.push_back(bin_highpt1);
  myVec.push_back(bin_highpt2);
  myVec.push_back(bin_highpt3);
  myVec.push_back(bin_highpt4);
  myVec.push_back(bin_highpt5);
  myVec.push_back(bin_highpt6);
  myVec.push_back(bin_highpt7);
  return myVec;
};
//HZbb
float bin_hzbb0[4] = {150,0,10000,1};
std::vector<float*> SetBinning_hzbb()
{
  std::vector<float*> myVec;
  myVec.push_back(bin_hzbb0);
  return myVec;
};
//HighRes
float bin_highres0[4] = {600,0.01,10000,1};
float bin_highres1[4] = {150,0.175,600,1};
float bin_highres2[4] = {150,0,400,0.175};
float bin_highres3[4] = {400,0,600,0.025};
float bin_highres4[4] = {400,0.025,600,0.175};
float bin_highres5[4] = {600,0,950,0.01};
float bin_highres6[4] = {950,0,10000,0.01};
std::vector<float*> SetBinning_highres()
{
  std::vector<float*> myVec;
  myVec.push_back(bin_highres0);
  myVec.push_back(bin_highres1);
  myVec.push_back(bin_highres2);
  myVec.push_back(bin_highres3);
  myVec.push_back(bin_highres4);
  myVec.push_back(bin_highres5);
  myVec.push_back(bin_highres6);
  return myVec;
};
//LowRes
float bin_lowres0[4] = {500,0.01,10000,1};
float bin_lowres1[4] = {150,0.15,500,1};
float bin_lowres2[4] = {150,0,400,0.15};
float bin_lowres3[4] = {400,0,500,0.015};
float bin_lowres4[4] = {400,0.015,500,0.15};
float bin_lowres5[4] = {500,0,800,0.01};
float bin_lowres6[4] = {800,0,10000,0.01};
std::vector<float*> SetBinning_lowres()
{
  std::vector<float*> myVec;
  myVec.push_back(bin_lowres0);
  myVec.push_back(bin_lowres1);
  myVec.push_back(bin_lowres2);
  myVec.push_back(bin_lowres3);
  myVec.push_back(bin_lowres4);
  myVec.push_back(bin_lowres5);
  myVec.push_back(bin_lowres6);
  return myVec;
};


//----------------
//Static Variables
//----------------
float HggRazorSystematics::Lumi  = 2300.0;
float HggRazorSystematics::NR_kf = 1.0;
int   HggRazorSystematics::n_PdfSys = 60;

void SetMapBinning( std::map<std::pair<float, float>, std::vector<float>>& myMap, TString category = "highpt" );

int main( int argc, char* argv[] )
{

  gROOT->Reset();
  std::string inputList = ParseCommandLine( argc, argv, "-inputList=" );
  if (  inputList == "" )
    {
      std::cerr << "[ERROR]: please provide an input list file using --inputList=<path_to_list_file>" << std::endl;
      return -1;
    }
  //-----------------
  //Selection TString
  //-----------------
  std::string categoryMode = ParseCommandLine( argc, argv, "-category=" );
  if (  categoryMode == "" )
    {
      std::cerr << "[ERROR]: please provide the category. Use --category=<highpt,hzbb,highres,lowres>" << std::endl;
      return -1;
    }
  
  TString cut = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 25. && pho2Pt>25.";
  TString categoryCutString;

  if (categoryMode == "highpt") categoryCutString = " && pTGammaGamma >= 110 ";
  else if (categoryMode == "hzbb") categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbH_L-125.) < 15. || ( abs(mbbH_L-125.) >= 15. && abs(mbbZ_L-91.) < 15 ) )";
  else if (categoryMode == "highres") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH_L-125.) >= 15 && abs(mbbZ_L-91.) >= 15 && sigmaMoverM < 0.0085";
  else if (categoryMode == "lowres") categoryCutString  = " && pTGammaGamma < 110 && abs(mbbH_L-125.) >= 15 && abs(mbbZ_L-91.) >= 15 && sigmaMoverM >= 0.0085 ";
  else if (categoryMode == "inclusive") categoryCutString = "";

  cut = cut + categoryCutString;

  std::ifstream ifs( inputList, std::ifstream::in );
  assert(ifs);

  std::vector<std::pair<float,float>> facScaleSys;
  std::vector<std::pair<float,float>> renScaleSys;
  std::vector<std::pair<float,float>> facRenScaleSys;

  std::vector<float*> myVectBinning;
  if ( categoryMode == "highpt")
    {
      myVectBinning = SetBinning_highpt();
    }
  else if ( categoryMode == "hzbb" )
    {
      myVectBinning = SetBinning_hzbb();
    }
  else if ( categoryMode == "highres" )
    {
      myVectBinning = SetBinning_highres();
    }
  else if ( categoryMode == "lowres" )
    {
      myVectBinning = SetBinning_lowres();
    }
  else
    {
      std::cerr << "[ERROR]: category is not <highpt/hzbb/highres/lowres>; quitting" << std::endl;
      return -1;
    }
  
  TH2Poly* nominal = new TH2Poly("nominal_SMH", "", 150, 10000, 0, 1 );

  TH2Poly* facScaleUp   = new TH2Poly("facScaleUp", "", 150, 10000, 0, 1 );
  TH2Poly* facScaleDown = new TH2Poly("facScaleDown", "", 150, 10000, 0, 1 );
  
  TH2Poly* renScaleUp   = new TH2Poly("renScaleUp", "", 150, 10000, 0, 1 );
  TH2Poly* renScaleDown = new TH2Poly("renScaleDown", "", 150, 10000, 0, 1 );
  
  TH2Poly* facRenScaleUp   = new TH2Poly("facRenScaleUp", "", 150, 10000, 0, 1 );
  TH2Poly* facRenScaleDown = new TH2Poly("facRenScaleDown", "", 150, 10000, 0, 1 );
  
  TH2Poly* JesUp   = new TH2Poly("JesUp", "", 150, 10000, 0, 1 );
  TH2Poly* JesDown = new TH2Poly("JesDown", "", 150, 10000, 0, 1 );
  
  TH2Poly* pdf[60];
  for ( int i = 0; i < 60; i++ )
    {
      TString pdfName = Form("pdf%d", i);
      pdf[i] = new TH2Poly( pdfName, "", 150, 10000, 0, 1 );
    }
  std::map< std::pair<float,float>, float > sysMap;
  for ( auto tmp : myVectBinning )
    {
      nominal->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      facScaleUp->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      facScaleDown->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      renScaleUp->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      renScaleDown->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      facRenScaleUp->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      facRenScaleDown->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      JesUp->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      JesDown->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
      for( int i = 0; i < 60; i++ ) pdf[i]->AddBin( tmp[0], tmp[1], tmp[2], tmp[3] );
    }
    
  std::string process, rootFileName;
  while ( ifs.good() )
    {
      ifs >> process >> rootFileName;
      if ( ifs.eof() ) break;
      if ( process.find("#") != std::string::npos ) continue;
      if ( _debug ) std::cout << process << " " << rootFileName << std::endl;
      TFile* fin = new TFile( rootFileName.c_str(), "READ");
      //std::cout << "[INFO]: checking file: " << rootFileName << std::endl;
      assert( fin );
      if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;
      
      //------------------------
      //Getting TTree and Histos
      //------------------------
      TTree* tree = (TTree*)fin->Get("HggRazor");
      assert( tree );
      TH1F* NEvents = (TH1F*)fin->Get("NEvents");
      assert( NEvents );
      TH1F* SumScaleWeights   = (TH1F*)fin->Get("SumScaleWeights");
      assert( SumScaleWeights );
      TH1F* SumPdfWeights   = (TH1F*)fin->Get("SumPdfWeights");
      assert( SumPdfWeights );
      
      TFile* tmp = new TFile("tmp.root", "RECREATE");
      TTree* cutTree = tree->CopyTree( cut );
      TString currentProcess = process.c_str();

      //-----------------------------
      //Map container for the binning
      //-----------------------------
      //std::map<std::pair<float, float>, std::vector<float>> binningMap;
      //SetMapBinning( binningMap, categoryMode );
      //---------------------------
      //Create HggSystematic object
      //---------------------------
      HggRazorSystematics* hggSys = new HggRazorSystematics( cutTree, currentProcess, categoryMode, false, false );
      //hggSys->PrintBinning();
      //hggSys->SetBinningMap( binningMap );
      //hggSys->PrintBinning();
      hggSys->SetBinningVector( myVectBinning );
      hggSys->InitMrRsqTH2Poly( 1 );
      hggSys->SetNeventsHisto( NEvents );
      hggSys->SetFacScaleWeightsHisto( SumScaleWeights );
      hggSys->SetPdfWeightsHisto( SumPdfWeights );
      hggSys->Loop();
      for ( auto tmp: myVectBinning )
	{
	  int bin = nominal->FindBin( tmp[0]+10, tmp[1]+0.0001 );
	  nominal->SetBinContent( bin, nominal->GetBinContent(bin) + hggSys->GetNominalYield( tmp[0], tmp[1] ) );
	  //facScale
	  std::pair<float, float> facSys = hggSys->GetFacScaleSystematic( tmp[0], tmp[1] );
	  facScaleUp->SetBinContent( bin, facScaleUp->GetBinContent(bin) + facSys.first );
	  facScaleDown->SetBinContent( bin, facScaleDown->GetBinContent(bin) + facSys.second );
	  //renScale
	  facSys = hggSys->GetRenScaleSystematic( tmp[0], tmp[1] );
	  renScaleUp->SetBinContent( bin, renScaleUp->GetBinContent(bin) + facSys.first );
	  renScaleDown->SetBinContent( bin, renScaleDown->GetBinContent(bin) + facSys.second );
	  //facRenScale
	  facSys = hggSys->GetFacRenScaleSystematic( tmp[0], tmp[1] );
	  facRenScaleUp->SetBinContent( bin, facRenScaleUp->GetBinContent(bin) + facSys.first );
	  facRenScaleDown->SetBinContent( bin, facRenScaleDown->GetBinContent(bin) + facSys.second );
	  //JES
	  facSys = hggSys->GetJesSystematic( tmp[0], tmp[1] );
	  JesUp->SetBinContent( bin, JesUp->GetBinContent(bin) + facSys.first );
	  JesDown->SetBinContent( bin, JesDown->GetBinContent(bin) + facSys.second );
	  //PDF
	  for ( int ipdf = 0; ipdf < 60; ipdf++ )
	    {
	      pdf[ipdf]->SetBinContent( bin, pdf[ipdf]->GetBinContent(bin) + hggSys->GetPdfSystematic( ipdf, tmp[0], tmp[1] ) );
	      //std::cout << "mr: " << tmp[0] << " rsq: " << tmp[1] << "; pdf: " << hggSys->GetPdfSystematic( ipdf, tmp[0], tmp[1] ) << std::endl;
	    }
	  
	}
      
      hggSys->WriteOutput( "histoMR_Rsq" );
      delete hggSys;
      if ( _debug ) std::cout << "deleted hggSys" << std::endl;
      //delete tmp;
      if ( _debug ) std::cout << "deleted tmp File" << std::endl;
    }

  float facScaleTotal[2] = {0,0};
  for ( auto tmp : facScaleSys )
    {
      //std::cout << "Up: " << tmp.first << " , Down: " << tmp.second << std::endl;
      facScaleTotal[0] += tmp.first;
      facScaleTotal[1] += tmp.second;
    } 
  //std::cout << "facScaleUp: " << facScaleTotal[0] << " facScaleDown: " << facScaleTotal[1] << std::endl;


  //std::cout << "#mr_l\tmr_h\trsq_l\trsq_h\tSMHY\tFSU\tFSD\tRSU\tRSD\tFRSU\t\FRSD\tJESU\tJESD\tPDF0\PDF1"<< std::endl;
  std::cout << "#category\t\tmr_l\tmr_h\trsq_l\trsq_h\tSMHY\t\tFSU\t\tFSD";
  for( int ipdf = 0; ipdf < 60; ipdf++ ) std::cout << "\t\tPDF" << ipdf;
  std::cout << "\t\tJESU\t\tJESD"<< std::endl;
  
   for ( auto tmp: myVectBinning )
     {
       int bin   = nominal->FindBin( tmp[0]+10, tmp[1]+0.0001 );
       float nom = nominal->GetBinContent( bin );
       facScaleUp->SetBinContent( bin, facScaleUp->GetBinContent(bin)/nom );
       facScaleDown->SetBinContent( bin, facScaleDown->GetBinContent(bin)/nom );
       
       renScaleUp->SetBinContent( bin, renScaleUp->GetBinContent(bin)/nom );
       renScaleDown->SetBinContent( bin, renScaleDown->GetBinContent(bin)/nom );
       facRenScaleUp->SetBinContent( bin, facRenScaleUp->GetBinContent(bin)/nom );
       facRenScaleDown->SetBinContent( bin, facRenScaleDown->GetBinContent(bin)/nom );
       JesUp->SetBinContent( bin, JesUp->GetBinContent( bin )/nom );
       JesDown->SetBinContent( bin, JesDown->GetBinContent( bin )/nom );

       std::cout << categoryMode << "\t" << tmp[0] << "\t" << tmp[2] << " \t" << tmp[1] << "\t" << tmp[3] << "\t" << nominal->GetBinContent( bin ) << "\t"
		 <<  facScaleUp->GetBinContent( bin ) << "\t" <<  facScaleDown->GetBinContent( bin ) << "\t";
       for( int ipdf = 0; ipdf < 60; ipdf++ )
	 {
	   pdf[ipdf]->SetBinContent( bin, pdf[ipdf]->GetBinContent( bin )/nom );
	   std::cout << "\t" << pdf[ipdf]->GetBinContent( bin );
	 }
       std::cout << "\t" << JesUp->GetBinContent( bin ) << "\t" <<  JesDown->GetBinContent( bin ) << "\t";
       
       std::cout <<  6*nominal->GetBinContent( bin ) << "\t" << facScaleUp->GetBinContent( bin ) << "\t" <<  facScaleDown->GetBinContent( bin ) << "\t";
       for( int ipdf = 0; ipdf < 60; ipdf++ )
	 {
	   pdf[ipdf]->SetBinContent( bin, pdf[ipdf]->GetBinContent( bin )/nom );
	   std::cout << "\t" << pdf[ipdf]->GetBinContent( bin );
	 }
       std::cout << "\t" << JesUp->GetBinContent( bin ) << "\t" <<  JesDown->GetBinContent( bin ) << "\t"
		 << std::endl;
       
       
     }

   
   
  TFile* sF = new TFile( "fullSys.root", "recreate" );
  nominal->Write("SMH_nominal");
  facScaleUp->Write("facScaleUp");
  facScaleDown->Write("facScaleDown");
  renScaleUp->Write("renScaleUp");
  renScaleDown->Write("renScaleDown");
  facRenScaleUp->Write("facRenScaleUp");
  facRenScaleDown->Write("facRenScaleDown");
  JesUp->Write("JesUp");
  JesDown->Write("JesDown");
  for( int ipdf = 0; ipdf < 60; ipdf++ ) pdf[ipdf]->Write();
  sF->Close();
  
  
  /*
  //-----------------------------
  //Map container for the binning
  //-----------------------------
  std::map<std::pair<float, float>, std::vector<float>> binningMap;
  SetMapBinning( binningMap, categoryMode );
  
  //---------------------------
  //Create HggSystematic object
  //---------------------------
  HggRazorSystematics* hggSys = new HggRazorSystematics();
  hggSys->PrintBinning();
  hggSys->SetBinningMap( binningMap );
  hggSys->PrintBinning();
  hggSys->Loop();
  */
  return 0;
}


void SetMapBinning( std::map<std::pair<float, float>, std::vector<float>>& myMap, TString category )
{
  if ( category == "highpt" )
    {
      for ( int i = 0; i < nMR; i++ )
	{
	  std::vector<float> vect;
	  for ( int j = 0; j < highpt_nRsq[i]; j++ )
	    {
	      switch (i)
		{
		case 0:
		  vect.push_back(highptRSQedges0[j]);
		  break;
		case 1:
		  vect.push_back(highptRSQedges1[j]);
		  break;
		case 2:
		  vect.push_back(highptRSQedges2[j]);
		  break;
		case 3:
		  vect.push_back(highptRSQedges3[j]);
		  break;
		default:
		  std::cerr << "MR bins exceeds number of defined Rsq binning!! Do nothing" << std::endl;
		  break;
		}
	    }
	  std::pair<float, float> mypair = std::make_pair( highptMRedges[i], highptMRedges[i+1] );
	  std::cout << mypair.first << " " << mypair.second << std::endl;
	  myMap[mypair] = vect;
	}
    }
  else if ( category == "hzbb" )
    {
      std::vector<float> vect;
      vect.push_back( 0.0 );
      vect.push_back( 1.0 );
      std::pair<float, float> mypair = std::make_pair( 150., 10000. );
      std::cout << mypair.first << " " << mypair.second << std::endl;
      myMap[mypair] = vect;
    }
  else if ( category == "highres" )
    {
      for ( int i = 0; i < nMR; i++ )
	{
	  std::vector<float> vect;
	  for ( int j = 0; j < highres_nRsq[i]; j++ )
	    {
	      switch (i)
		{
		case 0:
		  vect.push_back(highresRSQedges0[j]);
		  break;
		case 1:
		  vect.push_back(highresRSQedges1[j]);
		  break;
		case 2:
		  vect.push_back(highresRSQedges2[j]);
		  break;
		case 3:
		  vect.push_back(highresRSQedges3[j]);
		  break;
		default:
		  std::cerr << "MR bins exceeds number of defined Rsq binning!! Do nothing" << std::endl;
		  break;
		}
	    }
	  std::pair<float, float> mypair = std::make_pair( highresMRedges[i], highresMRedges[i+1] );
	  std::cout << mypair.first << " " << mypair.second << std::endl;
	  myMap[mypair] = vect;
	}
    }
  else if ( category == "lowres" )
    {
      for ( int i = 0; i < nMR; i++ )
	{
	  std::vector<float> vect;
	  for ( int j = 0; j < lowres_nRsq[i]; j++ )
	    {
	      switch (i)
		{
		case 0:
		  vect.push_back(lowresRSQedges0[j]);
		  break;
		case 1:
		  vect.push_back(lowresRSQedges1[j]);
		  break;
		case 2:
		  vect.push_back(lowresRSQedges2[j]);
		  break;
		case 3:
		  vect.push_back(lowresRSQedges3[j]);
		  break;
		default:
		  std::cerr << "MR bins exceeds number of defined Rsq binning!! Do nothing" << std::endl;
		  break;
		}
	    }
	  std::pair<float, float> mypair = std::make_pair( lowresMRedges[i], lowresMRedges[i+1] );
	  std::cout << mypair.first << " " << mypair.second << std::endl;
	  myMap[mypair] = vect;
	}
    }
  return;
};
