//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
//LOCAL INCLUDES
#include "HggRazorClass.hh"
/*
int   HggRazorClass::n_mgg;
float HggRazorClass::mgg_l;
float HggRazorClass::mgg_h;

int   HggRazorClass::n_ptgg;
float HggRazorClass::ptgg_l;
float HggRazorClass::ptgg_h;
*/
HggRazorClass::HggRazorClass( )
{
};
HggRazorClass::HggRazorClass( TTree* tree ) : HggTree( tree ), _info( false ), _debug( false )
{
  std::cout << "[DEBUG]: n_mgg: " << n_mgg << std::endl;
  InitMggHisto();
  InitPtggHisto();
};

HggRazorClass::HggRazorClass( TTree* tree, TString processName, TString boxName, bool info, bool debug ) : HggTree( tree ), _info( info ), _debug( debug )
{
  //processName
  if ( processName == "" )
    {
      this->processName = "dummy_process";
    }
  else
    {
      this->processName = processName;
    }
  //boxName
  if ( boxName == "" )
    {
      this->boxName = "dummy_box";
    }
  else
    {
      this->boxName = boxName;
    }
  
  InitMggHisto();
  InitPtggHisto();
  InitMrHisto();
  InitRsqHisto();
  InitMrRsqHisto();
};
HggRazorClass::~HggRazorClass()
{
  if ( _debug ) std::cout << "[DEBUG]: Entering Destructor" << std::endl;
  delete h_mgg;
  delete h_ptgg;
  if ( _debug ) std::cout << "[DEBUG]: Finishing Destructor" << std::endl;
};

bool HggRazorClass::InitMggHisto( )
{
  if ( n_mgg == 0 || mgg_l == mgg_h )
    {
      std::cerr << "[ERROR]: mgg histogram parameters not initialized" << std::endl;
      return false;
    }
  std::cout << "[INFO]: Initializing mgg histogram" << std::endl;
  h_mgg = new TH1F( this->processName + "_" + this->boxName +  "_mgg", "mgg", n_mgg, mgg_l, mgg_h );
  return true;
};

bool HggRazorClass::InitPtggHisto( )
{
  if ( n_ptgg == 0 || ptgg_l == ptgg_h )
    {
      std::cerr << "[ERROR]: ptgg histogram parameters not initialized" << std::endl;
      return false;
    }
  h_ptgg = new TH1F( this->processName + "_" + this->boxName +  "_ptgg", "ptgg", n_ptgg, ptgg_l, ptgg_h );
  return true;
};

bool HggRazorClass::InitMrHisto( )
{
  if ( n_mr == 0 || mr_l == mr_h )
    {
      std::cerr << "[ERROR]: mr histogram parameters not initialized" << std::endl;
      return false;
    }
  std::cout << "[INFO]: Initializing mr histogram" << std::endl;
  h_mr = new TH1F( this->processName + "_" + this->boxName +  "_mr", "mr", n_mr, mr_l, mr_h );
  return true;
};

bool HggRazorClass::InitRsqHisto( )
{
  if ( n_rsq == 0 || rsq_l == rsq_h )
    {
      std::cerr << "[ERROR]: rsq histogram parameters not initialized" << std::endl;
      return false;
    }
  std::cout << "[INFO]: Initializing rsq histogram" << std::endl;
  h_rsq = new TH1F( this->processName + "_" + this->boxName +  "_rsq", "rsq", n_rsq, rsq_l, rsq_h );
  return true;
};

//2D Histos
bool HggRazorClass::InitMrRsqHisto( )
{
  if ( ( n_rsq == 0 || rsq_l == rsq_h )  || ( n_mr == 0 || mr_l == mr_h ) )
    {
      std::cerr << "[ERROR]: mr_rsq histogram parameters not initialized" << std::endl;
      return false;
    }
  std::cout << "[INFO]: Initializing mr_rsq histogram" << std::endl;
  h_mr_rsq = new TH2F( this->processName + "_" + this->boxName +  "_mr_rsq", "mr_rsq", n_mr, mr_l, mr_h, n_rsq, rsq_l, rsq_h );
  return true; 
};
void HggRazorClass::Loop()
{
  if ( _debug ) std::cout << "[DEBUG]: Entering Loop" << std::endl;
  if (fChain == 0) return;
  if ( h_mgg == NULL || h_ptgg == NULL )
    {
      std::cerr << "[ERROR]: Histograms are nor initialized" << std::endl;
      return;
    }
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++ )
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      h_mgg->Fill( mGammaGamma, xsecSF );
      h_ptgg->Fill( pTGammaGamma, xsecSF );
      h_mr->Fill( MR, xsecSF );
      h_rsq->Fill( Rsq, xsecSF );
      h_mr_rsq->Fill( MR, Rsq, xsecSF );
    }
  if ( _debug ) std::cout << "[DEBUG]: Finishing Loop" << std::endl;
};

bool HggRazorClass::WriteOutput( TString outName )
{
  if ( _debug ) std::cout << "[DEBUG]: Entering WriteOutput" << std::endl;
  this->fout = new TFile( outName + "_" + this->processName + ".root", "recreate");
  h_mgg->Write("mgg");
  h_ptgg->Write("ptgg");
  h_mr->Write("mr");
  h_rsq->Write("rsq");
  h_mr_rsq->Write("mr_rsq");
  fout->Close();
  if ( _debug ) std::cout << "[DEBUG]: Finishing WriteOutput" << std::endl;
  return true;
};
