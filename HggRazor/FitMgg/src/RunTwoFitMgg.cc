//C++ INCLUDES
#include <vector>
#include <fstream>
//ROOT INCLUDES
//#include <TSYSTEM.h>
#include <TSystem.h>
#include <TTree.h>
#include <TLatex.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBox.h>
#include <TMath.h>
#include <TROOT.h>
#include <Math/GaussIntegrator.h>
#include <Math/IntegratorOptions.h>
//ROOFIT INCLUDES
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooExtendPdf.h>
#include <RooStats/SPlot.h>
#include <RooStats/ModelConfig.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooBernstein.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
//#include <RealVar.h>
//LOCAL INCLUDES
#include "RunTwoFitMgg.hh"
#include "DefinePdfs.hh"

RooWorkspace* DoubleGausFit( TTree* tree, float forceSigma, bool sameMu, float forceMu, TString mggName )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(57);
  //mgg.setRange( "signal", 110, 140. );
  mgg.setRange( "signal", 103, 160. );

  /*
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 100, 1000, "GeV" );
  mgg.setBins(180);
  mgg.setRange( "signal", 600, 850. );
  */
  
  //ws->var("dGauss_signal_gauss_Ns")->setVal( 1600 );
  
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import( *tree ) );
  int npoints = data.numEntries();
  //-----------------------
  //C r e a t e  doubleGaus
  //-----------------------
  TString tag;
  if( sameMu )
    {
      tag = MakeDoubleGauss( "dGauss_signal", mgg, *ws );
      ws->var("dGauss_signal_gauss_Ns")->setVal( (double)npoints );
    }
  else
    {
      tag = MakeFullDoubleGauss( "dGauss_signal", mgg, *ws );
      ws->var("dGauss_signal_DG_Ns")->setVal( (double)npoints );
    }
  
  //ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended( kTRUE ), RooFit::Range("signal") );
  RooFitResult* sres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("signal") );
  
  sres->SetName( tag + "_sres" );
  ws->import( *sres );
  
  RooPlot* frame = mgg.frame();
  data.plotOn( frame );
  ws->pdf( tag )->plotOn( frame, RooFit::LineColor( kBlue ), RooFit::Range("Full"), RooFit::NormRange("Full") );
  ws->import( mgg );
  ws->import( data );
  frame->SetName( tag + "_frame" );
  ws->import( *frame );

  return ws;
};


RooWorkspace* DoubleCBFit( TTree* tree, TString mggName, float muCB, float sigmaCB )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(57);
  mgg.setRange( "signal", 103, 160. );

  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import( *tree ) );
  int npoints = data.numEntries();

  //---------------------------------
  //C r e a t e  Double Crystall Ball
  //---------------------------------
  TString tag = MakeDoubleCB( "Signal", mgg, *ws );
  ws->var( tag + "_Ns" )->setVal( (double)npoints );
  
  RooFitResult* sres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("signal") );
  
  sres->SetName( tag + "_sres" );
  ws->import( *sres );
  
  gSystem->Load("~/Software/git/RazorFramework/HggRazor/FitMgg/CustomPdfs.so");
  RooPlot* frame = mgg.frame();
  data.plotOn( frame );
  ws->pdf( tag )->plotOn( frame, RooFit::LineColor( kBlue ), RooFit::Range("signal"), RooFit::NormRange("signal") );
  ws->import( mgg );
  ws->import( data );
  frame->SetName( tag + "_frame" );
  ws->import( *frame );
  
  return ws;
};



RooWorkspace* MakeSideBandFit( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(57);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  mgg.setRange("signal", 103, 160);
  
  RooRealVar w( "xsecSF", "w", 0, 10000 );

  //TString tag3 = MakeDoubleExp( "sideband_fit", mgg, *ws );
  TString tag3 = MakeDoubleExpN1N2( "sideband_fit", mgg, *ws ); 
  
  //RooDataSet data( "data", "", RooArgSet(mgg, w), RooFit::WeightVar(w), RooFit::Import(*tree) );
  //Sideband Fit
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  ws->pdf( tag3 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  RooFitResult* bres = ws->pdf( tag3 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  //Full Fit
  //ws->pdf( tag3 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //RooFitResult* bres = ws->pdf( tag3 )->fitTo( data, RooFit::Strategy(2), RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Range("Full") );

  mgg.setRange("sregion", 122.04, 128.96);
  RooAbsReal* sint = ws->pdf( tag3 )->createIntegral( mgg, RooFit::NormSet(mgg), RooFit::Range("sregion") ) ; 
  RooAbsReal* Nfit = ws->pdf( tag3 )->createIntegral( mgg, RooFit::Range("Full") ) ;
  std::cout << sint->getVal() << std::endl;
  float N_sideband = data.sumEntries(Form("(mgg>%0.2f && mgg <120) || (mgg>131 && mgg<%0.2f)",103.,160.));
  /*
  double n1 = ws->var("sideband_fitNbkg1")->getVal();
  double n2 = ws->var("sideband_fitNbkg2")->getVal();
  double Nbkg = n1*n1 + n2*n2; 

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Ntotal: " << data.sumEntries( ) << std::endl;
  std::cout << "Nfit: " << Nfit->getVal() << std::endl;
  std::cout << "Nbkg: " << Nbkg << std::endl;
  std::cout << "sf: "  << Nbkg*sint->getVal()/N_sideband << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  */
  bres->SetName( tag3 + "_b_fitres" );
  ws->import( *bres );

  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  ws->pdf( tag3 )->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  ws->pdf( tag3 )->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  //ws->pdf( tag3 )->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("Full"));
  fmgg->SetName( tag3 + "_frame" );
  ws->import( *fmgg );

  RooPlot* pdfFrame = mgg.frame();
  ws->pdf( tag3 )->plotOn( pdfFrame, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("Full") );
  pdfFrame->SetName( tag3+"_pdfframe" );
  ws->import( *pdfFrame );
  ws->import( mgg );
  return ws;
};

void MakePlot( TTree* tree,  RooWorkspace& w, TString pdfName, TString mggName )
{
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",100,400,"GeV");
  mgg.setBins(30);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  
  mgg.setRange("low_v2", 100, 120);
  mgg.setRange("high_v2", 131, 400);
  
  mgg.setRange("signal", 103, 160);
  mgg.setRange("sig", 120., 130.);

  TString tag3 = MakeDoubleExpN1N2( "sideband_fit_v2", mgg, w );
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  //fullFit
  //w.pdf( tag3 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //RooFitResult* bres = w.pdf( tag3 )->fitTo( data, RooFit::Strategy(2), RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //sidebandFit
  w.pdf( tag3 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low_v2,high_v2") );
  RooFitResult* bres = w.pdf( tag3 )->fitTo( data, RooFit::Strategy(2), RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Range("low_v2,high_v2") );
  
  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  //w.pdf( pdfName )->plotOn( fmgg, RooFit::LineColor(kBlue), RooFit::Range("low,high"), RooFit::NormRange("low,high") );
  
  //w.pdf( pdfName )->plotOn( fmgg, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("Full"), RooFit::NormRange("low_v2,high_v2") );
  w.pdf( pdfName )->plotOn( fmgg, RooFit::LineColor(kGreen), RooFit::Range("Full"), RooFit::NormRange("low,high") );
  //w.pdf( tag3 )->plotOn( fmgg, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("low_v2,high_v2") );
  fmgg->SetName( "pdf_goldenbin_test" );
  w.import( *fmgg );
  
  return;
};


double GetIntegral( RooWorkspace& w, TString pdfName, TString mggName )
{
  RooAbsPdf* NewModel = w.pdf( pdfName );
  RooRealVar* mgg = w.var( mggName );
  mgg->setRange("sig", 122.08, 128.92);
  //mgg->setRange("sig", 103.0, 160.);
  RooAbsReal* fIntegral = NewModel->createIntegral(*mgg);
  RooAbsReal* fIntegral2 = NewModel->createIntegral(*mgg, RooFit::NormSet(*mgg), RooFit::Range("sig") );
  std::cout << "test Int: " << fIntegral->getVal() << std::endl;
  std::cout << "test Int2: " << fIntegral2->getVal() << std::endl;
  mgg->setMax(400);
  mgg->setMin(100);
  mgg->setRange("sig", 131, 400);
  fIntegral2 = NewModel->createIntegral(*mgg, RooFit::NormSet(*mgg), RooFit::Range("sig") );
  std::cout << "test Int2': " << fIntegral2->getVal() << std::endl;
  return fIntegral2->getVal();
};

RooWorkspace* MakeSignalBkgFit( TTree* treeData, TTree* treeSignal, TTree* treeSMH, TString mggName )
{
  //------------------------------------------------
  // C r e a t e   s i g n a l  s h a p e from TTree
  //------------------------------------------------
  //TFile* ftmp = new TFile("tmp_output_OurID.root", "recreate");
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(38);
  mgg.setRange( "signal", 115, 135. );
  mgg.setRange( "high", 135, 160);
  mgg.setRange( "low", 103, 120);
  mgg.setRange( "full", 103, 160 );
  //--------------------------------
  //I m p or t i n g   D a t a
  //--------------------------------

  //Getting signal shape from signal MC
  //-----------------------
  //C r e a t e  doubleGaus
  //-----------------------
  bool sameMu = false;
  TString tagSignal, tagSMH;
  TString tag;
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*treeData) );
  RooDataSet dataSignal( "dataSignal", "", RooArgSet(mgg), RooFit::Import(*treeSignal) );
  RooDataSet dataSMH( "dataSMH", "", RooArgSet(mgg), RooFit::Import(*treeSMH) );
  //---------------------------------
  //D e f i n e   s i g n a l   P D F
  //---------------------------------
  int npoints = dataSignal.numEntries();
  if( sameMu )
    {
      tagSignal = MakeDoubleGauss( "dGauss_signal", mgg, *ws );
      ws->var("dGauss_signal_gauss_Ns")->setVal( (double)npoints );
    }
  else
    {
      tagSignal = MakeFullDoubleGauss( "dGauss_signal", mgg, *ws );
      ws->var("dGauss_signal_DG_Ns")->setVal( (double)npoints );
    }
  RooFitResult* sres = ws->pdf( tagSignal )->fitTo( dataSignal, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("full") );
  double gausFrac   =  ws->var("dGauss_signal_DG_frac")->getVal();
  double gausMu1    =  ws->var("dGauss_signal_DG_mu1")->getVal();
  double gausMu2    =  ws->var("dGauss_signal_DG_mu2")->getVal();
  double gausSigma1 =  ws->var("dGauss_signal_DG_sigma1")->getVal();
  double gausSigma2 =  ws->var("dGauss_signal_DG_sigma2")->getVal();

  //-------------------------------------
  //D e f i n e   S M - H i g g s   P D F
  //-------------------------------------
  npoints = dataSMH.numEntries();
  if( sameMu )
    {
      tagSMH = MakeDoubleGauss( "dGauss_SMH", mgg, *ws );
      ws->var("dGauss_SMH_gauss_Ns")->setVal( (double)npoints );
    }
  else
    {
      tagSMH = MakeFullDoubleGauss( "dGauss_SMH", mgg, *ws );
      ws->var("dGauss_SMH_DG_Ns")->setVal( (double)npoints );
    }
  RooFitResult* smhres = ws->pdf( tagSMH )->fitTo( dataSMH, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("full") );
  double gausFrac_SMH   =  ws->var("dGauss_SMH_DG_frac")->getVal();
  double gausMu1_SMH    =  ws->var("dGauss_SMH_DG_mu1")->getVal();
  double gausMu2_SMH    =  ws->var("dGauss_SMH_DG_mu2")->getVal();
  double gausSigma1_SMH =  ws->var("dGauss_SMH_DG_sigma1")->getVal();
  double gausSigma2_SMH =  ws->var("dGauss_SMH_DG_sigma2")->getVal();
  
  //------------------------------------
  // C r e a t e   b k g  s h a p e
  //------------------------------------
  
  TString tag_bkg = MakeSingleExp( "fullsb_fit_singleExp", mgg, *ws );
  RooFitResult* bres = ws->pdf( tag_bkg )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  npoints     = ws->var("fullsb_fit_singleExp_Nbkg")->getVal();
  float exp_a = ws->var("fullsb_fit_singleExp_a")->getVal();
  //creating blinded dataset
  //RooDataSet* data_toys  = GenerateToys( ws->pdf( tag_bkg ), mgg, npoints);
  //RooDataSet* data_toys  = ws->pdf( tag_bkg )->generate(mgg, npoints, kFALSE, kTRUE, "", kTRUE, kFALSE);
  RooAbsData* data_toys = ws->pdf( tag_bkg )->generateBinned( mgg, npoints, RooFit::ExpectedData() );
  
  data_toys->SetName("data_toys");
  
  //---------------------------
  //bkg_pdf for composite model
  //---------------------------
  TString tag_bkg2 = MakeSingleExpNE( "BkgModel", mgg, *ws );
  ws->var("BkgModel_sExp_a")->setVal( exp_a );
  //-----------------------------------------
  //Define double gaussian for signal and SMH
  //-----------------------------------------
  TString tag_signal2 = MakeFullDoubleGaussNE( "DG_Signal", mgg, *ws );
  TString tag_sg      = MakeSingleGaussNE( "Signal", mgg, *ws );
  TString tag_SMH2    = MakeFullDoubleGaussNE( "DG_SMH", mgg, *ws );

  std::cout << "creating s+b model" << std::endl;
  //---------------
  //S+B  m  o d e l
  //---------------
  //falling background + SM-Higgs + Signal
  //RooAddPdf* model = new RooAddPdf( "model", "model", RooArgList( *ws->pdf( tag_signal2 ), *ws->pdf( tag_SMH2 ), *ws->pdf( tag_bkg2 ) ) ) ;
  //RooAddPdf* model = new RooAddPdf( "model", "model", RooArgSet( *ws->pdf( tag_signal2 ), *ws->pdf( tag_bkg2 ) ) ) ;
  //RooAddPdf* model = new RooAddPdf( "model", "model", RooArgSet( *ws->pdf( tag_bkg2 ), *ws->pdf( tag_sg ) ) ) ;
  RooRealVar NbkgModel("NbkgModel","", npoints);
  NbkgModel.setConstant( kFALSE );
  RooRealVar NsModel("NsModel","",0.0);
  NsModel.setConstant( kFALSE );
  RooAddPdf* model = new RooAddPdf( "model", "model", RooArgSet( *ws->pdf( tag_bkg2 ), *ws->pdf( tag_sg ) ), RooArgSet(NbkgModel, NsModel) ) ;
  
  //SIGNAL
  NsModel.setRange(0, 1000.);
  ws->var("Signal_SG_mu")->setVal( gausMu1 );
  ws->var("Signal_SG_sigma")->setVal( gausSigma1 );
  //ws->var("Signal_SG_mu")->setConstant( kTRUE );
  ws->var("Signal_SG_sigma")->setConstant( kTRUE );


  //set intitial values for signal parameters
  //ws->var("DG_Signal_DG_Ns")->setVal( 15 );
  //ws->var("DG_Signal_DG_Ns")->setRange( 0, 200 );
  ws->var("DG_Signal_DG_frac")->setVal( gausFrac );
  ws->var("DG_Signal_DG_mu1")->setVal( gausMu1 );
  ws->var("DG_Signal_DG_mu2")->setVal( gausMu2 );
  ws->var("DG_Signal_DG_sigma1")->setVal( gausSigma1 );
  ws->var("DG_Signal_DG_sigma2")->setVal( gausSigma2 );
  //fix signal shape parameters
  ws->var("DG_Signal_DG_frac")->setConstant(kTRUE);
  ws->var("DG_Signal_DG_mu1")->setConstant(kTRUE);
  ws->var("DG_Signal_DG_mu2")->setConstant(kTRUE);
  ws->var("DG_Signal_DG_sigma1")->setConstant(kTRUE);
  ws->var("DG_Signal_DG_sigma2")->setConstant(kTRUE);
  //SM-Higgs
  //set intitial value for SM-Higgs parameters
  //ws->var("DG_SMH_DG_Ns")->setVal( 0.5 );
  //ws->var("DG_SMH_DG_Ns")->setMin( "", 0.0 );
  ws->var("DG_SMH_DG_frac")->setVal( gausFrac_SMH );
  ws->var("DG_SMH_DG_mu1")->setVal( gausMu1_SMH );
  ws->var("DG_SMH_DG_mu2")->setVal( gausMu2_SMH );
  ws->var("DG_SMH_DG_sigma1")->setVal( gausSigma1_SMH );
  ws->var("DG_SMH_DG_sigma2")->setVal( gausSigma2_SMH );
  //fix SM-Higgs shape parameters
  ws->var("DG_SMH_DG_frac")->setConstant(kTRUE);
  ws->var("DG_SMH_DG_mu1")->setConstant(kTRUE);
  ws->var("DG_SMH_DG_mu2")->setConstant(kTRUE);
  ws->var("DG_SMH_DG_sigma1")->setConstant(kTRUE);
  ws->var("DG_SMH_DG_sigma2")->setConstant(kTRUE);

  std::cout << "entering constraints" << std::endl;
  //--------------------------------------
  //H i g g s   C o n s t r a i n s
  //--------------------------------------
  RooRealVar HiggsYield("HiggsYield","",0.5);
  RooRealVar HiggsYieldUn("HiggsYieldUn","",0.1);
  //RooGaussian SMH_Constraint("SMH_Constraint", "SMH_Constraint", *ws->var("DG_SMH_DGF_Ns"), RooFit::RooConst(0.1), RooFit::RooConst(0.01) );
  //RooGaussian SMH_Constraint("SMH_Constraint", "SMH_Constraint", *ws->var("DG_SMH_DG_Ns"), HiggsYield, HiggsYieldUn );
  RooRealVar BkgYield ("BkgYield","", npoints);
  RooRealVar BkgYieldUn("BkgYieldUn","", 0.1*sqrt(npoints));
  //RooGaussian Bkg_Constraint("Bkg_Constraint", "Bkg_Constraint", *ws->var("BkgModel_singleExp_Nbkg"), BkgYield, BkgYieldUn );

  //---------------------------
  //SMH mass constraint
  //---------------------------
  //RooRealVar HiggsMass("HiggsMass","", ws->var("Signal_SG_mu")->getVal());
  RooRealVar HiggsMass("HiggsMass","", ws->var("Signal_SG_mu")->getVal() );
  RooRealVar HiggsMassUn("HiggsMassUn","", 0.05*ws->var("Signal_SG_mu")->getVal() );
  std::cout << "[INFO]: MC measured mass is:  " << HiggsMass.getVal() << " +/- " << HiggsMassUn.getVal() << std::endl;
  RooGaussian HiggsMass_Constraint("SMH_Constraint", "SMH_Constraint",  *ws->var("Signal_SG_mu"), HiggsMass, HiggsMassUn );
  std::cout << "pass constraints" << std::endl;
  std::cout << "pass forceSigma" << std::endl;

  //---------------------
  //F i t   t o   D a t a
  //---------------------
  //RooFitResult* sbres = model->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::ExternalConstraints(SMH_Constraint) ,RooFit::Range("low,high") );

  TH1F* h_1 = new TH1F("h_1", "h_1", 100, 0, 10);
  RooFitResult* sbres;
  //sbres = model->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("full") );
  //sbres = model->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("full") );
  sbres = model->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::ExternalConstraints(HiggsMass_Constraint) ,RooFit::Range("full") );
  
  
  for (int i = 0; i < 100; i++ )
    {
      //ws->var("Signal_SG_Ns")->setConstant( kFALSE );
      //ws->var("Signal_SG_Ns")->setVal( 0 + 0.1*i );
      //ws->var("Signal_SG_Ns")->setConstant( kTRUE );
      //sbres = model->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
      //h_1->Fill(0 + 0.1*i, sbres->minNll());
    }
  //RooFitResult* sbres = model->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
  //RooFitResult* sbres = model->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::ExternalConstraints(Bkg_Constraint),RooFit::Range("Full") );
  //Initialize fit parameter statuses
  int _status     = -1;
  int _status2    = -1;
  int _status3    = -1;
  int _covStatus  = -1;
  int _covStatus2 = -1;
  int _covStatus3 = -1;
  
  //--------------------------------------
  //Defining Negative Log Likelihood (NLL)
  //--------------------------------------
  //RooNLLVar* nll = (RooNLLVar*)model->createNLL( data, RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //Defininf RooMinimizer Object;
  //RooMinimizer m(*nll);
  //m.setStrategy(2);
  //m.setPrintLevel(3);
  //m.minimize("Minuit2", "Migrad");
  
  //m.migrad(); 
  //RooFitResult* r = m.save(); 
  //_status    = r->status();
  //_covStatus = r->covQual();
  
  //m.hesse();
  //m.minimize("Minuit2", "Hesse");
  //RooFitResult* r2 = m.save(); 
  //_status2    = r2->status();
  //_covStatus2 = r2->covQual();

  //m.setMinimizerType("Minuit2");
  //m.migrad();
  //m.hesse();
  //m.minos();
  //m.minimize("Minuit2", "Minos");
  //RooFitResult* r3 = m.save(); 
  //_status3    = r3->status();
  //_covStatus3 = r3->covQual();
  
  //------------------------------
  //E n d  M i n i m i z a t i o n
  //------------------------------
  
  sbres->SetName("sb_fit_result");
  //r->SetName("sb_fit_result_M2_migrad");
  //r2->SetName("sb_fit_result_M2_hesse");
  //r3->SetName("sb_fit_result_M2_minos");
  ws->import( *sbres );
  //ws->import( *r );
  //ws->import( *r2 );
  //ws->import( *r3 );


  //-----------------------------
  //S + B   p l o t t i n g
  //-----------------------------
  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  //data_toys->plotOn(fmgg);
  model->plotOn(fmgg, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));
  model->plotOn(fmgg, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg->SetName( "sb_fit_frame" );
  ws->import( *fmgg );
  
  //-----------------------------
  //S i g n a l   p l o t t i n g
  //-----------------------------
  RooPlot *fmgg2 = mgg.frame();
  dataSignal.plotOn(fmgg2);
  ws->pdf( tagSignal )->plotOn(fmgg2, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));
  ws->pdf( tagSignal )->plotOn(fmgg2, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg2->SetName( "signal_fit_frame" );
  ws->import( *fmgg2 );

  //---------------------------------
  //S M - H i g g s   p l o t t i n g
  //---------------------------------
  RooPlot *fmgg3 = mgg.frame();
  dataSMH.plotOn(fmgg3);
  ws->pdf( tagSMH )->plotOn(fmgg3, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));
  ws->pdf( tagSMH )->plotOn(fmgg3, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg3->SetName( "smh_fit_frame" );
  ws->import( *fmgg3 );
  

  //-------------------------------------
  //P r o f i l e d   L i k e l i h o o d 
  //-------------------------------------
  /*
    profile in Ns to obtain significance.
  */
  //RooArgSet poi   = RooArgSet( *ws->var("DG_Signal_DG_Ns") );
  //RooArgSet poi   = RooArgSet( *ws->var("Signal_SG_Ns") );
  RooArgSet poi   = RooArgSet( NsModel );
  //RooAbsReal* nll = model->createNLL( *data_toys, RooFit::Extended(kTRUE) );
  //RooAbsReal* nll = model->createNLL( data, RooFit::Extended(kTRUE) );
  RooAbsReal* nll = model->createNLL( data);
  RooFormulaVar n2ll = RooFormulaVar("n2ll", "2*@0", RooArgList(*nll) );
  RooAbsReal* p2ll = n2ll.createProfile( poi );
  //--------------------------------------
  // l i k e l i h o o d   p l o t t i n g
  //--------------------------------------
  float sVal = NsModel.getVal();
  float sErr = NsModel.getError();
  //RooPlot* fns = ws->var("DG_Signal_DG_Ns")->frame( RooFit::Range(0, 30, true) );
  //RooPlot* fns = ws->var("Signal_SG_Ns")->frame( RooFit::Range(-10, 10, true) );
  //RooPlot* fns = NsModel.frame( RooFit::Range(-10, 10, true) );
  RooPlot* fns = NsModel.frame( RooFit::Range(sVal-5.*sErr, sVal+5.*sErr, true) );
  //RooPlot* fns = ws->var("DG_Signal_DG_Ns")->frame( RooFit::Range( sVal-5.*sErr, sVal+5.*sErr, true) );
  fns->SetMinimum(0);
  fns->SetMaximum(12);
  n2ll.plotOn( fns, RooFit::ShiftToZero(), RooFit::LineColor(kBlue) );
  fns->SetMinimum(0);
  fns->SetMaximum(12);
  p2ll->plotOn( fns, RooFit::LineColor(kBlack) );
  fns->SetName("nll");
  ws->import( *fns );
  
  //ws->Write("w_sb");
  //h_1->Write("h_1");
  //ftmp->Close();
  return ws;
}

RooWorkspace* MakeDataCard( TTree* treeData, TTree* treeSignal, TTree* treeSMH, TString mggName, float SMH_Yield, std::string SMH_CF,
			    float Signal_Yield, std::string Signal_CF, TString binNumber, TString category )
{
  std::cout << "entering datacard: " << SMH_Yield << " " << Signal_Yield << std::endl;
  std::stringstream ss_smh, ss_signal;
  ss_smh << SMH_CF;
  ss_signal << Signal_CF;
  float tmp;
  std::vector<float> smh_sys, signal_sys;
  //---------------------------
  //SMH systematics into vector
  //---------------------------
  while ( ss_smh.good() )
    {
      ss_smh >> tmp;
      smh_sys.push_back( 1.0 + tmp );
      //std::cout << "tmp: " << tmp << std::endl;
      if ( ss_smh.eof() ) break;
    }
  //------------------------------
  //Signal systematics into vector
  //------------------------------
  while ( ss_signal.good() )
    {
      ss_signal >> tmp;
      signal_sys.push_back( 1.0 + tmp );
      //std::cout << "tmp: " << tmp << std::endl;
      if ( ss_signal.eof() ) break;
    }
  
  //------------------------------------------------
  // C r e a t e   s i g n a l  s h a p e from TTree
  //------------------------------------------------
  TString combinedRootFileName = "HggRazorWorkspace_bin" + binNumber + ".root";
  TFile* ftmp = new TFile( combinedRootFileName, "recreate");
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(38);
  mgg.setRange( "signal", 115, 135. );
  mgg.setRange( "high", 135, 160);
  mgg.setRange( "low", 103, 120);
  //--------------------------------
  //I m p or t i n g   D a t a
  //--------------------------------

  //Getting signal shape from signal MC
  //-----------------------
  //C r e a t e  doubleGaus
  //-----------------------
  bool sameMu = false;
  TString tagSignal, tagSMH;
  TString tag;
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*treeData) );
  RooDataSet dataSignal( "dataSignal", "", RooArgSet(mgg), RooFit::Import(*treeSignal) );
  RooDataSet dataSMH( "dataSMH", "", RooArgSet(mgg), RooFit::Import(*treeSMH) );
  //---------------------------------
  //D e f i n e   s i g n a l   P D F
  //---------------------------------
  int npoints = dataSignal.numEntries();
  if( sameMu )
    {
      tagSignal = MakeDoubleGauss( "DG_signal", mgg, *ws );
      ws->var("DG_signal_gauss_Ns")->setVal( (double)npoints );
    }
  else
    {
      tagSignal = MakeFullDoubleGauss( "DG_signal_bin"+binNumber , mgg, *ws );
      std::cout << tagSignal << std::endl;
      ws->var(tagSignal+"_Ns")->setVal( (double)npoints );
    }
  std::cout << tagSignal << std::endl;
  RooFitResult* sres = ws->pdf( tagSignal )->fitTo( dataSignal, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("Full") );
  double gausFrac    =  ws->var(tagSignal+"_frac")->getVal();
  double gausMu1     =  ws->var(tagSignal+"_mu1")->getVal();
  double gausMu2     =  ws->var(tagSignal+"_mu2")->getVal();
  double gausSigma1  =  ws->var(tagSignal+"_sigma1")->getVal();
  double gausSigma2  =  ws->var(tagSignal+"_sigma2")->getVal();

  //-------------------------------------
  //D e f i n e   S M - H i g g s   P D F
  //-------------------------------------
  npoints = dataSMH.numEntries();
  if( sameMu )
    {
      tagSMH = MakeDoubleGauss( "DG_SMH", mgg, *ws );
      ws->var("DG_SMH_gauss_Ns")->setVal( (double)npoints );
    }
  else
    {
      tagSMH = MakeFullDoubleGauss( "DG_SMH_bin"+binNumber, mgg, *ws );
      ws->var(tagSMH+"_Ns")->setVal( (double)npoints );
    }
  RooFitResult* smhres  = ws->pdf( tagSMH )->fitTo( dataSMH, RooFit::Strategy(2), RooFit::Extended( kTRUE ), RooFit::Save( kTRUE ), RooFit::Range("Full") );
  double gausFrac_SMH   =  ws->var(tagSMH+"_frac")->getVal();
  double gausMu1_SMH    =  ws->var(tagSMH+"_mu1")->getVal();
  double gausMu2_SMH    =  ws->var(tagSMH+"_mu2")->getVal();
  double gausSigma1_SMH =  ws->var(tagSMH+"_sigma1")->getVal();
  double gausSigma2_SMH =  ws->var(tagSMH+"_sigma2")->getVal();
  
  //------------------------------------
  // C r e a t e   b k g  s h a p e
  //------------------------------------
  TString tag_bkg = MakeSingleExp( "fullsb_fit_singleExp", mgg, *ws );
  //Initializing Nbkg
  npoints = data.numEntries();
  //set Nbkg Initial Value
  ws->var("fullsb_fit_singleExp_Nbkg")->setVal( npoints );
  std::cout << "entering constraints" << std::endl;
  //--------------------------------------
  //H i g g s   C o n s t r a i n s
  //--------------------------------------
  RooRealVar HiggsYield("HiggsYield","", SMH_Yield);
  //RooRealVar HiggsYieldUn("HiggsYieldUn","",SMH_YieldUn);
  //float SMH_NormUn = SMH_YieldUn/SMH_Yield;
  //RooGaussian SMH_Constraint("SMH_Constraint", "SMH_Constraint", *ws->var("DG_SMH_DGF_Ns"), RooFit::RooConst(0.1), RooFit::RooConst(0.01) );
  //RooGaussian SMH_Constraint("SMH_Constraint", "SMH_Constraint", *ws->var(tagSMH+"_Ns"), HiggsYield, HiggsYieldUn );
  std::cout << "pass constraints" << std::endl;
  std::cout << "pass forceSigma" << std::endl;

  //---------------------
  //F i t   t o   D a t a
  //---------------------
  RooFitResult* bres = ws->pdf( tag_bkg )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  float sExp_a = ws->var("fullsb_fit_singleExp_a")->getVal();
  float Nbkg   = ws->var("fullsb_fit_singleExp_Nbkg")->getVal();
  float NbkgUn = ws->var("fullsb_fit_singleExp_Nbkg")->getError();
  float BkgNormUn = 1.0 + NbkgUn/Nbkg;//input a lnN to combine
  
  //RooDataSet* data_toys = GenerateToys( ws->pdf( tag_bkg ), mgg, npoints);
  RooAbsData* data_toys = ws->pdf( tag_bkg )->generateBinned( mgg, npoints, RooFit::ExpectedData() );
  data_toys->SetName("data_bin"+binNumber);
  data.SetName("data_bin"+binNumber);
  //--------------------------------
  // m o d e l   1   p l o t t i n g
  //--------------------------------
  RooPlot *fmgg = mgg.frame();
  data_toys->plotOn(fmgg);
  ws->pdf( tag_bkg)->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  ws->pdf( tag_bkg)->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg->SetName( "fullsb_fit_frame" );
  //ws->import( *model );
  ws->import( *bres );
  ws->import( *fmgg );

  //-----------------------------
  //S i g n a l   p l o t t i n g
  //-----------------------------
  RooPlot *fmgg2 = mgg.frame();
  dataSignal.plotOn(fmgg2);
  ws->pdf( tagSignal )->plotOn(fmgg2, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));
  ws->pdf( tagSignal )->plotOn(fmgg2, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg2->SetName( "signal_fit_frame" );
  ws->import( *fmgg2 );

  //---------------------------------
  //S M - H i g g s   p l o t t i n g
  //---------------------------------
  RooPlot *fmgg3 = mgg.frame();
  dataSMH.plotOn(fmgg3);
  ws->pdf( tagSMH )->plotOn(fmgg3, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));
  ws->pdf( tagSMH )->plotOn(fmgg3, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg3->SetName( "smh_fit_frame" );
  ws->import( *fmgg3 );

  //-------------------------------------------------------
  // P r e p a r a t i o n   t o   C o m b i n e  I n p u t
  //-------------------------------------------------------
  ws->Write("w_sb");

  //--------------
  //SMH line shape
  //--------------
  RooWorkspace* combine_ws = new RooWorkspace( "combine_ws", "" );
  TString combineSMH;
  if ( category != "highres" ) combineSMH = MakeFullDoubleGaussNE( "SMH_bin"+binNumber, mgg, *combine_ws, true, true, category );//add global and cat scale uncertainty
  else combineSMH = MakeFullDoubleGaussNE( "SMH_bin"+binNumber, mgg, *combine_ws, true );//adding global scale uncertainty only
  combine_ws->var( combineSMH+"_frac")->setVal(gausFrac_SMH);
  combine_ws->var( combineSMH+"_mu1")->setVal(gausMu1_SMH);
  combine_ws->var( combineSMH+"_mu2")->setVal(gausMu2_SMH);
  combine_ws->var( combineSMH+"_sigma1")->setVal(gausSigma1_SMH);
  combine_ws->var( combineSMH+"_sigma2")->setVal(gausSigma2_SMH);
  combine_ws->var( combineSMH+"_frac")->setConstant(kTRUE);
  combine_ws->var( combineSMH+"_mu1")->setConstant(kTRUE);
  combine_ws->var( combineSMH+"_mu2")->setConstant(kTRUE);
  combine_ws->var( combineSMH+"_sigma1")->setConstant(kTRUE);
  combine_ws->var( combineSMH+"_sigma2")->setConstant(kTRUE);
  RooRealVar SMH_norm( combineSMH+"_norm" ,"", SMH_Yield);
  //SMH_norm.setConstant(kFALSE);
  combine_ws->import( SMH_norm );
  //-----------------
  //Signal line shape
  //-----------------
  TString combineSignal;
  if ( category != "highres" ) combineSignal = MakeFullDoubleGaussNE( "signal_bin"+binNumber, mgg, *combine_ws, true, true, category );//add global and cat scale uncertainty
  else combineSignal = MakeFullDoubleGaussNE( "signal_bin"+binNumber, mgg, *combine_ws, true );//add global scale uncertainty only
  combine_ws->var( combineSignal+"_frac")->setVal(gausFrac);
  combine_ws->var( combineSignal+"_mu1")->setVal(gausMu1);
  combine_ws->var( combineSignal+"_mu2")->setVal(gausMu2);
  combine_ws->var( combineSignal+"_sigma1")->setVal(gausSigma1);
  combine_ws->var( combineSignal+"_sigma2")->setVal(gausSigma2);
  combine_ws->var( combineSignal+"_frac")->setConstant(kTRUE);
  combine_ws->var( combineSignal+"_mu1")->setConstant(kTRUE);
  combine_ws->var( combineSignal+"_mu2")->setConstant(kTRUE);
  combine_ws->var( combineSignal+"_sigma1")->setConstant(kTRUE);
  combine_ws->var( combineSignal+"_sigma2")->setConstant(kTRUE);
  RooRealVar Signal_norm( combineSignal + "_norm", "", Signal_Yield );
  //Signal_norm.setConstant(kFALSE);
  combine_ws->import( Signal_norm );
  //---------
  //Bkg model
  //---------
  TString combineBkg    = MakeSingleExpNE( "Bkg_bin"+binNumber, mgg, *combine_ws );
  combine_ws->var( combineBkg + "_a" )->setVal( sExp_a );
  RooRealVar Bkg_norm(  combineBkg + "_norm", "", Nbkg );
  Bkg_norm.setConstant(kFALSE);
  combine_ws->import( Bkg_norm );
  combine_ws->import( *data_toys );
  //combine_ws->import( data );
    
  combine_ws->Write("combineWS");
  ftmp->cd();
  ftmp->Close();
  
  //std::string bNumber( binNumber );//TString to std::string
  TString dataCardName = "HggRazorCombinedCard_bin" + binNumber + ".txt";
  std::ofstream ofs( dataCardName , std::ofstream::out );
  ofs << "imax 1 number of bins\njmax 2 number of processes minus 1\nkmax * number of nuisance parameters\n";
  ofs << "----------------------------------------------------------------------------------------\n";
  ofs << "shapes Bkg\t\tbin"      << binNumber << "\t" << combinedRootFileName << " combineWS:" << combineBkg << "\n";
  ofs << "shapes SMH\t\tbin"      << binNumber << "\t" << combinedRootFileName << " combineWS:" << combineSMH << "\n";
  ofs << "shapes signal\t\tbin"   << binNumber << "\t" << combinedRootFileName << " combineWS:" << combineSignal << "\n";
  ofs << "shapes data_obs\t\tbin" << binNumber << "\t" << combinedRootFileName << " combineWS:" << "data_bin" << binNumber << "\n";
  ofs << "----------------------------------------------------------------------------------------\n";
  ofs << "bin\t\tbin" << binNumber << "\n";
  ofs << "observation\t-1.0\n";
  ofs << "----------------------------------------------------------------------------------------\n";
  ofs << "bin\t\t\t\t\t\tbin" << binNumber << "\t\tbin" << binNumber << "\t\tbin" << binNumber << "\n";
  ofs << "process\t\t\t\t\t\tsignal\t\tSMH\t\tBkg\n";
  ofs << "process\t\t\t\t\t\t0\t\t1\t\t2\n";
  ofs << "rate\t\t\t\t\t\t1\t\t1\t\t1\n";
  ofs << "----------------------------------------------------------------------------------------\n";
  ofs << "CMS_Lumi\t\t\tlnN\t\t1.027\t\t1.027\t\t-\n";
  ofs << "Photon_Trigger\t\t\tlnN\t\t1.05\t\t1.05\t\t-\n";
  
  ofs << "ScaleNorm\t\t\tlnN\t\t-\t\t0.931/1.065\t\t-\n";
  ofs << "PdfNorm\t\t\t\tlnN\t\t-\t\t0.948/1.062\t\t-\n";
  int totalSys = smh_sys.size();
  int ctr = 0;
  for( int isys = 0; isys < totalSys; isys++ )
    {
      if ( isys == 0 )
	{
	  ofs << "SMH_JES\t\t\t\tlnN\t\t-\t\t" << smh_sys.at(isys+1) << "/" << smh_sys.at(isys) << "\t\t-\n";
	}
      else if ( isys == 2 )
	{
	  ofs << "SMH_facScale\t\t\tlnN\t\t-\t\t" << smh_sys.at(isys+1) << "/" << smh_sys.at(isys) << "\t\t-\n";
	}
      else if ( isys == 4 )
	{
	  ofs << "SMH_renScale\t\t\tlnN\t\t-\t\t" << smh_sys.at(isys+1) << "/" << smh_sys.at(isys) << "\t\t-\n";
	}
      else if ( isys == 6 )
	{
	  ofs << "SMH_facRenScale\t\t\tlnN\t\t-\t\t" << smh_sys.at(isys+1) << "/" << smh_sys.at(isys) << "\t\t-\n";
	}
      else if ( isys > 7 )
	{
	  ofs << "SMH_pdf" << ctr << "\t\t\tlnN\t\t-\t\t" << smh_sys.at(isys) << "\t\t-\n";
	  ctr++;
	}
    }
  
  //ofs << "mu_Global\t\t\tparam\t\t 0 1.25\n";
  //if ( category != "highres" ) ofs << category << "_mu_Global\t\t\tparam\t\t 0 1.25\n";
  if ( category == "hzbb" )
    {
      ofs << "SMH_btag\t\t\tlnN\t\t-\t\t" << "0.961/1.04" "\t\t-\n";
      ofs << "SMH_misstag\t\t\tlnN\t\t-\t\t" << "0.992/1.008" << "\t\t-\n";
    }
  //ofs << "SMH_renScale\t\t\tlnN\t\t-\t\t" << SMH_renScale << "\t\t-\n";
  //ofs << "SMH_facRenScale\t\t\tlnN\t\t-\t\t" << SMH_facRenScale << "\t\t-\n";
  //ofs << "BkgNorm_bin" << binNumber << "\t\t\tlnN\t\t-\t\t-\t\t" << BkgNormUn << std::endl;
  ofs.close();
  return ws;
};


RooWorkspace* MakeSignalBkgFit( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName )
{

  //-------------------------------------
  // C r e a t e   s i g n a l  s h a p e
  //------------------------------------- 
  //Getting signal shape from MC
  TFile* fsignal = new TFile("signalModels/RSGraviton_750GeV.root", "read");
  RooWorkspace* wsignal = (RooWorkspace*)fsignal->Get("w_sFit");
  RooAbsPdf* signalPdf = wsignal->pdf("dGauss_signal_doublegauss");
  double gausFrac   =  wsignal->var("dGauss_signal_frac")->getVal();
  double gausMu     =  wsignal->var("dGauss_signal_gauss_mu")->getVal();
  double gausSigma1 =  wsignal->var("dGauss_signal_gauss_sigma1")->getVal();
  double gausSigma2 =  wsignal->var("dGauss_signal_gauss_sigma2")->getVal();

  delete fsignal;
  //fsignal->Close();
  
  TFile* ftmp = new TFile("tmp_output_OurID.root", "recreate");
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  /*  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  mgg.setRange("signal", 103, 160);*/
  
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 220, 1000, "GeV" );
  mgg.setBins(39);
  mgg.setRange( "signal", 650, 850. );
  mgg.setRange("high", 850, 1000);
  mgg.setRange("low", 220, 650);
  //--------------------------------
  //I m p or t i n g   D a t a
  //--------------------------------
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  
  //------------------------------------
  // C r e a t e   b k g  s h a p e
  //------------------------------------
  //TString tag_bkg = MakeDoubleExpN1N2( "fullsb_fit_bkg", mgg, *ws );
  TString tag_bkg = MakeDoubleExpN1N2( "fullsb_fit_doubleExp", mgg, *ws );
  //TString tag_bkg = MakeSingleExp( "fullsb_fit_singleExp", mgg, *ws );
  TString tag_bkg2 = MakeSingleExp( "fullsb_fit_singleExp2", mgg, *ws );
  
  
  std::cout << "reading model from file" << std::endl;
  //Define arbitrarily a double gaussian
  //TString tag_signal = MakeDoubleGauss( "doubleGaussSB1", mgg, *ws );
  TString tag_signal2 = MakeDoubleGauss( "doubleGaussSB", mgg, *ws );

  std::cout << "creating s+b model" << std::endl;
  //---------------
  //S+B  m  o d e l
  //---------------
  RooAddPdf* model = new RooAddPdf( "model", "model", RooArgSet( *ws->pdf( tag_signal2 ), *ws->pdf( tag_bkg ) ) ) ;
  //RooAddPdf* model2 = new RooAddPdf( "model2", "model2", RooArgSet( *ws->pdf( tag_signal2 ), *ws->pdf( tag_bkg2 ) ) ) ;
  //Fixing gaussian parameter for 2nd model
  ws->var("doubleGaussSB_gauss_Ns")->setVal( 0 );
  ws->var("doubleGaussSB_frac")->setVal( gausFrac );
  ws->var("doubleGaussSB_gauss_mu")->setVal( gausMu );
  ws->var("doubleGaussSB_gauss_sigma1")->setVal( gausSigma1 );
  ws->var("doubleGaussSB_gauss_sigma2")->setVal( gausSigma2 );
  ws->var("doubleGaussSB_gauss_Ns")->setConstant(kTRUE);
  ws->var("doubleGaussSB_frac")->setConstant(kTRUE);
  ws->var("doubleGaussSB_gauss_mu")->setConstant(kTRUE);
  ws->var("doubleGaussSB_gauss_sigma1")->setConstant(kTRUE);
  ws->var("doubleGaussSB_gauss_sigma2")->setConstant(kTRUE);

  std::cout << "entering constraints" << std::endl;
  //--------------------------------------
  //H i g g s   C o n s t r a i n s
  //--------------------------------------
  /*
  RooRealVar HiggsMass("HiggsMass","",128.7);
  RooRealVar HiggsMassError("HiggsMassError","",0.4);
  RooGaussian HiggsMassConstraint("HiggsMassConstraint","", *ws->var("fullsb_fit_signal_gauss_mu"),HiggsMass,HiggsMassError);
  std::cout << "pass constraints" << std::endl;
  std::cout << "entering forceSigma" << std::endl;
  if( forceSigma != -1 ) {
    ws->var("fullsb_fit_signal_gauss_sigma1")->setVal( forceSigma );
    ws->var("fullsb_fit_signal_gauss_sigma1")->setConstant(true);
    ws->var("fullsb_fit_signal_frac")->setVal( 1.0 );
    ws->var("fullsb_fit_signal_frac")->setConstant(true);
  }
  */
  std::cout << "pass forceSigma" << std::endl;

  //---------------------
  //F i t   t o   D a t a
  //---------------------
  
  model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  RooFitResult* bres = model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
  
  //model2->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //RooFitResult* bres2 = model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
 

  //-------------------------------------
  //P r o f i l e d   L i k e l i h o o d 
  //-------------------------------------
  /*
    profile in Ns to obtain significance.
  */
  RooArgSet poi   = RooArgSet( *ws->var("doubleGaussSB_gauss_Ns") );
  RooAbsReal* nll = model->createNLL(data);
  RooFormulaVar n2ll = RooFormulaVar("n2ll", "2*@0", RooArgList(*nll) );
  RooAbsReal* p2ll = n2ll.createProfile( poi );
    
  
  /*
  model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full"), RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)) );
  RooFitResult* bres = model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full"), RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)) );
  */
  //--------------------------------
  // m o d e l   1   p l o t t i n g
  //--------------------------------
  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  model->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  model->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg->SetName( "fullsb_fit_frame" );
  ws->import( *model );
  ws->import( *bres );
  ws->import( *fmgg );
  
  //--------------------------------
  // m o d e l   2   p l o t t i n g
  //--------------------------------
  /*
  RooPlot *fmgg2 = mgg.frame();
  data.plotOn(fmgg2);
  model2->plotOn(fmgg2,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  model2->plotOn(fmgg2,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg2->SetName( "fullsb_fit_frame2" );
  ws->import( *model2 );
  ws->import( *bres2 );
  ws->import( *fmgg2 );*/

  //--------------------------------------
  // l i k e l i h o o d   p l o t t i n g
  //--------------------------------------
  RooPlot* fns = ws->var("doubleGaussSB_gauss_Ns")->frame( RooFit::Range(0, 20, true) );
  fns->SetMinimum(0);
  fns->SetMaximum(6);
  n2ll.plotOn( fns, RooFit::ShiftToZero(), RooFit::LineColor(kBlue) );
  p2ll->plotOn( fns, RooFit::LineColor(kBlack) );
  fns->SetName("nll_trick");
  ws->import( *fns );
  ws->Write("w_sb");
  ftmp->Close();
  return ws;
};

RooWorkspace* MakeSideBandFitAIC( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName, double& AIC, TString ffName = "doubleExp" )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(57);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 135, 160);
  mgg.setRange("Full", 103, 160);

  TString tag;
  if ( ffName == "doubleExp" )
    {
      tag = MakeDoubleExpN1N2( "sideband_fit_doubleExp", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( ffName == "singleExp" )
    {
      tag = MakeSingleExp( "sideband_fit_singleExp", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( ffName == "modExp" )
    {
      tag = MakeModExp( "sideband_fit_modExp", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( ffName == "singlePow" )
    {
      tag = MakeSinglePow( "sideband_fit_singlePow", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( ffName == "doublePow" )
    {
      tag = MakeDoublePow( "sideband_fit_doublePow", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( ffName == "poly2" )
    {
      tag = MakePoly2( "sideband_fit_poly2", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( ffName == "poly3" )
    {
      tag = MakePoly3( "sideband_fit_poly3", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  
  std::cout << "[INFO]: ENTERING FIT" << std::endl;
  //Sideband Fit
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  //ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  //RooFitResult* bres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  RooFitResult* bres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
  //ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //RooFitResult* bres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
    
  bres->SetName( tag + "_b_fitres" );
  std::cout << "[INFO]: PASS FIT" << std::endl;
  //---------------------
  //g e t t i n g   n l l 
  //---------------------
  double minNll = bres->minNll();
  //RooAbsReal* nll = ws->pdf( tag )->createNLL(data,  RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  RooAbsReal* nll = ws->pdf( tag )->createNLL(data, RooFit::Extended(kTRUE), RooFit::Range("Full") );
  std::cout << "nll_nll->" << nll->getVal() << std::endl;
  std::cout << "minNll->" << minNll << std::endl;
  RooArgSet* floatPars = ws->pdf( tag )->getParameters(data);
  double K = floatPars->getSize() - 1.;
  std::cout << "K -> " << K << std::endl;
  double n = data.sumEntries(" (mgg>103 && mgg<120) || (mgg>131 && mgg<160)");
  std::cout << "n -> " << n << std::endl;
  AIC = 2*minNll + 2*K + 2*K*(K+1)/(n-K-1);
  std::cout << "AIC: " << AIC << std::endl;
  /*
  RooPlot* fns = ws->var("sideband_fit_doubleExp_a1")->frame( );
  nll->plotOn( fns, RooFit::LineColor(kBlue) );
  fns->SetName("nll_trick");
  ws->import( *fns );
  */
  ws->import( *bres );
  
  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  ws->pdf( tag )->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  ws->pdf( tag )->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  
  fmgg->SetName( tag + "_frame" );
  ws->import( *fmgg );
  
  RooPlot* pdfFrame = mgg.frame();
  ws->pdf( tag )->plotOn( pdfFrame, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("Full") );
  pdfFrame->SetName( tag+"_pdfframe" );
  ws->import( *pdfFrame );
  ws->import( mgg );

  return ws;
};

RooWorkspace* DoBiasTest( TTree* tree, TString mggName, TString f1, TString f2, int ntoys, int npoints )
{
  RooRandom::randomGenerator()->SetSeed( 0 );
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg( mggName,"m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 135, 160);
  mgg.setRange("sig", 122.08, 128.92);
  mgg.setRange("Full", 103, 160);
  /*
  TString tag1, tag2;
  if ( f1 == "doubleExp" )
    {
      tag1 = MakeDoubleExpN1N2( "sideband_fit_doubleExp", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( f1 == "singleExp" )
    {
      tag1 = MakeSingleExp( "sideband_fit_singleExp", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( f1 == "modExp" )
    {
      tag1 = MakeModExp( "sideband_fit_modExp", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( f1 == "singlePow" )
    {
      tag1 = MakeSinglePow( "sideband_fit_singlePow", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( f1 == "doublePow" )
    {
      tag1 = MakeDoublePow( "sideband_fit_doublePow", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( f1 == "poly2" )
    {
      tag1 = MakePoly2( "sideband_fit_poly2", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( f1 == "poly3" )
    {
      tag1 = MakePoly3( "sideband_fit_poly3", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  //------------------
  //f2
  //------------------
  if ( f2 == "doubleExp" )
    {
      tag2 = MakeDoubleExpN1N2( "sideband_fit_doubleExp", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( f2 == "singleExp" )
    {
      tag2 = MakeSingleExp( "sideband_fit_singleExp", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( f2 == "modExp" )
    {
      tag2 = MakeModExp( "sideband_fit_modExp", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( f2 == "singlePow" )
    {
      tag2 = MakeSinglePow( "sideband_fit_singlePow", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( f2 == "doublePow" )
    {
      tag2 = MakeDoublePow( "sideband_fit_doublePow", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( f2 == "poly2" )
    {
      tag2 = MakePoly2( "sideband_fit_poly2", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( f2 == "poly3" )
    {
      tag2 = MakePoly3( "sideband_fit_poly3", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  */
  TString tag1, tag2, tag2p;
  if ( f1 == "doubleExp" )
    {
      tag1 = MakeDoubleExpN1N2( "doubleExp_1", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( f1 == "singleExp" )
    {
      tag1 = MakeSingleExp( "singleExp_1", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( f1 == "modExp" )
    {
      tag1 = MakeModExp( "modExp_1", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( f1 == "singlePow" )
    {
      tag1 = MakeSinglePow( "singlePow_1", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( f1 == "doublePow" )
    {
      tag1 = MakeDoublePow( "doublePow_1", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( f1 == "poly2" )
    {
      tag1 = MakePoly2( "poly2_1", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( f1 == "poly3" )
    {
      tag1 = MakePoly3( "poly3_1", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  //------------------
  //f2
  //------------------
  if ( f2 == "doubleExp" )
    {
      tag2  = MakeDoubleExpN1N2( f2 + "_2", mgg, *ws );
      tag2p = MakeDoubleExpN1N2( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( f2 == "singleExp" )
    {
      tag2  = MakeSingleExp( f2 + "_2", mgg, *ws );
      tag2p = MakeSingleExp( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( f2 == "modExp" )
    {
      tag2  = MakeModExp( f2 + "_2", mgg, *ws );
      tag2p = MakeModExp( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( f2 == "singlePow" )
    {
      tag2  = MakeSinglePow( f2 + "_2", mgg, *ws );
      tag2p = MakeSinglePow( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( f2 == "doublePow" )
    {
      tag2  = MakeDoublePow( f2 + "_2", mgg, *ws );
      tag2p = MakeDoublePow( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( f2 == "poly2" )
    {
      tag2  = MakePoly2( f2 + "_2", mgg, *ws );
      tag2p = MakePoly2( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( f2 == "poly3" )
    {
      tag2  = MakePoly3( f2 + "_2", mgg, *ws );
      tag2p = MakePoly3( f2 + "_prime", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
    
  
  
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  npoints = data.numEntries();
  //Sideband Fit (not working with poly2 and poly3)
  //ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  //RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  //FullFit
  RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  //RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );

  RooPlot* f1_frame = mgg.frame();
  f1_frame->SetName("f1_frame");
  data.plotOn( f1_frame );
  ws->pdf( tag1 )->plotOn( f1_frame, RooFit::LineColor( kRed), RooFit::Range("low,high"), RooFit::NormRange("low,high") );
  ws->import( *f1_frame );
  RooAbsReal* f1Integral = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  std::cout << "f1 Int: " << f1Integral->getVal() << std::endl;
  
  RooDataSet* data_toys;
  RooFitResult* bres_toys;
  double n; 
  RooAbsReal* fIntegral;
  RooAbsReal* fIntegral2;
  RooRealVar bias("bias", "bias", -5.0, 5.0, "");
  RooDataSet data_bias( "data_bias", "bias data", bias);
  bias.setBins(100);
  
  
  //-----------------------
  //do a fit to bias data
  //-----------------------
  ws->pdf( tag2p )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  RooPlot* pFrame = mgg.frame();
  pFrame->SetName("pFitFrame");
  data.plotOn( pFrame );
  ws->pdf( tag2p )->plotOn( pFrame, RooFit::LineColor(kBlue), RooFit::Range("low,high"), RooFit::NormRange("low,high") );
  ws->pdf( tag2p )->plotOn( pFrame, RooFit::LineStyle(kDashed), RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("low,high") );
  ws->import( *pFrame );
  
  double dE_N1, dE_N2, dE_a1, dE_a2;//doubleExp
  double sE_N, sE_a;//singleExp
  double mE_N, mE_a, mE_m;//modExp
  double sP_N, sP_a;//singlePow
  double dP_N, dP_f, dP_a1, dP_a2;//doubleExp
  double pC, p0, p1, p2, pN;//poly2,pol3;
  if ( f2 == "doubleExp" )
    {
      dE_N1  = ws->var( f2 + "_prime_Nbkg1" )->getVal();
      dE_N2  = ws->var( f2 + "_prime_Nbkg2" )->getVal();
      dE_a1  = ws->var( f2 + "_prime_a1" )->getVal();
      dE_a2  = ws->var( f2 + "_prime_a2" )->getVal();
    }
  else if ( f2 == "singleExp" )
    {
      sE_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      sE_a  = ws->var( f2 + "_prime_a" )->getVal();
    }
  else if ( f2 == "modExp" )
    {
      mE_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      mE_a  = ws->var( f2 + "_prime_a" )->getVal();
      mE_m  = ws->var( f2 + "_prime_m" )->getVal();
    }
  else if ( f2 == "singlePow" )
    {
      sP_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      sP_a  = ws->var( f2 + "_prime_a" )->getVal();
    }
  else if ( f2 == "doublePow" )
    {
      dP_N   = ws->var( f2 + "_prime_Nbkg" )->getVal();
      dP_f   = ws->var( f2 + "_prime_f" )->getVal();
      dP_a1  = ws->var( f2 + "_prime_a1" )->getVal();
      dP_a2  = ws->var( f2 + "_prime_a2" )->getVal();
    }
  else if ( f2 == "poly2" )
    {
      pN = ws->var( f2 + "_prime_Nbkg" )->getVal();
      pC = ws->var( f2 + "_prime_pC" )->getVal();
      p0 = ws->var( f2 + "_prime_p0" )->getVal();
      p1 = ws->var( f2 + "_prime_p1" )->getVal();
    }
  else if ( f2 == "poly3" )
    {
      pN = ws->var( f2 + "_prime_Nbkg" )->getVal();
      pC = ws->var( f2 + "_prime_pC" )->getVal();
      p0 = ws->var( f2 + "_prime_p0" )->getVal();
      p1 = ws->var( f2 + "_prime_p1" )->getVal();
      p2 = ws->var( f2 + "_prime_p2" )->getVal();
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  for ( int i = 0; i < ntoys; i++ )
    {
      //-----------------------------
      //G e n e r a t i n g   t o y s
      //-----------------------------
      data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npoints );
      //ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
      //RooFitResult* bres_toys = ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
      if ( f2 == "doubleExp" )
      	{
	  ws->var( f2 + "_2_Nbkg1")->setVal( dE_N1 );
	  ws->var( f2 + "_2_Nbkg2")->setVal( dE_N2 );
	  ws->var( f2 + "_2_a1")->setVal( dE_a1 );
	  ws->var( f2 + "_2_a2")->setVal( dE_a2 );
	}
      else if ( f2 == "singleExp" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( sE_N );
	  ws->var( f2 + "_2_a" )->setVal( sE_a );
	}
      else if ( f2 == "modExp" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( mE_N );
	  ws->var( f2 + "_2_a" )->setVal( mE_a );
	  ws->var( f2 + "_2_m" )->setVal( mE_m );
	}
      else if ( f2 == "singlePow" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( sP_N );
	  ws->var( f2 + "_2_a" )->setVal( sP_a );
	}
      else if ( f2 == "doublePow" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( dP_N );
	  ws->var( f2 + "_2_f" )->setVal( dP_f );
	  ws->var( f2 + "_2_a1" )->setVal( dP_a1 );
	  ws->var( f2 + "_2_a2" )->setVal( dP_a2 );
	}
      else if ( f2 == "poly2" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( pN );
	  ws->var( f2 + "_2_pC" )->setVal( pC );
	  ws->var( f2 + "_2_p0" )->setVal( p0 );
	  ws->var( f2 + "_2_p1" )->setVal( p1 );
	}
      else if ( f2 == "poly3" )
	{
	  ws->var( f2 + "_2_Nbkg" )->setVal( pN );
	  ws->var( f2 + "_2_pC" )->setVal( pC );
	  ws->var( f2 + "_2_p0" )->setVal( p0 );
	  ws->var( f2 + "_2_p1" )->setVal( p1 );
	  ws->var( f2 + "_2_p2" )->setVal( p2 );
	}
      else
	{
	  ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
	}
      
      
      bres_toys = ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
      bres_toys->SetName("bres_toys");
      
      fIntegral2 = ws->pdf( tag2 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
      std::cout << "test Int2: " << fIntegral2->getVal() << std::endl;
      //std::cout << "signal events fit: " <<  (fIntegral2->getVal())*pow(ws->var("sideband_fit_singleExpse_Nbkg")->getVal(),2) << std::endl;
      bias =  (fIntegral2->getVal() - f1Integral->getVal())/f1Integral->getVal();
      data_bias.add(RooArgSet(bias));
    }
  RooPlot* f_mgg = mgg.frame();
  f_mgg->SetName("toys_plot");
  data_toys->plotOn( f_mgg );
  ws->pdf( tag2 )->plotOn( f_mgg, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full"));

  RooPlot* f_bias = bias.frame();
  f_bias->SetName("bias_plot");
  data_bias.plotOn( f_bias );
  ws->import( mgg );
  ws->import( bias );
  ws->import( data_bias );
  ws->import( *f_mgg );
  ws->import( *bres_toys );
  ws->import( *f_bias );

  return ws;
};

RooWorkspace* DoBiasTestSignal( TTree* tree, TString mggName, TString f1, TString f2, int ntoys, double frac, TString outName )
{
  RooRandom::randomGenerator()->SetSeed( 0 );
  RooWorkspace* ws = new RooWorkspace( "ws", "" );

  //Getting signal shape
  std::cout << "[INFO]: getting signal model" << std::endl;
  TFile* fsignal = new TFile("data/SM_MC_GluGluH_Fit.root", "read");
  RooWorkspace* wsignal = (RooWorkspace*)fsignal->Get("w_sFit");
  //Getting signal pdf and parementer except Ns
  RooAbsPdf* signalPdf = wsignal->pdf("dGauss_signal_doublegauss");
  double gausFrac   =  wsignal->var("dGauss_signal_frac")->getVal();
  double gausMu     =  wsignal->var("dGauss_signal_gauss_mu")->getVal();
  double gausSigma1 =  wsignal->var("dGauss_signal_gauss_sigma1")->getVal();
  double gausSigma2 =  wsignal->var("dGauss_signal_gauss_sigma2")->getVal();

  delete fsignal;

  //---------------------------------
  //Define output Tree for bias study
  //---------------------------------
  TTree* outTree = new TTree("biasTree", "tree containing bias tests");
  //define variables
  double _bias, _statUn, _biasNorm;
  double _Ns_hat, _Ns_sigma, _Ns, _alpha_hat, _alpha_sigma, _alpha, _Nbkg_hat, _Nbkg_sigma, _Nbkg;
  int _status, _status2, _status3, _covStatus, _covStatus2, _covStatus3;
  //define branches
  outTree->Branch("bias", &_bias, "bias/D");
  outTree->Branch("statUn", &_statUn, "statUn/D");
  outTree->Branch("biasNorm", &_biasNorm, "biasNorm/D");
  outTree->Branch("alpha_hat", &_alpha_hat, "alpha_hat/D");
  outTree->Branch("alpha_sigma", &_alpha_sigma, "alpha_sigma/D");
  outTree->Branch("alpha", &_alpha, "alpha/D");
  outTree->Branch("Ns_hat", &_Ns_hat, "Ns_hat/D");
  outTree->Branch("Ns_sigma", &_Ns_sigma, "Ns_sigma/D");
  outTree->Branch("Ns", &_Ns, "Ns/D");
  outTree->Branch("Nbkg_hat", &_Nbkg_hat, "Nbkg_hat/D");
  outTree->Branch("Nbkg_sigma", &_Nbkg_sigma, "Nbkg_sigma/D");
  outTree->Branch("Nbkg", &_Nbkg, "Nbkg/D");
  outTree->Branch("status", &_status, "status/I");
  outTree->Branch("covStatus", &_covStatus, "covStatus/I");
  outTree->Branch("status2", &_status2, "status2/I");
  outTree->Branch("covStatus2", &_covStatus2, "covStatus2/I");
  outTree->Branch("status3", &_status3, "status3/I");
  outTree->Branch("covStatus3", &_covStatus3, "covStatus3/I");

  //-------------------
  //Creating observable
  //-------------------
  RooRealVar mgg( mggName,"m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(38);
  mgg.setRange("low", 103, 120. );//low sideband
  mgg.setRange("high", 135, 160.);//high sideband
  mgg.setRange("sig", 122.08, 128.92);//highRes signal region
  mgg.setRange("Full", 103.0, 160.0);//Full range for this analysis

  //-----------------------------------------------------------------
  //Generate signal toys to check the signal shape in the output file
  //-----------------------------------------------------------------
  RooDataSet* signal_toys = GenerateToys( signalPdf, mgg, 1e4 );
  signal_toys->SetName("signal_toys_new");
  RooPlot* s_mgg = mgg.frame();
  s_mgg->SetName("signal_toys");
  signal_toys->plotOn( s_mgg );
  ws->import( *signal_toys );
  ws->import( *s_mgg );

  //---------------------
  //Define parent pdf #1
  //---------------------
  TString tag1, tag2, tag2p;
  if ( f1 == "doubleExp" )
    {
      //tag1 = MakeDoubleExpN1N2( f1 + "_1", mgg, *ws );
      tag1 = MakeDoubleExp( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "singleExp" )
    {
      tag1 = MakeSingleExp( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "modExp" )
    {
      tag1 = MakeModExp( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "singlePow" )
    {
      tag1 = MakeSinglePow( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "doublePow" )
    {
      tag1 = MakeDoublePow( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "poly2" )
    {
      tag1 = MakePoly2( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "poly3" )
    {
      tag1 = MakePoly3( f1 + "_1", mgg, *ws );
    }
  else if ( f1 == "poly4" )
    {
      tag1 = MakePoly4( f1 + "_1", mgg, *ws );
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  std::cout << "[INFO]: f1 is a " << f1 << std::endl;
  
  //------------------
  //Define test pdf #2
  //------------------
  if ( f2 == "doubleExp" )
    {
      //tag2  = MakeDoubleExpN1N2( f2 + "_2", mgg, *ws );
      //tag2p = MakeDoubleExpN1N2( f2 + "_prime", mgg, *ws );
      tag2  = MakeDoubleExpNE( f2 + "_2", mgg, *ws );
      tag2p = MakeDoubleExp( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "singleExp" )
    {
      tag2  = MakeSingleExpNE( f2 + "_2", mgg, *ws );
      tag2p = MakeSingleExp( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "modExp" )
    {
      tag2  = MakeModExpNE( f2 + "_2", mgg, *ws );
      tag2p = MakeModExp( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "singlePow" )
    {
      tag2  = MakeSinglePowNE( f2 + "_2", mgg, *ws );
      tag2p = MakeSinglePow( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "doublePow" )
    {
      tag2  = MakeDoublePowNE( f2 + "_2", mgg, *ws );
      tag2p = MakeDoublePow( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "poly2" )
    {
      tag2  = MakePoly2NE( f2 + "_2", mgg, *ws );
      tag2p = MakePoly2( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "poly3" )
    {
      tag2  = MakePoly3NE( f2 + "_2", mgg, *ws );
      tag2p = MakePoly3( f2 + "_prime", mgg, *ws );
    }
  else if ( f2 == "poly4" )
    {
      tag2  = MakePoly4NE( f2 + "_2", mgg, *ws );
      tag2p = MakePoly4( f2 + "_prime", mgg, *ws );
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  std::cout << "[INFO]: f2 is a " << f2 << std::endl;

  //-------------------------------------
  //Define data set, importing from TTree
  //-------------------------------------
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  //getting total n_entries and n_sideband
  //double n_sideband = data.sumEntries("(mgg > 103. && mgg < 120.) || (mgg > 135. && mgg < 160.)");
  TString sbCut = Form("(%s>%.2f && %s<%.2f) || (%s>%.2f && %s<%.2f)", mggName.Data(), 103., mggName.Data(), 120., mggName.Data(), 135., mggName.Data(), 160.);
  double n_sideband = data.sumEntries( sbCut );

  //----------------------------------------------
  //Sideband Fit to data set using parent function
  //----------------------------------------------
  RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  //_alpha = ws->var( f1 + "_1_a" )->getVal();
  //FullFit
  //RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );

  //---------------------------------------
  //P l o t t i n g   s i de b a n d  f i t 
  //---------------------------------------
  RooPlot* bFrame = mgg.frame();
  bFrame->SetName("bFitFrame");
  bFrame->SetTitle("");
  data.plotOn( bFrame );
  ws->pdf( tag1 )->plotOn( bFrame, RooFit::LineColor(kBlue), RooFit::Range("low,high"), RooFit::NormRange("low,high"));
  ws->import( *bFrame );
  //saving firResult
  bres->SetName( "fit_result_f1" );
  ws->import( *bres );

  //-----------------------------------
  //g e t t i n g   i n t e g t r a l s
  //-----------------------------------
  RooAbsReal* f1Integral = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  double f1Int = f1Integral->getVal();

  RooAbsReal* f1Integral_sb = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("low,high") );
  double f1Int_sb = f1Integral_sb->getVal();
  int npoints = (int)n_sideband/f1Int_sb;//re-scaling sideband to total bkg events (N_sideband/NORMALIZE_INTEGRAL_SIDEBAND)
  npoints = 4*npoints;
  //npoints = 350;//only use this to set the number of toys bkg;

  std::cout << "npoints: " << npoints << std::endl;
  //-------------------------------
  //S i g n a l   +   B k g   P d f
  //-------------------------------
  //bkg non-extended, see tag2 above
  //TString bkgTag = MakeSinglePowNE("bkg_model", mgg, *ws );
  //signal non-extended
  TString gaussTag = MakeDoubleGaussNE("signal", mgg, *ws );
  RooRealVar Ns( "sbModel_Ns", "N_{s}", 0, "" );
  Ns.setConstant(kFALSE);
  //Ns.setRange(-5,1e4);
  RooRealVar Nbkg( "sbModel_Nbkg", "N_{bkg}", 0, "" );
  Nbkg.setVal(npoints);
  Nbkg.setConstant(kFALSE);
  RooAddPdf* sbModel = new RooAddPdf( "sbModel", "sbModel", RooArgList( *ws->pdf(tag2), *ws->pdf(gaussTag) ), RooArgList( Nbkg, Ns ) );
  ws->import( *sbModel );
  
  RooDataSet*   data_toys;
  RooFitResult* bres_toys;
  double n; 

  //---------------------
  //RooFit bias variables
  //---------------------
  RooAbsReal* fIntegral;
  RooAbsReal* fIntegral2;
  RooRealVar bias("bias", "bias", -50.0, 50.0, "");
  RooDataSet data_bias( "data_bias", "bias data", bias);

  RooRealVar NsignalError("NsignalError", "NsignalError", -5.0, 5.0, "");
  RooDataSet data_Nse( "data_Nse", "data_Nse", NsignalError);
  
  bias.setBins(100);
  NsignalError.setBins(100);
  //----------------------------------
  //Setting initial parameter for Nbkg
  //----------------------------------
  if ( f2 == "doubleExp" )
    {
      //ws->var( f2 + "_prime_Nbkg1" )->setVal(npoints);
      //ws->var( f2 + "_prime_Nbkg2" )->setVal(.0);
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "singleExp" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "modExp" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "singlePow" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "doublePow" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "poly2" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal(npoints);
    }
  else if ( f2 == "poly3" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal( npoints );
    }
  else if ( f2 == "poly4" )
    {
      ws->var( f2 + "_prime_Nbkg" )->setVal( npoints );
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }

  //----------------------------------------
  //fit to data to obtain initial parameters
  //----------------------------------------
  //Generate dataset 
  data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npoints );
  data_toys->SetName("dataToysBkgOnly");
  RooFitResult* bres2p = ws->pdf( tag2p )->fitTo( *data_toys, RooFit::Save(kTRUE), RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  //RooFitResult* bres2p = ws->pdf( tag2p )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  bres2p->SetName("fit_result_f2_prime");
  ws->import( *bres2p );
  
  RooPlot* pFrame = mgg.frame();
  pFrame->SetName("primeFitFrame");
  data.plotOn( pFrame );
  //---------------------------------
  //Plotting sideband fit f1, f2prime(copy of f2)
  //---------------------------------
  ws->pdf( tag2p )->plotOn( pFrame, RooFit::LineColor(kRed), RooFit::Range("low,high"), RooFit::NormRange("low,high") );
  ws->pdf( tag1 )->plotOn( pFrame, RooFit::LineColor(kBlue), RooFit::Range("low,high"), RooFit::NormRange("low,high") );
  ws->import( *pFrame );

  //------------------------------------------------------------------------------
  //Define and obtain initial pdf parameters for f2, using sideband fit parameters
  //------------------------------------------------------------------------------
  double dE_N1, dE_N2, dE_a1, dE_a2;//doubleExp
  double sE_N, sE_a;//singleExp
  double mE_N, mE_a, mE_m;//modExp
  double sP_N, sP_a;//singlePow
  double dP_N, dP_f, dP_a1, dP_a2;//doubleExp
  double pC, p0, p1, p2, p3, pN;//poly2,pol3;
  if ( f2 == "doubleExp" )
    {
      //dE_N1  = ws->var( f2 + "_prime_Nbkg1" )->getVal();
      //dE_N2  = ws->var( f2 + "_prime_Nbkg2" )->getVal();
      dE_N1  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      dE_a1  = ws->var( f2 + "_prime_a1" )->getVal();
      dE_a2  = ws->var( f2 + "_prime_a2" )->getVal();
    }
  else if ( f2 == "singleExp" )
    {
      sE_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      sE_a  = ws->var( f2 + "_prime_a" )->getVal();
      //sE_N  = ws->var( f1 + "_1_Nbkg" )->getVal();
      //sE_a  = ws->var( f1 + "_1_a" )->getVal();
    }
  else if ( f2 == "modExp" )
    {
      mE_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      mE_a  = ws->var( f2 + "_prime_a" )->getVal();
      mE_m  = ws->var( f2 + "_prime_m" )->getVal();
    }
  else if ( f2 == "singlePow" )
    {
      sP_N  = ws->var( f2 + "_prime_Nbkg" )->getVal();
      sP_a  = ws->var( f2 + "_prime_a" )->getVal();
    }
  else if ( f2 == "doublePow" )
    {
      dP_N   = ws->var( f2 + "_prime_Nbkg" )->getVal();
      dP_f   = ws->var( f2 + "_prime_f" )->getVal();
      dP_a1  = ws->var( f2 + "_prime_a1" )->getVal();
      dP_a2  = ws->var( f2 + "_prime_a2" )->getVal();
    }
  else if ( f2 == "poly2" )
    {
      pN = ws->var( f2 + "_prime_Nbkg" )->getVal();
      pC = ws->var( f2 + "_prime_pC" )->getVal();
      p0 = ws->var( f2 + "_prime_p0" )->getVal();
      p1 = ws->var( f2 + "_prime_p1" )->getVal();
    }
  else if ( f2 == "poly3" )
    {
      pN = ws->var( f2 + "_prime_Nbkg" )->getVal();
      pC = ws->var( f2 + "_prime_pC" )->getVal();
      p0 = ws->var( f2 + "_prime_p0" )->getVal();
      p1 = ws->var( f2 + "_prime_p1" )->getVal();
      p2 = ws->var( f2 + "_prime_p2" )->getVal();
    }
  else if ( f2 == "poly4" )
    {
      pN = ws->var( f2 + "_prime_Nbkg" )->getVal();
      pC = ws->var( f2 + "_prime_pC" )->getVal();
      p0 = ws->var( f2 + "_prime_p0" )->getVal();
      p1 = ws->var( f2 + "_prime_p1" )->getVal();
      p2 = ws->var( f2 + "_prime_p2" )->getVal();
      p3 = ws->var( f2 + "_prime_p3" )->getVal();
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  
  bool _badFit = false;
  int _countPass = 0;
  TRandom3* r = new TRandom3(0);
  TRandom3* r2 = new TRandom3(0);
  std::cout << "ENTERING LOOP" << std::endl;
  for ( int i = 0; i < ntoys; i++ )
    {
      //-----------------------------------------------------------
      //Defining Number of signal toys; frac = S/B = 0,1,2,3,4,5...
      //-----------------------------------------------------------
      int stoys = int(frac*f1Int*npoints);
      double npointsP = r->PoissonD( npoints );
      Nbkg.setVal(npointsP);
      double stoysP   = r2->PoissonD( stoys );
      std::cout << "[INFO]: it-> " << i << " ======> stoys: " << stoysP << " Nbkg: " << npointsP << " nBkG SR: " << f1Int*npoints << std::endl;
      //-------------------------
      //Set Signal pdf paramaters
      //-------------------------
      Ns.setVal( stoysP );//Note initial paramenter is the true number of signal events
      ws->var(gaussTag + "_frac")->setVal( gausFrac );
      ws->var(gaussTag + "_mu")->setVal( gausMu );
      ws->var(gaussTag + "_sigma1")->setVal( gausSigma1 );
      ws->var(gaussTag + "_sigma2")->setVal( gausSigma2 );
      
      //-------------------------------------------
      //Fixing Signal pdf paramaters, float only Ns
      //-------------------------------------------
      ws->var(gaussTag + "_frac")->setConstant(kTRUE);
      ws->var(gaussTag + "_mu")->setConstant(kTRUE);
      ws->var(gaussTag + "_sigma1")->setConstant(kTRUE);
      ws->var(gaussTag + "_sigma2")->setConstant(kTRUE);

      //--------------------------------
      //Update f2 (test pdf) paramenters
      //--------------------------------
      //ws->var( f2 + "_2_Nbkg" )->setRange(0, 1e3);
      if ( f2 == "doubleExp" )
      	{
	  /*
	  ws->var( f2 + "_2_Nbkg1")->setVal( dE_N1 );
	  ws->var( f2 + "_2_Nbkg2")->setVal( dE_N2 );
	  ws->var( f2 + "_2_a1")->setVal( dE_a1 );
	  ws->var( f2 + "_2_a2")->setVal( dE_a2 );
	  */
	  /*
	  ws->var( f2 + "_2_Nbkg")->setVal( dE_N1 );
	  ws->var( f2 + "_2_a1")->setVal( dE_a1 );
	  ws->var( f2 + "_2_a2")->setVal( dE_a2 );
	  */
	  ws->var( tag2 + "_a1")->setVal( dE_a1 );
	  ws->var( tag2 + "_a2")->setVal( dE_a2 );
	
	}
      else if ( f2 == "singleExp" )
	{
	  ws->var( tag2 + "_a" )->setVal( sE_a );
	}
      else if ( f2 == "modExp" )
	{
	 // ws->var( f2 + "_2_Nbkg" )->setVal( mE_N );
	 // ws->var( f2 + "_2_a" )->setVal( mE_a );
	 // ws->var( f2 + "_2_m" )->setVal( mE_m );
	  ws->var( tag2 + "_a" )->setVal( mE_a );
	  ws->var( tag2 + "_m" )->setVal( mE_m );
	}
      else if ( f2 == "singlePow" )
	{
	  //ws->var( f2 + "_2_Nbkg" )->setVal( sP_N );
	  //ws->var( f2 + "_2_a" )->setVal( sP_a );
	  ws->var( tag2 + "_a" )->setVal( sP_a );
	}
      else if ( f2 == "doublePow" )
	{
	/*
	  ws->var( f2 + "_2_Nbkg" )->setVal( dP_N );
	  ws->var( f2 + "_2_f" )->setVal( dP_f );
	  ws->var( f2 + "_2_a1" )->setVal( dP_a1 );
	  ws->var( f2 + "_2_a2" )->setVal( dP_a2 );
	*/
	  ws->var( tag2 + "_f" )->setVal( dP_f );
	  ws->var( tag2 + "_a1" )->setVal( dP_a1 );
	  ws->var( tag2 + "_a2" )->setVal( dP_a2 );
	}
      else if ( f2 == "poly2" )
	{
	/*
	  ws->var( f2 + "_2_Nbkg" )->setVal( pN );
	  ws->var( f2 + "_2_pC" )->setVal( pC );
	  ws->var( f2 + "_2_p0" )->setVal( p0 );
	  ws->var( f2 + "_2_p1" )->setVal( p1 );
	*/
	  ws->var( tag2 + "_pC" )->setVal( pC );
	  ws->var( tag2 + "_p0" )->setVal( p0 );
	  ws->var( tag2 + "_p1" )->setVal( p1 );
	}
      else if ( f2 == "poly3" )
	{
	/*
	  ws->var( f2 + "_2_Nbkg" )->setVal( pN );
	  ws->var( f2 + "_2_pC" )->setVal( pC );
	  ws->var( f2 + "_2_p0" )->setVal( p0 );
	  ws->var( f2 + "_2_p1" )->setVal( p1 );
	  ws->var( f2 + "_2_p2" )->setVal( p2 );
	*/
	  ws->var( tag2 + "_pC" )->setVal( pC );
	  ws->var( tag2 + "_p0" )->setVal( p0 );
	  ws->var( tag2 + "_p1" )->setVal( p1 );
	  ws->var( tag2 + "_p2" )->setVal( p1 );
	
	}
      else if ( f2 == "poly4" )
	{
	  /*
	  ws->var( f2 + "_2_Nbkg" )->setVal( pN );
	  ws->var( f2 + "_2_pC" )->setVal( pC );
	  ws->var( f2 + "_2_p0" )->setVal( p0 );
	  ws->var( f2 + "_2_p1" )->setVal( p1 );
	  ws->var( f2 + "_2_p2" )->setVal( p2 );
	  ws->var( f2 + "_2_p3" )->setVal( p3 );
	*/
	  ws->var( tag2 + "_pC" )->setVal( pC );
	  ws->var( tag2 + "_p0" )->setVal( p0 );
	  ws->var( tag2 + "_p1" )->setVal( p1 );
	  ws->var( tag2 + "_p2" )->setVal( p1 );
	  ws->var( tag2 + "_p3" )->setVal( p1 );

	}
      else
	{
	  ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
	}

      delete data_toys;
      delete signal_toys;
      //------------------------------------
      //G e n e r a t i n g   B K G  t o y s
      //------------------------------------
      data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npointsP );
      //-------------------------------------------
      //G e n e r a t i n g   s i g n a l   t o y s
      //-------------------------------------------
      signal_toys = GenerateToys( signalPdf, mgg, stoysP );
      //appending signal toys to data_toys (bkg toys)
      data_toys->append( *signal_toys );
      data_toys->SetName("dataToysSignalAndBkg");
      //-------------------------------------------------------------
      //S+B FIT (Ns is the only parameter floated for the signal pdf)
      //-------------------------------------------------------------
      //bres_toys = sbModel->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );

      //Initialize fit parameter statuses
      _status     = -1;
      _status2    = -1;
      _status3    = -1;
      _covStatus  = -1;
      _covStatus2 = -1;
      _covStatus3 = -1;

      //--------------------------------------
      //Defining Negative Log Likelihood (NLL)
      //--------------------------------------
      RooNLLVar* nll = (RooNLLVar*)sbModel->createNLL( *data_toys, RooFit::Extended(kTRUE), RooFit::Range("Full") );
      //Defininf RooMinimizer Object;
      RooMinimizer m(*nll);
      m.setStrategy(2);
      m.setPrintLevel(-1);
      m.setMaxFunctionCalls(5000);
      m.setMinimizerType("Minuit2");
      
      m.minimize("Minuit2", "Migrad");
      //m.migrad(); 
      RooFitResult* r = m.save(); 
      _status    = r->status();
      _covStatus = r->covQual();
      r->SetName("rMigrad");
      ws->import( *r );

      //m.hesse();
      if(_status!=0)
	{
      m.minimize("Minuit2", "Hesse");
	}
      RooFitResult* r2 = m.save(); 
      _status2    = r2->status();
      _covStatus2 = r2->covQual();
      r2->SetName("rHesse");
      ws->import( *r2 );
	

     // m.minimize("Minuit2", "Migrad");
     // m.minimize("Minuit2", "Hesse");
     // m.minimize("Minuit2", "Migrad");
     // m.minimize("Minuit2", "Hesse");
     /* 
      m.minimize("Minuit2", "Migrad");
      r = m.save();
      _status    = r->status();
      _covStatus = r->covQual();
*/
	
      int loop_Max = 100;
      int status2_before = _status2;
      
      while((status2_before!=0)&&(loop_Max>0))
	{
	  m.minimize("Minuit2", "Migrad");
      	  r = m.save();
          _status    = r->status();
          _covStatus = r->covQual();
	  
	  if ( _status ==0 ) break;
	
	  m.minimize("Minuit2", "Hesse");
	  r2 = m.save(); 
	  _status2    = r2->status();
	  _covStatus2 = r2->covQual();
	  if ( _status2 == 0 ) break;
	  
          loop_Max --;
	}

      //if ( _status2 == -1 ) break;
      
      //m.minos();
      RooFitResult* r3 = m.save(); 
      _status3    = r3->status();
      _covStatus3 = r3->covQual();

      /*
      if( !( _status == 0 && _covStatus == 3 && _status2 == 0 && _covStatus == 3 ) )
	{
	  //_countPass++;
	  //continue;
      }
      */
      _countPass++;

      //---------------------------
      //Saving paramenters to TTree
      //---------------------------
      //_alpha_hat = ws->var( f2 + "_2_a" )->getVal();
      //_alpha_sigma = ws->var( f2 + "_2_a" )->getError();
      std::cout << "[INFO]: Getting parameters after fit" << std::endl;
      _Nbkg       = npointsP;
      _Nbkg_hat   = Nbkg.getVal();
      _Nbkg_sigma = Nbkg.getError();
      _Ns         = stoysP;
      _Ns_hat     = Ns.getVal();
      _Ns_sigma   = Ns.getError();
      
      //-----------------------
      //G e t t i n g   b i a s
      //-----------------------
      double Nsignal  = _Ns_hat;
      double Ns_Error = _Ns_sigma/Nsignal;
      bias =  (_Ns - _Ns_hat)/_Ns_sigma;
      _bias = _Ns - _Ns_hat;
      _statUn = _Ns_sigma;//tree variable
      _biasNorm = bias.getVal();//tree variable
      
      std::cout << "DEBUG DEBUG" << std::endl;
      NsignalError.setVal( Ns_Error );

      //appending values to bias dataset outputs
      data_bias.add( RooArgSet(bias) );
      data_Nse.add( RooArgSet(NsignalError) );

      //------------------------------------
      //F i l l i n g  o u t p u t   T r e e
      //------------------------------------
      std::cout << "before filling tree" << std::endl;
      //std::cout << "iteration:" << i << std::endl;
      outTree->Fill();
      //delete data_toys;
      //delete signal_toys;
      delete nll;
    }

  RooPlot* ffPlot = mgg.frame();
  ffPlot->SetName("SBFit");
  ffPlot->SetTitle("");
  data_toys->plotOn( ffPlot );
  sbModel->plotOn( ffPlot, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("low,high"));      
  ws->import( *ffPlot );
  
  std::cout << "[DEBUG]: out of toys loop" << std::endl; 
  RooPlot* f_mgg = mgg.frame();
  f_mgg->SetName("toys_plot");
  f_mgg->SetTitle("");
  data_toys->plotOn( f_mgg );
  //ws->pdf( tag2 )->plotOn( f_mgg, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("low,high"));
  
  //ws->pdf( tag2p )->plotOn( f_mgg, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Normalization(npoints, RooAbsReal::NumEvent) );  
  //ws->pdf( tag1 )->plotOn( f_mgg, RooFit::LineColor(kGreen), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Normalization(npoints, RooAbsReal::NumEvent) );
  //Plotting background only component of s+b fit
  sbModel->plotOn( f_mgg, RooFit::LineStyle(kDashed), RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Components(tag2));
  sbModel->plotOn( f_mgg, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full"));
  
  RooPlot* f_bias = bias.frame();
  f_bias->SetName("bias_plot");
  data_bias.plotOn( f_bias );
  
  RooPlot* f_Nse = NsignalError.frame();
  f_Nse->SetName("NsignalError_plot");
  data_Nse.plotOn( f_Nse );
  
  ws->import( bias );
  ws->import( data_bias );
  ws->import( data_Nse );
  ws->import( *f_mgg );
  //ws->import( *bres_toys );
  ws->import( *f_bias );
  ws->import( *f_Nse );
  
  outName = outName + "_Tree.root";
  TFile* _fout = new TFile( outName, "recreate");
  outTree->Write();
  std::cout << "[INFO]: Wrote Bias TTree " << std::endl;
  _fout->Close();
  std::cout << "[INFO]: Bias Tree File Closed " << std::endl;
  return ws;
};

RooDataSet* GenerateToys( RooAbsPdf* pdf, RooRealVar x, int ntoys = 100 )
{
  return pdf->generate( x, ntoys);
};

RooWorkspace* MakeSideBandFitAIC_2( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName, double* Nbkg,  double& AIC, double& AIC_2, double& AIC_3,double& fitStatus_1, double& fitStatus_2, double& fitStatus_3, double& fitStatus_4,  TString ffName = "doubleExp")
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 135, 160);
  mgg.setRange("Full", 103, 160);

  TString tag;
  if ( ffName == "doubleExp" )
    {
      //tag = MakeDoubleExpN1N2( "sideband_fit_doubleExp", mgg, *ws );
      tag = MakeDoubleExp( "sideband_fit_doubleExp", mgg, *ws );
      std::cout << "[INFO]: Running double exponential fit" << std::endl; 
    }
  else if ( ffName == "singleExp" )
    {
      tag = MakeSingleExp( "sideband_fit_singleExp", mgg, *ws );
      std::cout << "[INFO]: Running single exponential fit" << std::endl; 
    }
  else if ( ffName == "modExp" )
    {
      tag = MakeModExp( "sideband_fit_modExp", mgg, *ws );
      std::cout << "[INFO]: Running modified exponential fit" << std::endl; 
    }
  else if ( ffName == "singlePow" )
    {
      tag = MakeSinglePow( "sideband_fit_singlePow", mgg, *ws );
      std::cout << "[INFO]: Running single pow fit" << std::endl; 
    }
  else if ( ffName == "doublePow" )
    {
      tag = MakeDoublePow( "sideband_fit_doublePow", mgg, *ws );
      //tag = MakeDoublePowN1N2( "sideband_fit_doublePow", mgg, *ws );
      std::cout << "[INFO]: Running double pow fit" << std::endl; 
    }
  else if ( ffName == "poly2" )
    {
      tag = MakePoly2( "sideband_fit_poly2", mgg, *ws );
      std::cout << "[INFO]: Running poly2 fit" << std::endl; 
    }
  else if ( ffName == "poly3" )
    {
      tag = MakePoly3( "sideband_fit_poly3", mgg, *ws );
      std::cout << "[INFO]: Running poly3 fit" << std::endl; 
    }
  else if ( ffName == "poly4" )
    {
      tag = MakePoly4( "sideband_fit_poly4", mgg, *ws );
      std::cout << "[INFO]: Running poly4 fit" << std::endl; 
    }

  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }

  std::cout << "====================" << std::endl;
  std::cout << "[INFO]: ENTERING FIT" << std::endl;
  std::cout << "====================" << std::endl;
  //Sideband Fit
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  double n = data.sumEntries(" (mGammaGamma>103 && mGammaGamma<120) || (mGammaGamma>135 && mGammaGamma<160)");
  *Nbkg = n;
  ws->var( "sideband_fit_"+ffName+"_Nbkg")->setVal( n * 1.2);
  
  //RooFitResult* bres = ws->pdf( tag )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
   RooNLLVar* nll = (RooNLLVar*)ws->pdf( tag )->createNLL( data, RooFit::Extended(kTRUE), RooFit::Range("low,high") );
   //Defining RooMinimizer Object;
   RooMinimizer m(*nll);
   m.setStrategy(2);
   //m.setEps(1e-4);
   //m.setMaxIterations(10000);
   //m.setMaxFunctionCalls(10000);
   m.setPrintLevel(2);
   m.minimize("Minuit2", "Migrad");

   //m.migrad(); 
   RooFitResult* bres = m.save() ; 
   //_status    = r->status();
   //_covStatus = r->covQual();
  
   fitStatus_1 = bres->status();
   fitStatus_2 = bres->covQual();
 
   m.minimize("Minuit2", "Hesse");
   RooFitResult* r2 = m.save() ; 
   //_status2    = r2->status();
   //_covStatus2 = r2->covQual();
 
//print out the status here.... 
   fitStatus_3 = r2->status();
   fitStatus_4 = r2->covQual();

	int max_loop_AIC = 100;	
   while(fitStatus_3!=0 && max_loop_AIC>0)
	{
   	m.minimize("Minuit2", "Migrad");
	bres = m.save();
	fitStatus_1 = bres->status();
   	fitStatus_2 = bres->covQual();
	if(fitStatus_1==0)	break;

   	m.minimize("Minuit2", "Hesse");
	r2 = m.save();
	fitStatus_3 = r2->status();
   	fitStatus_4 = r2->covQual();
	if(fitStatus_3==0) 	break;
	//{
	//	fitStatus_1 = fitStatus_3;
	//	fitStatus_2 = fitStatus_4;
	//	break;
	//}
	max_loop_AIC --;
	}
 	
	if(fitStatus_3==0)
	{
		fitStatus_1 = fitStatus_3;
		fitStatus_2 = fitStatus_4;
	
	} 
  std::cout << "===================" << std::endl;
  std::cout << "[INFO]: LEAVING FIT" << std::endl;
  std::cout << "===================" << std::endl;
  
  bres->SetName( tag + "_b_fitres" );
  std::cout << "[INFO]: PASS FIT" << std::endl;
  //---------------------
  //g e t t i n g   n l l 
  //---------------------
  double minNll = bres->minNll();
  //RooAbsReal* nll = ws->pdf( tag )->createNLL(data, RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  //RooAbsReal* nll = ws->pdf( tag )->createNLL(data, RooFit::Extended(kTRUE), RooFit::Range("Full") );
  std::cout << "nll_nll->" << nll->getVal() << std::endl;
  std::cout << "minNll->" << minNll << std::endl;
  RooArgSet* floatPars = ws->pdf( tag )->getParameters(data);
  double K = floatPars->getSize() - 1.;
  std::cout << "K -> " << K << std::endl;
  //double n = data.sumEntries(" (mgg>103 && mgg<120) || (mgg>135 && mgg<160)");

  if ( n-K-1 > 0 )
    {
      std::cout << "n -> " << n << std::endl;
      AIC = 2*minNll + 2*K + 2*K*(K+1)/(n-K-1);
      AIC_2 = 2*minNll + 2*K;// + 2*K*(K+1)/(n-K-1);
      AIC_3 = 2*minNll;// + 2*K + 2*K*(K+1)/(n-K-1);
      std::cout << "AIC: " << AIC << std::endl;
      std::cout << "AIC_2: " << AIC_2 << std::endl;
      std::cout << "AIC_3: " << AIC_3 << std::endl;
    }
  else
    {
      std::cerr << "[ERROR]: (n-K-1) zero or negative!: " << (n-K-1) << "; setting AIC to 999 " << std::endl;
      AIC   = 999.;
      AIC_2 = 999.;
      AIC_3 = 999.;
    }
  /*
  RooPlot* fns = ws->var("sideband_fit_doubleExp_a1")->frame( );
  nll->plotOn( fns, RooFit::LineColor(kBlue) );
  fns->SetName("nll_trick");
  ws->import( *fns );
  */
  ws->import( *bres );
  
  RooPlot *fmgg = mgg.frame();
  RooDataSet* dataCut = (RooDataSet*) data.reduce(RooFit::Name("dataCut"),RooFit::SelectVars(RooArgSet(mgg)),RooFit::CutRange("low"));
  RooDataSet* dataHigh = (RooDataSet*) data.reduce(RooFit::Name("dataHigh"),RooFit::SelectVars(RooArgSet(mgg)),RooFit::CutRange("high"));
  dataCut->append(*dataHigh);
  data.plotOn(fmgg,RooFit::Invisible());
  dataCut->plotOn(fmgg);
  ws->pdf( tag )->plotOn(fmgg,RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("low,high"));
  ws->pdf( tag )->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  
  /*
  double top_y = fmgg->GetMaximum();
  int n_g = 20;
  TGraph *grshade = new TGraph(2*n_g);
  for ( int ig = 0; ig < n_g; ig++ )
    {
      grshade->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
      grshade->SetPoint(n_g + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
    }
  
  grshade->SetFillStyle(3013);
  grshade->SetFillColor(16);
  grshade->Draw("f same");
  */
  fmgg->SetName( tag + "_frame" );
  ws->import( *fmgg );
  
  RooPlot* pdfFrame = mgg.frame();
  ws->pdf( tag )->plotOn( pdfFrame, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("Full") );
  pdfFrame->SetName( tag+"_pdfframe" );
  ws->import( *pdfFrame );
  ws->import( mgg );

  std::cout <<ffName<<" AIC fit status:  "<< fitStatus_1 << "," << fitStatus_2 << "," << fitStatus_3 << "," << fitStatus_4 << std::endl;
  return ws;
};

RooWorkspace* SelectBinning( TH1F* mggData, TString mggName, TString f1, TString f2, int ntoys, int npoints, TString outName )
{
  TFile* fout = new TFile( outName + "_biasTree.root", "RECREATE");
  RooRandom::randomGenerator()->SetSeed( 0 );
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 135, 160);
  mgg.setRange("Full", 103, 160);
  //HighRes Signal Region
  mgg.setRange("sig", 122.08, 128.92);
  
  TString tag1, tag2;
  if ( f1 == "doubleExp" )
    {
      tag1 = MakeDoubleExpN1N2( f1 + "_1", mgg, *ws );
      tag2 = MakeDoubleExpN1N2( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "singleExp" )
    {
      tag1 = MakeSingleExp( f1 + "_1", mgg, *ws );
      tag2 = MakeSingleExp( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "modExp" )
    {
      tag1 = MakeModExp( f1 + "_1", mgg, *ws );
      tag2 = MakeModExp( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "singlePow" )
    {
      tag1 = MakeSinglePow( f1 + "_1", mgg, *ws );
      tag2 = MakeSinglePow( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "doublePow" )
    {
      tag1 = MakeDoublePow( f1 + "_1", mgg, *ws );
      tag2 = MakeDoublePow( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "poly2" )
    {
      tag1 = MakePoly2( f1 + "_1", mgg, *ws );
      tag2 = MakePoly2( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "poly3" )
    {
      tag1 = MakePoly3( f1 + "_1", mgg, *ws );
      tag2 = MakePoly3( f1 + "_1_clone", mgg, *ws );
    }
  else if ( f1 == "poly4" )
    {
      tag1 = MakePoly4( f1 + "_1", mgg, *ws );
      tag2 = MakePoly4( f1 + "_1_clone", mgg, *ws );
    }
  else
    {
      std::cout << "[ERROR]: fit option not recognized. QUITTING PROGRAM" << std::endl;
      exit (EXIT_FAILURE);
    }
  
  //-------------------
  //Creatin RooDataHist
  //-------------------
  RooDataHist data( "data", "my_data", RooArgList(mgg), mggData );
  TString totalEntriesStr = Form("(%s>%.2f && %s<%.2f)", mggName.Data(), 103., mggName.Data(), 160.);
  TString sbEntriesStr = Form("(%s>%.2f && %s<%.2f) || (%s>%.2f && %s<%.2f)", mggName.Data(), 103., mggName.Data(), 120., mggName.Data(), 135., mggName.Data(), 160.);
  //npoints = data.sumEntries( totalEntriesStr );
  int sbpoints =  data.sumEntries( sbEntriesStr );

  // Represent data in dh as pdf in x
  RooHistPdf histpdf("histpdf", "histpdf", mgg, data, 0);
  RooAbsReal* fHistInt   = histpdf.createIntegral( mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  //mc toys dataHist
  RooDataHist* data_toys;
  
  std::cout << "=====================" << std::endl;
  std::cout << "[INFO]: total Str: " << totalEntriesStr << std::endl;
  std::cout << "[INFO]: sb Str: " << sbEntriesStr << std::endl;
  std::cout << "[INFO]: mgg variable name: " << mggName << std::endl;
  std::cout << "[INFO]: sideband entries: " << sbpoints << std::endl;
  std::cout << "[INFO]: total entries: "    << npoints << std::endl;
  std::cout << "=====================" << std::endl;
  
  TTree* biasTree = new TTree("biasTree", "Tree containing bias test information");
  double alpha_hat, alpha_true, n_hat, n_true, sigInt_hat, sigInt_true;
  double alpha_sigma, n_sigma, sigInt_sigma, sigInt_sigma2;
  double _status, _covStatus, _status2, _covStatus2, _status3, _covStatus3;
  double intTF1, intErr, intTF1_tot;
  biasTree->Branch("alpha_hat", &alpha_hat, "alpha_hat/D");
  biasTree->Branch("alpha_true", &alpha_true, "alpha_true/D");
  biasTree->Branch("alpha_sigma", &alpha_sigma, "alpha_sigma/D");
  biasTree->Branch("n_hat", &n_hat, "n_hat/D");
  biasTree->Branch("n_true", &n_true, "n_true/D");
  biasTree->Branch("n_sigma", &n_sigma, "n_sigma/D");
  biasTree->Branch("sigInt_hat", &sigInt_hat, "sigInt_hat/D");
  biasTree->Branch("sigInt_true", &sigInt_true, "sigInt_true/D");
  biasTree->Branch("sigInt_sigma", &sigInt_sigma, "sigInt_sigma/D");
  biasTree->Branch("sigInt_sigma2", &sigInt_sigma2, "sigInt_sigma2/D");
  biasTree->Branch("intTF1", &intTF1, "intTF1/D");
  biasTree->Branch("intErr", &intErr, "intErr/D");
  biasTree->Branch("intTF1_tot", &intTF1_tot, "intTF1_tot/D");
  biasTree->Branch("status", &_status,"status/D");
  biasTree->Branch("covStatus", &_covStatus,"covStatus/D");
  biasTree->Branch("status2", &_status2,"status2/D");
  biasTree->Branch("covStatus2", &_covStatus2,"covStatus2/D");
  biasTree->Branch("status3", &_status3,"status3/D");
  biasTree->Branch("covStatus3", &_covStatus3,"covStatus3/D");
  
  
  RooFitResult* bres;
  RooAbsReal* fIntegral2;
  //npoints = 1000;
  RooDataHist* data_toys2 = histpdf.generateBinned( mgg, npoints );
  ws->var("singleExp_1_clone_Nbkg")->setVal( npoints );
  RooFitResult* bres2 = ws->pdf( tag2 )->fitTo( *data_toys2, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  RooAbsReal* fIntTrue = ws->pdf( tag2 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  double alpha_clone = ws->var( f1 + "_1_clone_a")->getVal();
  
  TF1* myPdf;
  TRandom3* rnd = new TRandom3( 0 );
  //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator( "ADAPTIVE" );
  
  RooArgSet* paramSet = ws->pdf( tag2 )->getParameters( RooArgSet(mgg) );
  myPdf = ws->pdf( tag2 )->asTF( RooArgList(mgg), RooArgList(*paramSet) );
  //TMatrixDSym covMatrix = bres2->covarianceMatrix();
  TMatrixDSym covMatrix = bres2->correlationMatrix();
  std::cout << "Number Of Elements: " << covMatrix.GetNoElements() << std::endl;
  /*const double* pData = covMatrix.GetMatrixArray();
  for ( int j = 0; j < covMatrix.GetNoElements(); j++ )
    {
      std::cout << "par " << j << " = " << pData[j] << std::endl;
    }
  
  double* params = myPdf->GetParameters();
  */
  //sigInt_true = myPdf->Integral( 103, 160, 1e-15 );
  //sigInt_true = fIntTrue->getVal();
  alpha_true  = alpha_clone;
  double normTrueInt = fIntTrue->getVal();
  double* params;
  for ( int i = 0; i < ntoys; i++ )
    {
      n_true = rnd->PoissonD( (double)npoints );
      //n_true = (double)npoints;//avoid poisson fluctuations
      sigInt_true = n_true*normTrueInt;
      data_toys = ws->pdf( tag2 )->generateBinned( mgg, n_true );
      /*
      std::cout << "===========================" << std::endl;
      std::cout << "[INFO]: Setting N_bkg-> " << npoints << " , ntrue: " << n_true << "; and alpha-> " << alpha_clone << std::endl;
      std::cout << "===========================" << std::endl;
      */
      /*
      ws->var( f1 + "_1_Nbkg" )->setVal( n_true );
      ws->var( f1 + "_1_a" )->setVal( alpha_clone );
      bres = ws->pdf( tag1 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
      */
      _status = -1;
      _covStatus = -1;
      _status2 = -1;
      _covStatus2 = -1;
      _status3 = -1;
      _covStatus3 = -1;

      ws->var( f1 + "_1_Nbkg" )->setVal( n_true );
      ws->var( f1 + "_1_a" )->setVal( alpha_clone );
      //RooNLLVar* nll = (RooNLLVar*)ws->pdf( tag1 )->createNLL( *data_toys, RooFit::Extended(kTRUE), RooFit::Range("low,high") );
      RooNLLVar* nll = (RooNLLVar*)ws->pdf( tag1 )->createNLL( *data_toys, RooFit::Extended(kTRUE), RooFit::Range("Full") );
      RooMinimizer m(*nll);
      //m.fit("mhvr0");
      //ws->var( f1 + "_1_Nbkg" )->setVal( n_true );
      //ws->var( f1 + "_1_a" )->setVal( alpha_clone );
      m.migrad(); 
      RooFitResult* r = m.save() ; 
      _status    = r->status();
      _covStatus = r->covQual();
      
      m.hesse(); 
      RooFitResult* r2 = m.save() ; 
      _status2    = r2->status();
      _covStatus2 = r2->covQual();
      
      m.minos();
      RooFitResult* r3 = m.save() ; 
      _status3    = r3->status();
      _covStatus3 = r3->covQual();
      
      fIntegral2 = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
      RooAbsReal* f2int = ws->pdf( tag1 )->createIntegral(mgg, RooFit::Range("sig") );
      RooAbsReal* f3int = ws->pdf( tag1 )->createIntegral(mgg, RooFit::Range("Full") );
      RooAbsReal* f4int = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("Full") );

      paramSet = ws->pdf( tag1 )->getParameters( RooArgSet(mgg) );
      myPdf = ws->pdf( tag1 )->asTF( RooArgList(mgg), RooArgList(*paramSet) );
      /*covMatrix = bres->covarianceMatrix();*/
      covMatrix = r3->covarianceMatrix();
      params = myPdf->GetParameters();
      //Getting Integrals
      intTF1     = fIntegral2->getVal();
      //intErr     = myPdf->IntegralError( 103., 160., params, covMatrix.GetMatrixArray(), 1e-14 );
      //intTF1_tot = myPdf->Integral( 103., 160., 1e-15 );

      /*
	std::cout << "=======================" << i << "==========================" << std::endl;
      std::cout << "Int RooFit Normalized --> " << fIntegral2->getVal() << std::endl;
      std::cout << "Int TF1 normalized --> " << intTF1/intTF1_tot << std::endl;
      std::cout << "Int RooFit --> " << f2int->getVal() << std::endl;
      std::cout << "Int TF1 --> " << intTF1 << std::endl;
      std::cout << "Int RooFit Total --> " << f3int->getVal() << std::endl;
      std::cout << "Int RooFit Total Normalized --> " << f4int->getVal() << std::endl;
      std::cout << "Int TF1 Total --> " << intTF1_tot << std::endl;
      std::cout << "=================================================" << std::endl;
      */
      //alpha
      alpha_hat   = ws->var( f1 + "_1_a" )->getVal();
      alpha_sigma = ws->var( f1 + "_1_a" )->getError();
      //n
      n_hat       = ws->var( f1 + "_1_Nbkg" )->getVal();
      n_sigma     = ws->var( f1 + "_1_Nbkg" )->getError();
      //signal integral
      sigInt_hat = intTF1*n_hat;
      //Single Exp Uncertainty Calculation
      //122.08, 128.92
      //getting correlation matrix parameters from correlation matrix
      TMatrixDSym corrMatrix = r3->correlationMatrix();
      const double* pData = corrMatrix.GetMatrixArray();
      //propagation of uncertainty
      double eta    = exp(-alpha_hat*122.08) - exp(-alpha_hat*128.92);
      double etaMod = 128.92*exp(-alpha_hat*128.92) - 122.08*exp(-alpha_hat*122.08);
      double xi     = exp(-alpha_hat*103.) - exp(-alpha_hat*160.);
      double xiMod  = 103.*exp(-alpha_hat*103.) - 160.*exp(-alpha_hat*160.);
      double dIdN = eta/xi;
      double dIdAlpha = n_hat*( eta*xiMod/(xi*xi) + etaMod/xi );
      double sigmaI = sqrt(dIdN*dIdN*n_sigma*n_sigma + dIdAlpha*dIdAlpha*alpha_sigma*alpha_sigma + dIdN*dIdAlpha*pData[1]*n_sigma*alpha_sigma);

      //std::cout << "dIdN: " << dIdN << " , dIdAlpha: " << dIdAlpha << std::endl;
      //std::cout << "Int: " << sigInt_hat << " +/- " << sigmaI << std::endl;
      sigInt_sigma = sigmaI;
      sigInt_sigma2 = sqrt(dIdN*dIdN*n_sigma*n_sigma + dIdAlpha*dIdAlpha*alpha_sigma*alpha_sigma);
      biasTree->Fill();
    }
  //--------
  //INTEGRAL
  //--------
  RooAbsReal* fIntegral  = ws->pdf( tag1 )->createIntegral(mgg);
  
  
  std::cout << "test Int: " << fIntegral->getVal() << std::endl;
  std::cout << "test Int2: " << fIntegral2->getVal() << std::endl;
  std::cout << "hist Int: " << fHistInt->getVal() << std::endl;
  bres2->SetName( tag2 + "_b_fitRes");
  //ws->import( *bres );

  RooPlot* fmgg = mgg.frame();
  data_toys2->plotOn( fmgg );
  ws->pdf( tag2 )->plotOn( fmgg, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("low,high") );
  fmgg->SetName( tag2 + "frame" );
  ws->import( *fmgg );

  
  biasTree->Write();
  ws->Write("my_ws");
  fout->Close();
  return ws;
};
