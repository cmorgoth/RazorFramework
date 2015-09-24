//C++ INCLUDES
#include <vector>
//ROOT INCLUDES
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
//ROOFIT INCLUDES
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
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
//#include <RealVar.h>
//LOCAL INCLUDES
#include "RunTwoFitMgg.hh"
#include "DefinePdfs.hh"

RooWorkspace* DoubleGausFit( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  RooRealVar mgg( mggName, "m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(57);
  mgg.setRange( "signal", 120, 130. );
  //C r e a t e  doubleGaus
  TString tag = MakeDoubleGauss( "dGauss_signal", mgg, *ws );
  ws->var("dGauss_signal_gauss_Ns")->setVal( 1600 );
  
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import( *tree ) );
  ws->pdf( tag )->fitTo( data, RooFit::Strategy(0), RooFit::Extended( kTRUE ), RooFit::Range("signal") );
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
  double n1 = ws->var("sideband_fitNbkg1")->getVal();
  double n2 = ws->var("sideband_fitNbkg2")->getVal();
  double Nbkg = n1*n1 + n2*n2; 

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Ntotal: " << data.sumEntries( ) << std::endl;
  std::cout << "Nfit: " << Nfit->getVal() << std::endl;
  std::cout << "Nbkg: " << Nbkg << std::endl;
  std::cout << "sf: "  << Nbkg*sint->getVal()/N_sideband << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  
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
};

RooWorkspace* MakeSignalBkgFit( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName )
{
  RooWorkspace* ws = new RooWorkspace( "wsb", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  mgg.setRange("signal", 103, 160);

  //--------------------------------
  //I m p or t i n g   D a t a
  //--------------------------------
  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );
  
  //------------------------------------
  // C r e a t e   b k g  s h a p e
  //------------------------------------
  TString tag_bkg = MakeDoubleExpN1N2( "fullsb_fit_bkg", mgg, *ws );
  //------------------------------------
  // C r e a t e   b k g  s h a p e
  //------------------------------------
  TString tag_signal = MakeDoubleGauss( "fullsb_fit_signal", mgg, *ws );

  //RooAddPdf* model = new ( "model", "model", RooFit::RooArgList( *ws->pdf( tag_signal ), *ws->pdf( tag_bkg ) ) ) ;
  RooAddPdf* model = new RooAddPdf( "model", "model", RooArgSet( *ws->pdf( tag_signal ), *ws->pdf( tag_bkg ) ) ) ;

  //--------------------------------------
  //H i g g s   C o n s t r a i n s
  //--------------------------------------
  RooRealVar HiggsMass("HiggsMass","",128.7);
  RooRealVar HiggsMassError("HiggsMassError","",0.4);
  RooGaussian HiggsMassConstraint("HiggsMassConstraint","", *ws->var("fullsb_fit_signal_gauss_mu"),HiggsMass,HiggsMassError);

  if( forceSigma != -1 ) {
    ws->var("fullsb_fit_signal_gauss_sigma1")->setVal( forceSigma );
    ws->var("fullsb_fit_signal_gauss_sigma1")->setConstant(true);
    ws->var("fullsb_fit_signal_frac")->setVal( 1.0 );
    ws->var("fullsb_fit_signal_frac")->setConstant(true);
  }

  //---------------------
  //F i t   t o   D a t a
  //---------------------
  
  model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  RooFitResult* bres = model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );

  //-------------------------------------
  //P r o f i l e d   L i k e l i h o o d 
  //-------------------------------------
  /*
    profile in Ns to obtain significance.
  */
  RooArgSet poi   = RooArgSet( *ws->var("fullsb_fit_signal_gauss_Ns") );
  RooAbsReal* nll = model->createNLL(data);
  RooFormulaVar n2ll = RooFormulaVar("n2ll", "2*@0", RooArgList(*nll) );
  RooAbsReal* p2ll = n2ll.createProfile( poi );
    
  
  /*
  model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full"), RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)) );
  RooFitResult* bres = model->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full"), RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)) );
  */
  RooPlot *fmgg = mgg.frame();
  data.plotOn(fmgg);
  model->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
  model->plotOn(fmgg,RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed), RooFit::Range("low,high"),RooFit::NormRange("low,high"));
  fmgg->SetName( "fullsb_fit_frame" );
  ws->import( *model );
  ws->import( *bres );
  ws->import( *fmgg );

  RooPlot* fns = ws->var("fullsb_fit_signal_gauss_Ns")->frame( RooFit::Range(0, 20, true) );
  fns->SetMinimum(0);
  fns->SetMaximum(6);
  n2ll.plotOn( fns, RooFit::ShiftToZero(), RooFit::LineColor(kBlue) );
  p2ll->plotOn( fns, RooFit::LineColor(kBlack) );
  fns->SetName("nll_trick");
  ws->import( *fns );
  return ws;
};

RooWorkspace* MakeSideBandFitAIC( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName, double& AIC, TString ffName = "doubleExp" )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(57);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);

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
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg( mggName,"m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(57);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  mgg.setRange("sig", 122.08, 128.92);
  
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

  RooDataSet data( "data", "", RooArgSet(mgg), RooFit::Import(*tree) );

  //Sideband Fit (not working with poly2 and poly3)
  //ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
  //RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
  //FullFit
  ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );

  RooAbsReal* f1Integral = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  std::cout << "f1 Int: " << f1Integral->getVal() << std::endl;
  
  RooDataSet* data_toys;
  RooFitResult* bres_toys;
  double n; 
  RooAbsReal* fIntegral;
  RooAbsReal* fIntegral2;
  RooRealVar bias("bias", "bias", -1.0, 1.0, "");
  RooDataSet data_bias( "data_bias", "bias data", bias);
  bias.setBins(100);
  for ( int i = 0; i < ntoys; i++ )
    {
      //-----------------------------
      //G e n e r a t i n g   t o y s
      //-----------------------------
      data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npoints );
      //ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("low,high") );
      //RooFitResult* bres_toys = ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
      ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
      bres_toys = ws->pdf( tag2 )->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
      bres_toys->SetName("bres_toys");
      
      //n = data_toys->sumEntries(" (mgg>122.08 && mgg<128.92)");
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

RooWorkspace* DoBiasTestSignal( TTree* tree, TString mggName, TString f1, TString f2, int ntoys, int npoints )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );

  //Getting signal shape
  TFile* fsignal = new TFile("data/SM_MC_GluGluH_Fit.root", "read");
  RooWorkspace* wsignal = (RooWorkspace*)fsignal->Get("w_sFit");
  RooAbsPdf* signalPdf = wsignal->pdf("dGauss_signal_doublegauss");
  double gausFrac   =  wsignal->var("dGauss_signal_frac")->getVal();
  double gausMu     =  wsignal->var("dGauss_signal_gauss_mu")->getVal();
  double gausSigma1 =  wsignal->var("dGauss_signal_gauss_sigma1")->getVal();
  double gausSigma2 =  wsignal->var("dGauss_signal_gauss_sigma2")->getVal();

  delete fsignal;
  //fsignal->Close();
  
  RooRealVar mgg( mggName,"m_{#gamma#gamma}", 103, 160, "GeV" );
  mgg.setBins(38);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);
  mgg.setRange("sig", 122.08, 128.92);

  RooDataSet* signal_toys = GenerateToys( signalPdf, mgg, npoints );
  signal_toys->SetName("signal_toys");
  RooPlot* s_mgg = mgg.frame();
  s_mgg->SetName("signal_toys_plot");
  signal_toys->plotOn( s_mgg );
  ws->import( *signal_toys );
  ws->import( *s_mgg );
  
  
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
      tag1 = MakeDoublePow( "doublePow_2", mgg, *ws );
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
  ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  RooFitResult* bres = ws->pdf( tag1 )->fitTo( data, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
  
  RooAbsReal* f1Integral = ws->pdf( tag1 )->createIntegral(mgg, RooFit::NormSet(mgg), RooFit::Range("sig") );
  double f1Int = f1Integral->getVal();

  RooPlot* bFrame = mgg.frame();
  bFrame->SetName("bFitFrame");
  data.plotOn( bFrame );
  ws->pdf( tag1 )->plotOn( bFrame, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full") );
  ws->import( *bFrame );

  //-------------------------------
  //S i g n a l   +   B k g   P d f
  //-------------------------------
  TString gaussTag  = MakeDoubleGauss("doubleGauseSB", mgg, *ws );
  std::cout << "NS: " << ws->var("doubleGauseSB_gauss_Ns")->getVal() << std::endl;
  RooAddPdf* sbModel = new RooAddPdf( "sbModel", "sbModel", RooArgList( *ws->pdf(tag2), *ws->pdf(gaussTag) ) );
  ws->import( *sbModel );
  
  RooDataSet* data_toys;
  RooFitResult* bres_toys;
  double n; 
  RooAbsReal* fIntegral;
  RooAbsReal* fIntegral2;
  RooRealVar bias("bias", "bias", -2.0, 2.0, "");
  RooDataSet data_bias( "data_bias", "bias data", bias);

  RooRealVar NsignalError("NsignalError", "NsignalError", -.02, .02, "");
  RooDataSet data_Nse( "data_Nse", "data_Nse", NsignalError);
  
  bias.setBins(100);
  NsignalError.setBins(100);
  data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npoints );
  ws->pdf( tag2p )->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
  RooPlot* pFrame = mgg.frame();
  pFrame->SetName("pFitFrame");
  data_toys->plotOn( pFrame );
  ws->pdf( tag2p )->plotOn( pFrame, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full") );
  ws->pdf( tag2p )->plotOn( pFrame, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Normalization(npoints, RooAbsReal::NumEvent) );
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
  bool _badFit = false;
  for ( int i = 0; i < ntoys; i++ )
    {
      //-----------------------------
      //G e n e r a t i n g   t o y s
      //-----------------------------
      data_toys = GenerateToys( ws->pdf( tag1 ), mgg, npoints );
      double frac = 1.0;
      int stoys = int(frac*f1Int*npoints);
      std::cout << "[INFO]:======> stoys: " << stoys << std::endl;
      ws->var("doubleGauseSB_gauss_Ns")->setVal( sqrt(stoys) );
      ws->var("doubleGauseSB_frac")->setVal( gausFrac );
      ws->var("doubleGauseSB_gauss_mu")->setVal( gausMu );
      ws->var("doubleGauseSB_gauss_sigma1")->setVal( gausSigma1 );
      ws->var("doubleGauseSB_gauss_sigma2")->setVal( gausSigma2 );
      
      
      ws->var("doubleGauseSB_frac")->setConstant(kTRUE);
      ws->var("doubleGauseSB_gauss_mu")->setConstant(kTRUE);
      ws->var("doubleGauseSB_gauss_sigma1")->setConstant(kTRUE);
      ws->var("doubleGauseSB_gauss_sigma2")->setConstant(kTRUE);
      
      //ws->var("singleExpse_Nbkg")->setVal( sqrt(npoints) );
      //if ( f2 == "doubleExp" ) ws->var("doubleExpNbkg2")->setVal(0.01);
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
      std::cout << "debug git" << std::endl;
      signal_toys = GenerateToys( signalPdf, mgg, stoys );
      data_toys->append( *signal_toys );
      
      //sbModel->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Range("Full") );
      //bres_toys = sbModel->fitTo( *data_toys, RooFit::Strategy(2), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
      bres_toys = sbModel->fitTo( *data_toys, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("Full") );
      bres_toys->SetName("bres_toys");
      
      double Nsignal  = pow(ws->var("doubleGauseSB_gauss_Ns")->getVal(), 2);
      double Ns_Error = 2.0*ws->var("doubleGauseSB_gauss_Ns")->getError()/Nsignal;
      bias =  (Nsignal - double(stoys))/(double)stoys;
      NsignalError.setVal( Ns_Error );
      //std::cout << "bias: " << bias.getVal() << std::endl;
      //if ( fabs( bias.getVal() ) > .8  )
      if ( bres_toys->status() != 0  )
	{
	  std::cout << "FIT STATUS: " << bres_toys->status() << ", covQual: " << bres_toys->covQual()
		    << " Bias: " << bias.getVal() << std::endl;
	  _badFit = true;
	  break;
	}
      data_bias.add( RooArgSet(bias) );
      data_Nse.add( RooArgSet(NsignalError) );
    }
  
  RooPlot* f_mgg = mgg.frame();
  f_mgg->SetName("toys_plot");
  data_toys->plotOn( f_mgg );
  //ws->pdf( tag2 )->plotOn( f_mgg, RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("low,high"));
  
  ws->pdf( tag2p )->plotOn( f_mgg, RooFit::LineColor(kBlue), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Normalization(npoints, RooAbsReal::NumEvent) );  
  ws->pdf( tag1 )->plotOn( f_mgg, RooFit::LineColor(kGreen), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Normalization(npoints, RooAbsReal::NumEvent) );
  //Plotting background only component of s+b fit
  sbModel->plotOn( f_mgg, RooFit::LineStyle(kDashed), RooFit::LineColor(kViolet), RooFit::Range("Full"), RooFit::NormRange("Full"), RooFit::Components(tag2));
  sbModel->plotOn( f_mgg, RooFit::LineColor(kRed), RooFit::Range("Full"), RooFit::NormRange("Full"));

  RooPlot* f_bias = bias.frame();
  f_bias->SetName("bias_plot");
  data_bias.plotOn( f_bias );

  RooPlot* f_Nse = NsignalError.frame();
  f_Nse->SetName("NsignalError");
  data_Nse.plotOn( f_Nse );
  
  ws->import( mgg );
  ws->import( bias );
  ws->import( data_bias );
  ws->import( *f_mgg );
  ws->import( *bres_toys );
  ws->import( *f_bias );
  ws->import( *f_Nse );
  std::cout << "[INFO]: npoints-> " << npoints << std::endl;
  
  return ws;
};
RooDataSet* GenerateToys( RooAbsPdf* pdf, RooRealVar x, int ntoys = 100 )
{
  return pdf->generate( x, ntoys );
};

RooWorkspace* MakeSideBandFitAIC_2( TTree* tree, float forceSigma, bool constrainMu, float forceMu, TString mggName, double& AIC, double& AIC_2, double& AIC_3, TString ffName = "doubleExp" )
{
  RooWorkspace* ws = new RooWorkspace( "ws", "" );
  
  RooRealVar mgg(mggName,"m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(57);
  mgg.setRange("low", 103, 120);
  mgg.setRange("high", 131, 160);

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
  AIC_2 = 2*minNll + 2*K;// + 2*K*(K+1)/(n-K-1);
  AIC_3 = 2*minNll;// + 2*K + 2*K*(K+1)/(n-K-1);
  std::cout << "AIC: " << AIC << std::endl;
  std::cout << "AIC_2: " << AIC_2 << std::endl;
  std::cout << "AIC_3: " << AIC_3 << std::endl;
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
