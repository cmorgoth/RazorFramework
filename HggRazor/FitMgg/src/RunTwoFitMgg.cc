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
//#include <RealVar.h>
//LOCAL INCLUDES
#include "RunTwoFitMgg.hh"
#include "DefinePdfs.hh"

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
  RooFitResult* bres = ws->pdf( tag3 )->fitTo( data, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("low,high") );
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


  //PROFILE
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

  RooPlot* fns = ws->var("fullsb_fit_signal_gauss_Ns")->frame( RooFit::Range(0,20, true) );
  fns->SetMinimum(0);
  fns->SetMaximum(6);
  n2ll.plotOn( fns, RooFit::ShiftToZero(), RooFit::LineColor(kBlue) );
  p2ll->plotOn( fns, RooFit::LineColor(kBlack) );
  fns->SetName("nll_trick");
  ws->import( *fns );
  return ws;
};
