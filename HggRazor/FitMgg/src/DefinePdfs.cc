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
#include <RooChebychev.h>
#include <RooMinuit.h>
//LOCAL INCLUDES
#include "DefinePdfs.hh"

TString MakeDoubleGauss( TString tag, RooRealVar& mgg, RooWorkspace& w )
{
  //------------------------------
  //C r e a t e  V a r i a b l e s
  //------------------------------
  RooRealVar* mu     = new RooRealVar( tag+"_gauss_mu", "#mu", 125, 0., 1000. );
  RooRealVar* sigma1 = new RooRealVar( tag+"_gauss_sigma1", "#sigma_{1}", 1.0, .0, 5000. );
  RooRealVar* sigma2 = new RooRealVar( tag+"_gauss_sigma2", "#sigma_{2}", 1.0, .0, 5000. );
  RooRealVar* frac   = new RooRealVar( tag+"_frac", "frac", 0.1, .0, 1.0 );
  RooRealVar* Ns     = new RooRealVar( tag+"_gauss_Ns", "N_{s}", 1, -1e8, 1e8 );
  //------------------
  //C r e a t e  p.d.f
  //------------------
  RooGaussian* gauss1    = new RooGaussian( tag+"_gauss1", "", mgg, *mu, *sigma1 );
  RooGaussian* gauss2    = new RooGaussian( tag+"_gauss2", "", mgg, *mu, *sigma2 );
  RooAddPdf* doublegauss = new RooAddPdf( tag + "_doublegauss", "", RooArgList( *gauss1, *gauss2 ), *frac );

  //------------------------------------
  //C r e a t e   E x t e n d e d  p.d.f
  //------------------------------------
  RooFormulaVar* Ns_Sq = new RooFormulaVar( tag + "Ns_Sq","","@0*@0", *Ns );
  TString ex_pdf_name          = tag+"_doublegauss_ext";
  RooExtendPdf* ex_doublegauss = new RooExtendPdf( ex_pdf_name, "extDgauss", *doublegauss, *Ns_Sq );
  w.import( *ex_doublegauss );
  
  return ex_pdf_name;
};

TString MakeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w)
{
  //------------------------------
  //C r e a t e  V a r i a b l e s
  //------------------------------
  RooRealVar* alpha1 = new RooRealVar( tag + "_dexp_alpha1", "#alpha_{1}", 0.6, -10, 10 );
  RooRealVar* alpha2 = new RooRealVar( tag + "_dexp_alpha2", "#alpha_{2}", 0.00001, -10, 10 );
  //--------------------------------------------
  //Square variables to avoid rising exponential
  //--------------------------------------------
  RooFormulaVar* asq1 = new RooFormulaVar( tag + "_dexp_alpha1Sq", "#alpha^{2}_{1}", "-1*@0*@0", *alpha1);
  RooFormulaVar* asq2 = new RooFormulaVar( tag + "_dexp_alpha2Sq", "#alpha^{2}_{2}", "-1*@0*@0", *alpha2);
  RooRealVar* frac    = new RooRealVar( tag+"_dexp_frac", "frac", 0.99, .0, 1.0 );
  RooRealVar* Nbkg    = new RooRealVar( tag+"_dexp_Nbkg", "N_{bkg}", 10, 1e-10, 1e8 );
  
  //------------------
  //C r e a t e  p.d.f
  //------------------
  RooExponential* exp1    = new RooExponential( tag+"_exp1", "", mgg, *asq1 );
  RooExponential* exp2    = new RooExponential( tag+"_exp2", "", mgg, *asq2 );
  RooAddPdf* doubleExp = new RooAddPdf( tag + "_doubleExp", "", RooArgList( *exp1, *exp2 ), *frac );

  //------------------------------------
  //C r e a t e   E x t e n d e d  p.d.f
  //------------------------------------
  TString ex_pdf_name          = tag+"_doubleExp_ext";
  RooExtendPdf* ex_doubleExp = new RooExtendPdf( ex_pdf_name, "extDexp", *doubleExp, *Nbkg);
  w.import( *ex_doubleExp );
  
  return ex_pdf_name;
};

TString MakeDoubleExpN1N2( TString tag, RooRealVar& mgg, RooWorkspace& w )
{
  RooRealVar* a1 = new RooRealVar( tag + "_a1", "", 0.6,-1.,1.);
  RooRealVar* a2 = new RooRealVar( tag + "_a2", "", 0.4/150,-1/150.,1/150.);
  
  RooFormulaVar* a1sq = new RooFormulaVar( tag + "a1sq","","-1*@0*@0", *a1);
  RooFormulaVar* a2sq = new RooFormulaVar( tag + "a2sq","","-1*@0*@0", *a2);
  
  RooRealVar* NBkg1 = new RooRealVar( tag + "Nbkg1","",10,-1e5,1e5);
  RooRealVar* NBkg2 = new RooRealVar( tag + "Nbkg2","",1,-1e5,1e5);
  
  RooFormulaVar* NBkg1Sq = new RooFormulaVar( tag + "Nbkg1Sq","","@0*@0", *NBkg1);
  RooFormulaVar* NBkg2Sq = new RooFormulaVar( tag + "Nbkg2Sq","","@0*@0", *NBkg2);

  RooExponential* e1 = new RooExponential( tag + "_e1","", mgg, *a1sq);
  RooExponential* e2 = new RooExponential( tag + "_e2","", mgg, *a2sq);

  TString pdfName = tag+"_pdf_dExp_N1N2";
  RooAddPdf* pdf = new RooAddPdf( pdfName,"", RooArgSet( *e1, *e2), RooArgSet( *NBkg1Sq, *NBkg2Sq) );
  w.import( *pdf );
  return pdfName;
};

TString MakeSingleExp( TString tag, RooRealVar& mgg, RooWorkspace& w )
{
  RooRealVar* a = new RooRealVar( tag + "_se_a", "", 0.6,-1.,1.);
  RooFormulaVar* asq = new RooFormulaVar( tag + "_se_asq","","-1*@0*@0", *a);
  RooRealVar* NBkg = new RooRealVar( tag + "se_Nbkg","",10,-1e5,1e5);
  RooFormulaVar* NBkg1Sq = new RooFormulaVar( tag + "se_Nbkg1Sq","","@0*@0", *NBkg);
  
  RooExponential* se = new RooExponential( tag + "_se","", mgg, *asq);
  
  TString pdfName = tag+"_pdf_sExp_ext";
  RooExtendPdf* ext_singleExp = new RooExtendPdf( pdfName, "extSexp", *se, *NBkg1Sq);
  w.import( *ext_singleExp );
  return pdfName;
};

TString MakeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_mexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *m1 = new RooRealVar(tag+"_mexp_m1","m_{1}",1.,0.,10.);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_mexp_Nbkg","N_{bkg}",10, -1e8, 1e8);
  RooGenericPdf *mexp1 = new RooGenericPdf(tag+"_mexp","mod_exp","exp(@0*(@1^@2))",RooArgList(*alpha1,mgg,*m1));
  
  TString pdfName = tag+"_mexp_ext";
  RooFormulaVar* NBkg1Sq = new RooFormulaVar( tag + "modExp_Nbkg1Sq","","@0*@0", *Nbkg );
  RooExtendPdf* modExp_Ext = new RooExtendPdf( pdfName, "modExp", *mexp1, *NBkg1Sq );
  w.import( *modExp_Ext );
  return pdfName;
};

TString MakeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w) 
{
  RooRealVar *alpha1  = new RooRealVar(tag+"_spow_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *Nbkg    = new RooRealVar(tag+"_spow_Nbkg","N_{bkg}",10,1,1E9);
  RooGenericPdf *pow1 = new RooGenericPdf(tag+"_spow","","@0^@1",RooArgList(mgg,*alpha1));
  TString pdfName     = tag+"_spow_ext";
  RooFormulaVar* NBkg1Sq = new RooFormulaVar( tag + "singlePow_Nbkg1Sq","","@0*@0", *Nbkg );
  w.import( *(new RooExtendPdf( pdfName,"",*pow1,*NBkg1Sq)) );
  return pdfName;
};

TString MakeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w)
{
  RooRealVar *alpha1  = new RooRealVar(tag+"_dpow_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *alpha2  = new RooRealVar(tag+"_dpow_alpha2","#alpha_{2}",-2,-10,-0.0001);
  RooRealVar *f       = new RooRealVar(tag+"_dpow_f","f",0.3,0.0001,0.9999);
  RooRealVar *Nbkg    = new RooRealVar(tag+"_dpow_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *pow1 = new RooGenericPdf(tag+"_dpow_pow1","","@0^@1",RooArgList(mgg,*alpha1));
  RooGenericPdf *pow2 = new RooGenericPdf(tag+"_dpow_pow2","","@0^@1",RooArgList(mgg,*alpha2));

  TString pdfName = tag+"_dpow_ext";
  RooAddPdf *add      = new RooAddPdf( tag+"_dpow","",*pow1,*pow2,*f);  
  w.import( *(new RooExtendPdf( pdfName,"",*add,*Nbkg)) );

  return pdfName;
};

TString MakePoly2(TString tag, RooRealVar& mgg,RooWorkspace& w)
{
  RooRealVar *pC = new RooRealVar(tag+"_pol2_pC","C",1,-100,100);
  RooRealVar *p0 = new RooRealVar(tag+"_pol2_p0","p_0",1,-100,100);
  RooRealVar *p1 = new RooRealVar(tag+"_pol2_p1","p_1",1,-100,100);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_pol2_Nbkg","N_{bkg}",100,1,1E9);

  RooFormulaVar *pCmod = new RooFormulaVar(tag+"_pol2_pCmod","@0*@0",*pC);
  RooFormulaVar *p0mod = new RooFormulaVar(tag+"_pol2_p0mod","@0*@0",*p0);
  RooFormulaVar *p1mod = new RooFormulaVar(tag+"_pol2_p1mod","@0*@0",*p1);

  RooBernstein* bern = new RooBernstein(tag+"_pol2","",mgg,RooArgList(*pCmod,*p0mod,*p1mod));
  //RooChebychev* chev = new RooChebychev(tag+"_pol2","",mgg,RooArgList(*pCmod,*p0mod,*p1mod));
  TString pdfName = tag+"_pol2_ext";
  w.import(*(new RooExtendPdf( pdfName,"",*bern,*Nbkg)));
  //w.import(*(new RooExtendPdf( pdfName,"",*chev,*Nbkg)));
  

  return pdfName;
};

TString MakePoly3(TString tag, RooRealVar& mgg,RooWorkspace& w)
{
  RooRealVar *pC = new RooRealVar(tag+"_pol3_pC","C",1,-100,100);
  RooRealVar *p0 = new RooRealVar(tag+"_pol3_p0","p_0",1.0,-100,100);
  RooRealVar *p1 = new RooRealVar(tag+"_pol3_p1","p_1",0,-100,100);
  RooRealVar *p2 = new RooRealVar(tag+"_pol3_p2","p_2",5,-100,100);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_pol3_Nbkg","N_{bkg}",2000,1,1E9);
  
  RooFormulaVar *pCmod = new RooFormulaVar(tag+"_pol3_pCmod","@0*@0",*pC);
  RooFormulaVar *p0mod = new RooFormulaVar(tag+"_pol3_p0mod","@0*@0",*p0);
  RooFormulaVar *p1mod = new RooFormulaVar(tag+"_pol3_p1mod","@0*@0",*p1);
  RooFormulaVar *p2mod = new RooFormulaVar(tag+"_pol3_p2mod","@0*@0",*p2);

  RooBernstein* bern = new RooBernstein(tag+"_pol3","",mgg,RooArgList(*pCmod,*p0mod,*p1mod,*p2mod));

  TString pdfName = tag+"_pol3_ext";
  w.import( *(new RooExtendPdf( pdfName,"",*bern,*Nbkg)) );
  
  return pdfName;
};
