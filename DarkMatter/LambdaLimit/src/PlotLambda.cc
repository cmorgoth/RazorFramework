/*
Computes lambda cutoff scale limit for the different couplings (g)
and using R_lambda = 0.8
*/
//C++ INCLUDES
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
//ROOT INCLUDES
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TPad.h"
#include "TMath.h"
//LOCAL INCLUDES
#include "PlotLambda.hh"

const double max_y_sd = 2e6;
const double min_y_sd = 1e2;

const double max_y_si = 2e6;
const double min_y_si = 1e2;

void  PlotLambda( )
{
  
  TString RootFileName;
  TString outputDir = "Plots";
  
  ///////////////////////
  ////CreatingPLOT AV///
  //////////////////////
  
  //Getting Razor DM TGraphs
  TFile* f3 = new TFile("inputData/no_btag/SD_Graphs.root");//Root file containing Razor TGraphs
  TGraph* Lc = (TGraph*)f3->Get("RazorDM_AV_8TeV_Combined_Toys");
  TGraph* Lcom = (TGraph*)f3->Get("RazorDM_AV_8TeV_Combined_Expected");
  TGraph* Lup = (TGraph*)f3->Get("RazorDM_AV_8TeV_Combined_Expected_up");
  TGraph* Ldown = (TGraph*)f3->Get("RazorDM_AV_8TeV_Combined_Expected_down");
  
  //Getting Direct Detection SI TGraphs
  //TFile* f4 = new TFile("inputData/no_btag/DirectDetection_SD_Lambda_Graphs.root");
  TFile* f4 = new TFile("inputData/no_btag/SD_Lambda_Graphs_FINAL.root");
  TGraph* COUPP_2012_SD = (TGraph*)f4->Get("COUPP_2012_SD");
  TGraph* SIMPLE_2012_SD = (TGraph*)f4->Get("SIMPLE_2012_SD");
  TGraph* SUPER_K_SD = (TGraph*)f4->Get("SUPER-K_SD");
  TGraph* ICE_CUBE_WW = (TGraph*)f4->Get("ICE_CUBE_WW");
  TGraph* XENON_SD = (TGraph*)f4->Get("XENON_SD_FINAL");
  TGraph* PICO_SD = (TGraph*)f4->Get("PICO");
  TGraph* PICASSO_SD = (TGraph*)f4->Get("PICASSO");
  
  TString labelX = "M_{#chi} (GeV)";
  TString labelY = "#Lambda (GeV)";
  int n = Lcom->GetN();
  double* x = Lup->GetX();
  double* yup = Lup->GetY();
  double* ydown = Ldown->GetY();
  TGraph *grshade = new TGraph(2*n);
  for(int i = 0; i < n; i++) {
    grshade->SetPoint(i,x[i],ydown[i]);
    grshade->SetPoint(n+i,x[n-i-1],yup[n-i-1]);
  }
  grshade->SetFillStyle(1001);
  grshade->SetFillColor(kOrange);
  grshade->SetLineColor(kOrange);
  grshade->SetTitle("");
  
  //TCanvas *c = new TCanvas("c", "c",2119,33,800,700);
  TCanvas *c = new TCanvas("c", "c", 0, 0, 800, 700);
  c->SetHighLightColor(2);
  c->Range(-0.6543224,-1290.871,3.177829,8696.391);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLogx();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.12);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);  
  c->SetTickx();
  c->SetTicky();
  
  TMultiGraph *mg0 = new TMultiGraph();
  mg0->Add(grshade);
  
  mg0->Draw("af");
  mg0->GetYaxis()->SetRangeUser( min_y_sd, max_y_sd );
  mg0->GetXaxis()->SetRangeUser( 0.001, 1000 );
  mg0->GetXaxis()->SetTitle(labelX);
  mg0->GetYaxis()->SetTitle(labelY);
  mg0->GetYaxis()->SetTitleSize( 0.05 );
  mg0->GetXaxis()->SetTitleSize( 0.05 );
  mg0->GetYaxis()->SetLabelSize( 0.04 );
  mg0->GetXaxis()->SetLabelSize( 0.04 );
  mg0->GetXaxis()->SetTitleOffset(1.1);
  mg0->GetYaxis()->SetTitleOffset(1.05);
  
  Lc->SetLineWidth(3);
  Lc->SetLineStyle(1);
  Lc->SetLineColor(kBlue);
  
  Lcom->SetLineColor(kRed-7);
  Lcom->SetLineWidth(3);
  Lcom->SetLineStyle(1);
  
  Lup->SetLineColor(kRed-7);
  Lup->SetLineWidth(1);
  Lup->SetLineStyle(1);
  
  Ldown->SetLineColor(kRed-7);
  Ldown->SetLineWidth(1);
  Ldown->SetLineStyle(1);
  
  COUPP_2012_SD->SetLineColor( kBlue-7 );
  COUPP_2012_SD->SetLineWidth( 3 );
  COUPP_2012_SD->SetLineStyle( 6 );
  
  SIMPLE_2012_SD->SetLineColor( kOrange-3 );
  SIMPLE_2012_SD->SetLineWidth( 3 );
  SIMPLE_2012_SD->SetLineStyle( 2 );
  
  SUPER_K_SD->SetLineColor( kViolet-2 );
  SUPER_K_SD->SetLineWidth( 3 );
  SUPER_K_SD->SetLineStyle( 6 );
  
  ICE_CUBE_WW->SetLineColor( kAzure+1 );
  ICE_CUBE_WW->SetLineWidth( 3 );
  ICE_CUBE_WW->SetLineStyle( 4 );

  XENON_SD->SetLineColor( kAzure+6 );
  XENON_SD->SetLineWidth( 3 );
  XENON_SD->SetLineStyle( 5 );
  
  PICO_SD->SetLineColor( kGreen-6 );
  PICO_SD->SetLineWidth( 3 );
  PICO_SD->SetLineStyle( 2 );

  PICASSO_SD->SetLineColor( kBlue-2 );
  PICASSO_SD->SetLineWidth( 3 );
  PICASSO_SD->SetLineStyle( 2 );
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add( Lup );
  mg->Add( Lcom );
  mg->Add( Ldown ); 
  mg->Add( Lc );
  mg->Add( COUPP_2012_SD );
  mg->Add( SIMPLE_2012_SD );
  mg->Add( SUPER_K_SD );
  mg->Add( ICE_CUBE_WW );
  mg->Add( XENON_SD );
  mg->Add( PICO_SD );
  mg->Add( PICASSO_SD );
  mg->Draw("L");
  
  TLegend *leg = new TLegend(0.16,0.61,0.5,0.85,NULL,"brNDC");
  leg->AddEntry("NULL", "Razor-0#mu 90% CL limit: AV EFT operator", "h");
  grshade->SetLineColor(kRed-7);
  grshade->SetLineWidth(3);
  grshade->SetLineStyle(1);
  leg->AddEntry(grshade, "Expected limit, with uncertainty" ,"lf");
  leg->AddEntry(Lc, "Observed limit" ,"l");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->Draw();

  TLegend *leg2 = new TLegend( 0.67, 0.61, 0.92, 0.9, NULL, "brNDC" );
  leg2->AddEntry( ICE_CUBE_WW, "IceCube W^{+}W^{-}", "l" );
  leg2->AddEntry( COUPP_2012_SD, "COUPP 2012", "l" );
  leg2->AddEntry( SIMPLE_2012_SD, "SIMPLE 2012", "l" );
  leg2->AddEntry( SUPER_K_SD, "Super-K W^{+}W^{-}", "l" );
  leg2->AddEntry( XENON_SD, "XENON 100", "l" );
  leg2->AddEntry( PICO_SD, "PICO", "l" );
  leg2->AddEntry( PICASSO_SD, "PICASSO", "l" );
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetLineColor(1);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->Draw();
  
  float lumix = 0.955;
  float lumiy = 0.945;
  float lumifont = 42;
  
  float cmsx = 0.33;
  float cmsy = 0.875;
  TString CMSText = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.078;
  float extray = cmsy - 0.04;
  TString extraText   = "Preliminary";
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TString lumiText = "18.8 fb^{-1} (8 TeV)";
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.75*cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);
  
  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.75*cmsSize);
  //latex.DrawLatex(cmsx-0.12, lumiy, CMSText);
  latex.DrawLatex(cmsx-0.085, cmsy-0.01, CMSText);
  
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  //latex.DrawLatex(extrax, extray, extraText);
  
  std::cout << "here" << std::endl;
  c->cd();
  c->RedrawAxis();
  c->SetLogx();
  c->SetLogy();
  c->Update();
  
  
  c->SaveAs(outputDir+"/Final_av_Lambda_VarCoupling_40Percent.pdf");
  c->SaveAs(outputDir+"/Final_av_Lambda_VarCoupling_40Percent.C");
  
  ///////////////////////
  ////CreatingPLOT V///
  //////////////////////
  delete f3;
  delete Lc, Lcom, Lup, Ldown, grshade, mg0, mg;
  f3 = new TFile("inputData/no_btag/SI_Graphs.root");//Root file containing Razor TGraphs
  Lc = (TGraph*)f3->Get("RazorDM_V_8TeV_Combined_Toys");
  Lcom = (TGraph*)f3->Get("RazorDM_V_8TeV_Combined_Expected");
  
  Lup = (TGraph*)f3->Get("RazorDM_V_8TeV_Combined_Expected_up");
  Ldown = (TGraph*)f3->Get("RazorDM_V_8TeV_Combined_Expected_down");
  
  //Loading direct detection limits
  delete f4;
  //f4 = new TFile("inputData/no_btag/DirectDetection_SI_Lambda_Graphs.root");
  f4 = new TFile("inputData/no_btag/SI_Lambda_Graphs_FINAL.root");
  TGraph* CDMSII_SI = (TGraph*)f4->Get("CDMSII_SI");
  TGraph* CRESTII_NEW = (TGraph*)f4->Get("CRESTII_NEW");
  TGraph* CRESTII_LT = (TGraph*)f4->Get("CRESTII_LowThreshold");
  TGraph* COUPP_2012_SI = (TGraph*)f4->Get("COUPP_2012_SI");
  TGraph* LUX2013 = (TGraph*)f4->Get("LUX2013");
  TGraph* SIMPLE_2012_SI = (TGraph*)f4->Get("SIMPLE_2012_SI");
  TGraph* XENON100_SI = (TGraph*)f4->Get("XENON100_SI");
  TGraph* superCDMS_SI = (TGraph*)f4->Get("superCDMS_SI");
  
  n = Lcom->GetN();
  x = Lup->GetX();
  yup = Lup->GetY();
  ydown = Ldown->GetY();
  grshade = new TGraph(2*n);
  for(int i = 0; i < n; i++) {
    grshade->SetPoint(i,x[i],ydown[i]);
    grshade->SetPoint(n+i,x[n-i-1],yup[n-i-1]);
  }
  grshade->SetFillStyle(1001);
  grshade->SetFillColor(kOrange);
  grshade->SetLineColor(kOrange);
  grshade->SetTitle("");
  
  mg0 = new TMultiGraph();
  mg0->Add(grshade);
  
  mg0->Draw("af");
  mg0->GetYaxis()->SetRangeUser( min_y_si, max_y_si );
  mg0->GetXaxis()->SetRangeUser( 0.001, 1000 );
  mg0->GetXaxis()->SetTitle(labelX);
  mg0->GetYaxis()->SetTitle(labelY);
  mg0->GetYaxis()->SetTitleSize( 0.05 );
  mg0->GetXaxis()->SetTitleSize( 0.05 );
  mg0->GetYaxis()->SetLabelSize( 0.04 );
  mg0->GetXaxis()->SetLabelSize( 0.04 );
  mg0->GetXaxis()->SetTitleOffset(1.1);
  mg0->GetYaxis()->SetTitleOffset(1.05);
  
  Lc->SetLineWidth(3);
  Lc->SetLineStyle(1);
  Lc->SetLineColor(kBlue);
  
  Lcom->SetLineColor(kRed-7);
  Lcom->SetLineWidth(3);
  Lcom->SetLineStyle(1);
  
  Lup->SetLineColor(kRed-7);
  Lup->SetLineWidth(1);
  Lup->SetLineStyle(1);
  
  Ldown->SetLineColor(kRed-7);
  Ldown->SetLineWidth(1);
  Ldown->SetLineStyle(1);

  CDMSII_SI->SetLineColor( kPink-3 );
  CDMSII_SI->SetLineWidth( 3 );
  CDMSII_SI->SetLineStyle( 2 );

  CRESTII_NEW->SetLineColor( kGreen - 3 );
  CRESTII_NEW->SetLineWidth( 3 );
  CRESTII_NEW->SetLineStyle( 2 );

  CRESTII_LT->SetLineColor( kGreen-3 );
  CRESTII_LT->SetLineWidth( 3 );
  CRESTII_LT->SetLineStyle( 3 );
  
  COUPP_2012_SI->SetLineColor( kBlue-7 );
  COUPP_2012_SI->SetLineWidth( 3 );
  COUPP_2012_SI->SetLineStyle( 6 );
  
  LUX2013->SetLineColor( kGreen+3 );
  LUX2013->SetLineWidth( 3 );
  LUX2013->SetLineStyle( 2 );
  
  
  SIMPLE_2012_SI->SetLineColor( kOrange-3 );
  SIMPLE_2012_SI->SetLineWidth( 3 );
  SIMPLE_2012_SI->SetLineStyle( 2 );
    
  XENON100_SI->SetLineColor( kAzure+6 );
  XENON100_SI->SetLineWidth( 3 );
  XENON100_SI->SetLineStyle( 5 ); 

  superCDMS_SI->SetLineColor( kViolet-2 );
  superCDMS_SI->SetLineWidth( 3 );
  superCDMS_SI->SetLineStyle( 6 );
  
  mg = new TMultiGraph();
  mg->Add( Lup );
  mg->Add( Lcom );
  mg->Add( Ldown ); 
  mg->Add( Lc );
  mg->Add( CDMSII_SI );
  mg->Add( CRESTII_NEW );
  mg->Add( CRESTII_LT );
  mg->Add( COUPP_2012_SI );
  mg->Add( LUX2013 );
  mg->Add( SIMPLE_2012_SI );
  mg->Add( XENON100_SI );
  mg->Add( superCDMS_SI );
  mg->Draw("L");
    
  //Legend
  delete leg;
  leg = new TLegend(0.16,0.61,0.5,0.85,NULL,"brNDC");
  leg->AddEntry("NULL", "Razor-0#mu 90% CL limit: V EFT operator", "h");
  grshade->SetLineColor(kRed-7);
  grshade->SetLineWidth(3);
  grshade->SetLineStyle(1);
  leg->AddEntry(grshade, "Expected limit, with uncertainty" ,"lf");
  leg->AddEntry(Lc, "Observed limit" ,"l");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->Draw();

  delete leg2;
  leg2 = new TLegend( 0.67, 0.61, 0.92, 0.9, NULL, "brNDC" );
  leg2->AddEntry( XENON100_SI, "XENON 100", "l" );
  leg2->AddEntry( SIMPLE_2012_SI, "SIMPLE 2012", "l" );
  leg2->AddEntry( COUPP_2012_SI, "COUPP 2012", "l" );
  leg2->AddEntry( superCDMS_SI, "superCDMS", "l" );
  leg2->AddEntry( CDMSII_SI, "CDMSII", "l" );
  leg2->AddEntry( CRESTII_NEW, "CRESSTII", "l" );
  leg2->AddEntry( CRESTII_LT, "CRESSTII LT", "l" );
  leg2->AddEntry( LUX2013, "LUX 2013", "l" );
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetLineColor(1);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->Draw();
  
  
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.75*cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);
  
  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.75*cmsSize);
  latex.DrawLatex(cmsx-0.085, cmsy-0.01, CMSText);
  //latex.DrawLatex(cmsx-0.12, lumiy, CMSText);
  
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  //latex.DrawLatex(extrax, extray, extraText);
  
  c->cd();
  c->RedrawAxis();
  c->SetLogx();
  c->SetLogy();
  c->Update();
  
  c->SaveAs(outputDir+"/Final_v_Lambda_VarCoupling_40Percent.pdf");
  c->SaveAs(outputDir+"/Final_v_Lambda_VarCoupling_40Percent.C");
  
  return;
};
