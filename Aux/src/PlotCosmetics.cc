//C++ INCLUDES
#include <fstream>
//ROOT INCLUDES
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
//LOCAL INCLUDES
#include "PlotCosmetics.hh"

const float lumi = 5;
//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.84;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
TString extraText   = "Preliminary";
TString lumiText = "1.26 fb^{-1} (13 TeV)";
//TString lumiText = "19.8 fb^{-1} (8 TeV)";

bool MakeCustomMrRsq( TH2F* h, TString outName )
{
  TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLeftMargin( leftMargin );
  c->SetRightMargin( rightMargin );
  c->SetTopMargin( topMargin );
  c->SetBottomMargin( bottomMargin );
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);

  gStyle->SetPaintTextFormat("4.3f");
  
  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetTitleSize( axisTitleSize );
  h->GetXaxis()->SetTitleOffset( axisTitleOffset );
  h->GetYaxis()->SetTitleSize( axisTitleSize );
  h->GetYaxis()->SetTitleOffset( axisTitleOffset );
  h->SetXTitle("M_{R} (GeV)");
  h->SetYTitle("R^{2} (GeV)");
  h->SetMarkerSize( 1.5 );
  h->Draw("colzTEXTE");
  c->SetLogy();
  c->SetLogx();
  c->SaveAs( outName+".pdf" );
  c->SaveAs( outName+".C" );
  return true;
};

bool MakeStackPlot( THStack* s, TString var, TString outName, TLegend* leg )
{
  TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLeftMargin( leftMargin );
  c->SetRightMargin( rightMargin );
  c->SetTopMargin( topMargin );
  c->SetBottomMargin( bottomMargin );
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  
  s->SetTitle("");
  s->Draw();
  s->GetXaxis()->SetTitleSize( axisTitleSize );
  s->GetXaxis()->SetTitleOffset( axisTitleOffset );
  s->GetYaxis()->SetTitleSize( axisTitleSize );
  s->GetYaxis()->SetTitleOffset( axisTitleOffset );
  if( var == "rsq" )
    {
      s->GetXaxis()->SetRangeUser(0.0, 1.0);
      s->GetXaxis()->SetTitle("R^{2}");
      s->SetMinimum( 1e-1 );
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->GetYaxis()->SetTitle("events / 0.04");
      c->SetLogy();
      c->Update();
    }
  else if ( var == "mr" )
    {
      s->GetXaxis()->SetRangeUser(0.0, 2000.0);
      s->GetXaxis()->SetTitle("M_{R} (GeV)");
      s->GetYaxis()->SetTitle("events / 50 (GeV)");
      s->SetMinimum( 1e-1 );
      s->SetMaximum( 1e1*s->GetMaximum() );
      c->SetLogy();
      c->Update();
    }
  else if ( var == "mgg" )
    {
      s->GetXaxis()->SetRangeUser( 50., 150. );
      s->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
      s->GetYaxis()->SetTitle("events / 1.5 (GeV)");
      s->SetMinimum( 0 );
      s->SetMaximum( 2*s->GetMaximum() );
    }
  else if ( var == "ptgg" )
    {
      //s->GetXaxis()->SetRangeUser( 103., 160. );
      s->GetXaxis()->SetTitle("p^{T}_{#gamma#gamma} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      c->SetLogy();
      c->Update();
    }
    
  if ( leg != NULL )
    {
      leg->SetBorderSize(0);
      //leg->SetTextSize(0.03);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->Draw();
    }
  AddCMS( c );
  c->SaveAs( outName+".pdf" );
  c->SaveAs( outName+".C" );
  delete c;
  return true;
};

bool MakeStackPlot( THStack* s, TH1D* data, TH1D* mc, TString var, TString outName, TLegend* leg )
{
  TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLeftMargin( leftMargin );
  c->SetRightMargin( rightMargin );
  c->SetTopMargin( topMargin );
  c->SetBottomMargin( bottomMargin );
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);

  TPad *pad1 = new TPad("pad1","pad1", .0, 0.25, 1., 1.);
  pad1->SetBottomMargin(0);
  pad1->SetRightMargin( rightMargin );
  pad1->SetLeftMargin( leftMargin );
  pad1->Draw();
  pad1->cd();
  
  s->SetTitle("");
  s->Draw();
  s->GetXaxis()->SetTitleSize( axisTitleSize );
  s->GetXaxis()->SetTitleOffset( axisTitleOffset );
  s->GetYaxis()->SetTitleSize( axisTitleSize );
  s->GetYaxis()->SetTitleOffset( axisTitleOffset );
    
  if ( leg != NULL )
    {
      leg->SetBorderSize(0);
      //leg->SetTextSize(0.03);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->Draw();
    }
  
  data->Draw("sameE");
  AddCMS( c );

  TPad *pad2 = new TPad("pad2","pad2", .0, .0, 1., 0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin( rightMargin );
  pad2->SetLeftMargin( leftMargin );
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  std::cout << "deb" << std::endl;
  //TH1D* MC = (TH1D*)s->GetHistogram();
  TH1D* ratio = new TH1D( *data );
  ratio->Divide( mc );
  ratio->SetMarkerStyle( 20 );
  ratio->GetXaxis()->SetTitleSize( axisTitleSizeRatioX );
  ratio->GetXaxis()->SetLabelSize( axisLabelSizeRatioX );
  ratio->GetXaxis()->SetTitleOffset( axisTitleOffsetRatioX );
  ratio->GetYaxis()->SetTitleSize( axisTitleSizeRatioY );
  ratio->GetYaxis()->SetLabelSize( axisLabelSizeRatioY );
  ratio->GetYaxis()->SetTitleOffset( axisTitleOffsetRatioY );  
  ratio->SetMarkerColor( kBlue );
  ratio->SetLineColor( kBlue );
  ratio->GetYaxis()->SetRangeUser( 0.5, 1.5 );
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("data / mc");
  ratio->GetYaxis()->CenterTitle( true );
  //ratio->GetYaxis()->SetNdivisions( 6, false );
  ratio->SetStats( 0 );
  ratio->Draw("E");

  if( var == "rsq" )
    {
      s->GetXaxis()->SetRangeUser(0.0, 1.0);
      ratio->GetXaxis()->SetRangeUser(0.0, 1.0);
      ratio->GetXaxis()->SetTitle("R^{2}");
      s->SetMinimum( 1e-1 );
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->GetYaxis()->SetTitle("events / 0.04");
      pad1->SetLogy();
      pad1->Update();
      pad2->Update();
    }
  else if ( var == "mr" )
    {
      s->GetXaxis()->SetRangeUser(0.0, 2000.0);
      ratio->GetXaxis()->SetRangeUser(0.0, 2000.0);
      ratio->GetXaxis()->SetTitle("M_{R} (GeV)");
      s->GetYaxis()->SetTitle("events / 50 (GeV)");
      s->SetMinimum( 1e-1 );
      s->SetMaximum( 1e1*s->GetMaximum() );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "mgg" )
    {
      s->GetXaxis()->SetRangeUser( 50., 160. );
      ratio->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
      s->GetYaxis()->SetTitle("events / 1.0 (GeV)");
      s->SetMinimum( 0 );
      s->SetMaximum( 1.2*s->GetMaximum() );
    }
  else if ( var == "ptgg" )
    {
      ratio->GetXaxis()->SetTitle("p_{T}^{#gamma#gamma} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1pt" )
    {
      ratio->GetXaxis()->SetTitle("p^{#gamma_{1}}_{T} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1eta" )
    {
      ratio->GetXaxis()->SetTitle("|#eta|^{#gamma^{1}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 2.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "pho1phi" )
    {
      ratio->GetXaxis()->SetTitle("#phi^{#gamma^{1}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "pho1sigmaIetaIeta" )
    {
      ratio->GetXaxis()->SetTitle("#sigma_{i#etai#eta}^{#gamma^{1}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1r9" )
    {
      ratio->GetXaxis()->SetTitle("R^{#gamma^{1}}_{9}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1HoverE" )
    {
      ratio->GetXaxis()->SetTitle("H/E^{#gamma^{1}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1sumChargedHadronPt" )
    {
      ratio->GetXaxis()->SetTitle("sumChargedHadronPt^{#gamma^{1}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1sumNeutralHadronEt" )
    {
      ratio->GetXaxis()->SetTitle("sumNeutralHadronEt^{#gamma^{1}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1sumPhotonEt" )
    {
      ratio->GetXaxis()->SetTitle("sumPhotonEt^{#gamma^{1}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho1sigmaEoverE" )
    {
      ratio->GetXaxis()->SetTitle("(#sigma_{E}/E)^{#gamma_{1}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "pho2pt" )
    {
      ratio->GetXaxis()->SetTitle("p^{#gamma_{2}}_{T} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2eta" )
    {
      ratio->GetXaxis()->SetTitle("|#eta|^{#gamma^{2}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 2.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "pho2phi" )
    {
      ratio->GetXaxis()->SetTitle("#phi^{#gamma^{2}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "pho2sigmaIetaIeta" )
    {
      ratio->GetXaxis()->SetTitle("#sigma_{i#etai#eta}^{#gamma^{2}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2r9" )
    {
      ratio->GetXaxis()->SetTitle("R^{#gamma^{2}}_{9}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2HoverE" )
    {
      ratio->GetXaxis()->SetTitle("H/E^{#gamma^{2}}");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2sumChargedHadronPt" )
    {
      ratio->GetXaxis()->SetTitle("sumChargedHadronPt^{#gamma^{2}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2sumNeutralHadronEt" )
    {
      ratio->GetXaxis()->SetTitle("sumNeutralHadronEt^{#gamma^{2}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2sumPhotonEt" )
    {
      ratio->GetXaxis()->SetTitle("sumPhotonEt^{#gamma^{2}} (GeV)");
      s->GetYaxis()->SetTitle("events / 10 (GeV)");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }
  else if ( var == "pho2sigmaEoverE" )
    {
      ratio->GetXaxis()->SetTitle("(#sigma_{E}/E)^{#gamma_{2}}");
      s->GetYaxis()->SetTitle("events");
      s->SetMaximum( 1.3e0*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
    }
  else if ( var == "njets" )
    {
      ratio->GetXaxis()->SetTitle("N_{jets}");
      s->GetYaxis()->SetTitle("events");
      s->SetMaximum( 1e1*s->GetMaximum() );
      s->SetMinimum( 1e-1 );
      pad1->SetLogy();
      pad1->Update();
    }

  pad1->SetLogy();
  pad1->Update();
  c->SaveAs( outName+".pdf" );
  c->SaveAs( outName+".C" );
  delete c;
  return true;
};

TH1D* GetStyledHisto( TH1D* h, Process process )
{
  TH1D* h_n = new TH1D( *h );
  SetHistoStyle( h_n, process );
  return h_n;
};

bool SetHistoStyle( TH1D* h, Process process )
{
  if ( process == Process::gammaJet )
    {
      h->SetFillColor( kAzure + 7 );
      h->SetLineColor( kAzure + 7 );
    }
  else if ( process == Process::diphoton )
    {
      h->SetFillColor( kAzure - 9 );
      h->SetLineColor( kAzure - 9 );
    }
  else if ( process == Process::qcd )
    {
      h->SetFillColor( kOrange - 9 );
      h->SetLineColor( kOrange - 9 );
    }
  else if ( process == Process::ttH )
    {
      h->SetFillColor( kViolet + 6 );
      h->SetLineColor( kViolet + 6 );
    }
  else if ( process == Process::ggH )
    {
      h->SetFillColor( kGreen + 1 );
      h->SetLineColor( kGreen + 1 );
    }
  else if ( process == Process::vH )
    {
      h->SetFillColor( kOrange + 7 );
      h->SetLineColor( kOrange + 7 );
    }
  else if ( process == Process::vbfH )
    {
      h->SetFillColor( kOrange - 2 );
      h->SetLineColor( kOrange - 2 );
    }
  else if ( process == Process::w )
    {
      h->SetFillColor( kGreen - 2 );
      h->SetLineColor( kGreen - 2 );
    }
  else if ( process == Process::znunu )
    {
      h->SetFillColor( kAzure - 2 );
      h->SetLineColor( kAzure - 2 );
    }
  else if ( process == Process::dy )
    {
      h->SetFillColor( kViolet - 2 );
      h->SetLineColor( kViolet - 2 );
    }
  else if ( process == Process::tt )
    {
      h->SetFillColor( kPink - 2 );
      h->SetLineColor( kPink - 2 );
    }
  else if ( process == Process::top )
    {
      h->SetFillColor( kBlue - 2 );
      h->SetLineColor( kBlue - 2 );
    }
  else if ( process == Process::vv )
    {
      h->SetFillColor( kRed - 2 );
      h->SetLineColor( kRed - 2 );
    }
  else if ( process == Process::ww )
    {
      h->SetFillColor( kRed - 2 );
      h->SetLineColor( kRed - 2 );
    }
  else if ( process == Process::zz )
    {
      h->SetFillColor( kRed + 4 );
      h->SetLineColor( kRed + 4 );
    }
  else if ( process == Process::wz )
    {
      h->SetFillColor( kRed + 6 );
      h->SetLineColor( kRed + 6 );
    }
  else if ( process == Process::data )
    {
      h->SetMarkerStyle( 20 );
      h->SetMarkerColor( kBlack );
      h->SetLineColor( kBlack );
    }
  else
    {
      std::cout << "[WARNING]: Setting Style, process not found" << std::endl;
      std::cout << "[WARNING]: Setting Fill Color to kGray" << std::endl;
      h->SetFillColor( kGray );
    }
  return true;
};


bool AddLegend( TH1D* h, TLegend* leg, Process process )
{
  if ( leg == NULL ) return false;

  if ( process == Process::gammaJet )
    {
      //leg->AddEntry( h, "#gamma + jets", "f" );
      //leg->AddEntry( h, "Z/#gamma(ll) + jets", "f" );
      leg->AddEntry( h, "prompt-fake", "f" );
      return true;
    }
  else if ( process == Process::diphoton )
    {
      leg->AddEntry( h, "#gamma#gamma + jets", "f" );
      return true;
    }
  else if ( process == Process::qcd )
    {
      leg->AddEntry( h, "fake-fake", "f" );
      return true;
    }
  else if ( process == Process::ttH )
    {
      leg->AddEntry( h, "t#bar{t}H + jets", "f" );
      return true;
    }
  else if ( process == Process::ggH )
    {
      leg->AddEntry( h, "ggH + jets", "f" );
      return true;
    }
  else if ( process == Process::vH )
    {
      leg->AddEntry( h, "VH + jets", "f" );
      return true;
    }
  else if ( process == Process::vbfH )
    {
      leg->AddEntry( h, "VBFH + jets", "f" );
      return true;
    }
   else if ( process == Process::w )
    {
      leg->AddEntry( h, "W(l#nu) + jets", "f" );
      return true;
    }
  else if ( process == Process::znunu )
    {
      leg->AddEntry( h, "Z(#nu#nu) + jets", "f" );
      return true;
    }
  else if ( process == Process::dy )
    {
      leg->AddEntry( h, "Z/#gamma^{*} + jets", "f" );
      return true;
    }
  else if ( process == Process::tt )
    {
      leg->AddEntry( h, "t#bar{t} + jets", "f" );
      return true;
    }
  else if ( process == Process::top )
    {
      leg->AddEntry( h, "t + jets", "f" );
      return true;
    }
  else if ( process == Process::vv )
    {
      leg->AddEntry( h, "vv + jets", "f" );
      return true;
    }
  else if ( process == Process::ww )
    {
      leg->AddEntry( h, "WW + jets", "f" );
      return true;
    }
  else if ( process == Process::zz )
    {
      leg->AddEntry( h, "ZZ + jets", "f" );
      return true;
    }
  else if ( process == Process::wz )
    {
      leg->AddEntry( h, "WZ + jets", "f" );
      return true;
    }
  else if ( process == Process::data )
    {
      leg->AddEntry( h, "data", "lep" );
      return true;
    }
  else
    {
      std::cout << "[WARNING]: unrecognized process, setting unknown legend" << std::endl;
      leg->AddEntry( h, "unknown", "f" );
      return false;
    }
  
  return true;
};

bool AddCMS( TCanvas* C )
{
  C->cd();
  float lumix = 0.955;
  float lumiy = 0.945;
  float lumifont = 42;
  
  float cmsx = 0.33;
  float cmsy = 0.875;
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.078;
  float extray = cmsy - 0.04;
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(cmsx, cmsy, CMSText);
   
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  latex.DrawLatex(extrax, extray, extraText);
  return true;
};
