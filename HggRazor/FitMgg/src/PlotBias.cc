//C++
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
//ROOT
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
//RooFit
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooPlot.h>
//LOCAL
#include "PlotBias.hh"
#include "HggAux.hh"


//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.84;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.52;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

void PlotBias( std::string fname = "", std::string outDir = "bias_plots")
{
  //-----------------------
  //G e t t i n g   F i l e
  //-----------------------
  
  TFile* f = new TFile( fname.c_str() , "READ" );
  //Getting Workspace
  RooWorkspace* w = (RooWorkspace*)f->Get("w_biasSignal");
  //Getting bias RooRealVar
  RooRealVar* bias = w->var("bias");
  //Getting data_bias RooDataSet
  RooDataSet* data_bias = (RooDataSet*)w->obj("data_bias");

  //-------------------------------------
  //D e f i n i n g   p l o t   r a n g e
  //-------------------------------------
  double mean = data_bias->mean( *bias );
  double rms  = data_bias->rmsVar( *bias )->getVal();
  if (rms == .0 ) rms = 0.01;
  double bias_min = mean - 5.*rms;
  double bias_max = mean + 5.*rms;
  bias->setMin( bias_min );
  bias->setMax( bias_max );
  bias->setBins(50);

  //---------------------------------
  //S e t t i n g   u p   C a n v a s 
  //---------------------------------
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
  
  RooPlot* frame = bias->frame();
  frame->SetTitle("");
  frame->GetXaxis()->SetTitleSize( axisTitleSize );
  frame->GetXaxis()->SetTitleOffset( axisTitleOffset );
  frame->GetYaxis()->SetTitleSize( axisTitleSize );
  frame->GetYaxis()->SetTitleOffset( axisTitleOffset );
  data_bias->plotOn(frame)->Draw();

  //------------------------------
  //P l ot  M e a n  a n d   R m s 
  //------------------------------
  TString s_mean = Form("mean = %.2f %%", mean*100.);
  TString s_rms;
  if ( mean < .0 )
    {
      s_rms = Form("rms  =  %.2f %%", rms*100.);
    }
  else
    {
      s_rms = Form("rms  = %.2f %%", rms*100.);
    }
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.04);    
  latex.DrawLatex( 0.9, 0.87, s_mean );
  latex.DrawLatex( 0.9, 0.81, s_rms );

  //-----------------------
  //S a v i n g   P l o t s
  //-----------------------
  //std::string outName = fname;
  std::string findName = "biasTest_";
  int begin_s = fname.find( findName ) + findName.size();
  int end_s   = fname.find( ".root" );
  std::string outName = outDir + "/" +fname.substr( begin_s, end_s - begin_s);
  TString figName = outName;
  c->SaveAs(figName+".pdf");
  c->SaveAs(figName+".png");
  c->SaveAs(figName+".C");
  delete c;
  return;
};

//std::pair<double,double> GetMeanRms( std::string fname = "", std::string dataSet = "bias_plots", std::string var = "bias")
std::pair<double,double> GetMeanRms( std::string fname, std::string dataSet, std::string var )
{ 
  //-----------------------
  //G e t t i n g   F i l e
  //-----------------------
  TFile* f = new TFile( fname.c_str() , "READ" );
  //Getting Workspace
  //RooWorkspace* w = (RooWorkspace*)f->Get("w_bias");
  RooWorkspace* w = (RooWorkspace*)f->Get("w_biasSignal");
  //Getting bias RooRealVar
  RooRealVar* bias = w->var( var.c_str() );
  //Getting data_bias RooDataSet
  RooDataSet* data_bias = (RooDataSet*)w->obj( dataSet.c_str() );
  
  //-------------------------------------
  //D e f i n i n g   p l o t   r a n g e
  //-------------------------------------
  double mean = data_bias->meanVar( *bias )->getVal();
  double rms  = data_bias->rmsVar( *bias )->getVal();

  //-----------------------------------
  //G e t t i n g   p a i r   n a m e s 
  //-----------------------------------
  std::string findName = "biasTest_";
  int begin_s = fname.find( findName ) + findName.size();
  int end_s   = fname.find_last_of( "_" );
  std::string f1 = fname.substr( begin_s, end_s - begin_s);
  begin_s = end_s + 1;
  end_s   = fname.find( ".root" );
  std::string f2 = fname.substr( begin_s, end_s - begin_s);
  std::cout << "f1: "<< f1 << ", f2: " << f2 << ", mean = " << mean <<    std::endl;
  f->Close();
  return std::make_pair( mean, rms );
};

void MakeTable( std::map< std::pair<std::string,std::string>, double > mymap, TString type )
{
  int row_ctr = 0;
  std::vector< std::string > row_text;
  std::stringstream ss_fl;//first line for table
  ss_fl << "--";
  for( const auto& fitf : FitFunction() )
    //for( const auto& fitf : mymap )
    {
      std::string f1 = GetFitFunctionString( fitf );
      //std::string f1 = fitf.first.first;
      std::stringstream ss;
      ss << f1; 
      ss.precision(1);
      ss.setf ( std::ios::showpoint );
      ss.setf ( std::ios::fixed ); 
      for( const auto& fitf2 : FitFunction() )
	{
	  std::string f2 = GetFitFunctionString( fitf2 );
	  std::pair<std::string, std::string> mypair = std::make_pair( f1, f2 );
	  if ( mymap.find( mypair ) == mymap.end() )
	    {
	      continue;
	    }
	  if ( row_ctr == 0 ) ss_fl << " & " << f2; 
	  std::pair< std::string, std::string > pair = std::make_pair( f1, f2 );
	  //std::cout << f1 << " " << f2 << " " << mymap[pair] << std::endl;
	  ss << " & " << mymap[pair]*100. << "\\%"; 
	  //std::cout << GetFitFunctionString( fitf )  << " "  << GetFitFunctionString( fitf2 ) << std::endl;
	}
      if ( row_ctr == 0 ) row_text.push_back( ss_fl.str() );
      row_text.push_back( ss.str() );
      row_ctr++;
    }

  std::cout << "\\begin{table}[htb]\n\\caption{Bias test performed in HighRes box. using the set of seven fit functions}\n\\centering\n\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n\\hline" << std::endl;
  int ctr = 0;
  for ( auto tmp : row_text)
    {
      std::cout << tmp << "\\\\" << std::endl;
      if ( ctr == 0 ) std::cout << "\\hline\n";
      ctr++;
    }
  std::cout << "\\hline\n\n\\hline\n\\hline\n\\end{tabular}\n\\end{table}" << std::endl;
   return;
};

double FitBias( TString fname = "", TString f1 = "dumm1", TString f2 = "dummy2" )
{
  TFile* f = new TFile( fname , "READ" );
  TTree* tree = (TTree*)f->Get("biasTree");
  
  //tree->Draw("biasNorm>>hbias(200,-50, 50)", "status==0 && covStatus==3 && status2==0 && covStatus2==3", "goff");
  tree->Draw("biasNorm>>hbias(200,-50, 50)", "", "goff"); 
  TH1F* hbias = (TH1F*)gDirectory->Get("hbias");
  double mean_val = hbias->GetMean();
  
  double _max = hbias->GetBinContent( hbias->GetMaximumBin() );
  //double _xlow = hbias->GetBinCenter( hbias->FindFirstBinAbove( 0.5*_max ) );
  //double _xhigh = hbias->GetBinCenter( hbias->GetMaximumBin() ) + ( hbias->GetBinCenter( hbias->GetMaximumBin() ) - _xlow );

  double _xlow = hbias->GetBinCenter(  hbias->GetMaximumBin() - 2 );
  double _xhigh = hbias->GetBinCenter(  hbias->GetMaximumBin() + 3 );

  double _xlow2sig = hbias->GetBinCenter( hbias->FindFirstBinAbove( 0.03*_max ) );
  double _xhigh2sig = hbias->GetBinCenter( hbias->GetMaximumBin() ) + ( hbias->GetBinCenter( hbias->GetMaximumBin() ) - _xlow2sig );
  
 
  TCanvas *C = new TCanvas("c", "c",2119,33,800,700);
  C->SetHighLightColor(2);
  //C->Range(-0.6543224,-1290.871,3.177829,8696.391);
  C->SetFillColor(0);
  C->SetBorderMode(0);
  C->SetBorderSize(2);
  C->SetLeftMargin(0.16);
  C->SetRightMargin(0.05);
  C->SetTopMargin(0.07);
  C->SetBottomMargin(0.12);
  C->SetFrameBorderMode(0);
  C->SetFrameBorderMode(0);
  
  TF1* gaus = new TF1( "gaus", "gaus(0)", _xlow, _xhigh );
  hbias->Fit( gaus, "R");
  //hbias->GetXaxis()->SetRangeUser( _xlow2sig, _xhigh2sig );
  hbias->GetXaxis()->SetRangeUser( -4, 4 );
  hbias->GetYaxis()->SetRangeUser( 0, 1.2*_max );
  std::cout << "Gaussian Mean: " << gaus->GetParameter(1) << std::endl;
  hbias->SetStats( kFALSE );
  hbias->Draw("e");
  

  //-----------------
  //LATEX
  //-----------------
  float lumix = 0.92;
  float lumiy = 0.85;
  float lumifont = 42;
   
  float cmsx = 0.28;
  float cmsy = 0.875;
  TString CMSText = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.067;
  float extray = cmsy - 0.04;
  TString extraText   = "Preliminary";
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.05;
  
  TLatex latex;
  TString lumiText;
  TString lumiText2;
  lumiText = Form("mean: %.2f", mean_val );
  lumiText2 = Form("#mu: %.2f #pm %.2f", gaus->GetParameter(1), gaus->GetParError(1) );
  
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);    
  latex.DrawLatex(lumix, lumiy, lumiText);
  latex.DrawLatex(lumix, lumiy-0.1, lumiText2);
  
  C->SaveAs( "biasFits/" + f1 + "_" + f2 + ".pdf" );
  return gaus->GetParameter(1);
};
