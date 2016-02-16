#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TDirectory.h>

struct AnaBin
{
  double MR_high;
  double MR_low;
  double Rsq_high;
  double Rsq_low;
  double nevents;
};

void SelectBinning( TString fname, TString categoryMode = "highres" )
{
  TFile* f = new TFile( fname, "READ");
  assert( f );

  TTree* tree = (TTree*)f->Get("HggRazor");
  double MR_max   = 2000.;
  double Rsq_max  = 1.000;
  double MR_step  = 1000.;
  double Rsq_step = 0.500;
  double MR_i     = MR_max - MR_step;
  double Rsq_i    = Rsq_max - Rsq_step;
  
  TString cut = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 25. && pho2Pt>25.";
  TString categoryCutString;
  
  if (categoryMode == "highpt") categoryCutString = " && pTGammaGamma >= 110 ";
  else if (categoryMode == "hbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.)<25";
  else if (categoryMode == "zbb") categoryCutString = " && pTGammaGamma < 110 && abs(mbbZ-91.2)<25 ";
  else if (categoryMode == "highres") categoryCutString = " && pTGammaGamma < 110 && abs(mbbH-125.)>=25 && abs(mbbZ-91.2)>=25 && pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015 ";
  else if (categoryMode == "lowres") categoryCutString = " && pTGammaGamma < 110  && abs(mbbH-125.)>=25 && abs(mbbZ-91.2)>=25 && !(pho1sigmaEOverE < 0.015 && pho2sigmaEOverE < 0.015) ";
  else if (categoryMode == "inclusive") categoryCutString = "";
  
  bool _firstBin = false;
  bool _reachedZero = false;
  int n_stripe_bins = 0;
  std::vector<AnaBin> binVect;
  for ( int i = 0; i < 300; i++ )
    {
      if ( MR_i <= 1000 ) MR_step = 50.0;
      
      TString razorCut = Form("&& MR > %.1f && MR < %.1f && t1Rsq > %.5f && t1Rsq < %.5f", MR_i, MR_max, Rsq_i, Rsq_max );
      TString finalCut = "weight*2300.*1.37*(" + cut + categoryCutString + razorCut + ")";
      tree->Draw("mGammaGamma>>tmp1(38,103,160)", finalCut, "goff");
      TH1F* h = (TH1F*)gDirectory->Get("tmp1");
        
      //check if last Rsq bin can be merged
      bool _ignoreThisBin = false;
      if ( _reachedZero ) 
	{
	  //check if a bin collection is larget than zero
	  if ( binVect.size() > 0 ) 
	    {
	      //check if last bin was not promoted to be bin
	      if ( h->Integral() < binVect.at(binVect.size()-1).nevents )
		{
		  binVect.at(binVect.size()-1).Rsq_low = Rsq_i;
		  binVect.at(binVect.size()-1).nevents = h->Integral() + binVect.at(binVect.size()-1).nevents;
		  _ignoreThisBin = true;
		}
	    }
	}
      
      //adapting Rsq binnning to converge faster
      if( Rsq_i <= 0.5 ) Rsq_step = 0.1;
      if( Rsq_i <= 0.3 ) Rsq_step = 0.05;
      
      //very fine binning to find the first bin in the MR stripe
      if( n_stripe_bins == 0 && h->Integral() > 1.0 )  Rsq_step = 0.002;
      
      if ( h->Integral() < 18. )
	{
	  Rsq_i = Rsq_i - Rsq_step;
	}
      else
	{
	  AnaBin thisbin;
	  thisbin.MR_high = MR_max;
	  thisbin.MR_low = MR_i;
	  thisbin.Rsq_high = Rsq_max;
	  thisbin.Rsq_low  = Rsq_i;
	  thisbin.nevents  = h->Integral();
	  if( !_ignoreThisBin ) binVect.push_back( thisbin );
	  _firstBin = true;
	  Rsq_max = Rsq_i;
	  Rsq_step = 0.05;
	  Rsq_i = Rsq_i - Rsq_step;
	  n_stripe_bins++;
	}
      delete h;
      
      //Check if bottom of Rsq is reached for the first time
      if ( Rsq_i <= 0.0 && !_reachedZero )
	{
	  _reachedZero = true;
	  Rsq_i = 0.0;
	  continue;
	}
      
      //After first MR boundary is reached 
      if ( _firstBin  ) MR_step  = 2.0*MR_i*0.25;

      //Reached the lower left corner, finish binning
      if ( MR_i == 150.0 && Rsq_i < 0.0 ) break;
      
      //If Rsq bottom is reached, do some magic
      if ( Rsq_i < 0.0 )
	{
	  Rsq_max = 1.0; 
	  if ( n_stripe_bins > 0 ) MR_max = MR_i;
	  MR_i = MR_i - MR_step;
	  if ( MR_i <= 150. ) MR_i = 150.;
	  Rsq_step = 0.5;
	  Rsq_i = Rsq_max - Rsq_step;
	  n_stripe_bins = 0;
	  _reachedZero = false;
	}
      
    }

  std::cout << "VECTOR SIZE: " << binVect.size() << std::endl;
  for ( auto& tmp : binVect )
    {
      std::cout << "( " << tmp.MR_low << "-" << tmp.MR_high << ", " << tmp.Rsq_low << "-" << tmp.Rsq_high 
		<< " ) --> " << tmp.nevents << std::endl;
    }

  return;
};
