#include <assert.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2Poly.h>
//LOCAL INCLUDES
#include "AssortedFunctions.hh"
#include "CommandLineInput.hh"

struct binEdges
{
  int xl;
  int yl;
  int xh;
  int yh;
};

struct finalBin
{
  int xl;
  int yl;
  int xh;
  int yh;
  float b_nr;
  float b_smh;
  float b_full;
  float s;
  bool _final;
};

const float minNevents = 10.;

int main( int argc, char* argv[] )
{
  gROOT->Reset();

  //-----------------
  //Selection TString
  //-----------------
  std::string categoryMode = ParseCommandLine( argc, argv, "-category=" );
  if (  categoryMode == "" )
    {
      std::cerr << "[ERROR]: please provide the category. Use --category=<highpt,hzbb,highres,lowres>" << std::endl;
      return -1;
    }

  //-----------------
  //Selection TString
  //-----------------
  std::string lumi = ParseCommandLine( argc, argv, "-lumi=" );
  if (  lumi == "" )
    {
      std::cerr << "[ERROR]: please provide the integrated luminosity. Use --lumi=<lumi(pb-1)>" << std::endl;
      return -1;
    }
  
  TString cut = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40)  && pho1Pt> 25. && pho2Pt>25. && MR > 150. && t1Rsq > 0.0";
  TString categoryCutString;
  TString massWindowCut;
  if (categoryMode == "highpt")
    {
      categoryCutString = " && pTGammaGamma >= 110 ";
      massWindowCut     = " && mGammaGamma > 122.96 && mGammaGamma < 127.04";
    }
  else if (categoryMode == "hzbb")
    {
      categoryCutString = " && pTGammaGamma < 110 && ( abs(mbbH_L-125.) < 15. || ( abs(mbbH_L-125.) >= 15. && abs(mbbZ_L-91.) < 15 ) )";
      massWindowCut     = " && mGammaGamma > 122.5 && mGammaGamma < 127.5";
    }
  else if (categoryMode == "highres")
    {
      categoryCutString = " && pTGammaGamma < 110 && abs(mbbH_L-125.) >= 15 && abs(mbbZ_L-91.) >= 15 && sigmaMoverM < 0.0085";
      massWindowCut     = " && mGammaGamma > 122.98 && mGammaGamma < 127.02";
    }
  else if (categoryMode == "lowres")
    {
      categoryCutString  = " && pTGammaGamma < 110 && abs(mbbH_L-125.) >= 15 && abs(mbbZ_L-91.) >= 15 && sigmaMoverM >= 0.0085 ";
      massWindowCut      = " && mGammaGamma > 120.0 && mGammaGamma < 130.00";
    }
  else
    {
      std::cerr << "[ERROR]: please provide a valid category: <highpt,hzbb,highres, lowres>" << std::endl;
      return -1;
    }

  cut = cut + categoryCutString;
  massWindowCut = cut + categoryCutString + massWindowCut;
  
  TFile* fb = new TFile("/Users/cmorgoth/Work/data/HggRazorRun2/MC/CMSSW_7_6_March15_Ntuples/NonResonant_CMSSW_7_6_March15_1pb_weighted.root");
  assert(fb);
  TTree* bkgTree = (TTree*)fb->Get("HggRazor");
  assert(bkgTree);

  TFile* fsmh = new TFile("/Users/cmorgoth/Work/data/HggRazorRun2/MC/CMSSW_7_6_March15_Ntuples/SystematicsNtuples/SMH_CMSSW_7_6_X_systematics_1pb_weighted.root");
  assert(fsmh);
  TTree* smhTree = (TTree*)fsmh->Get("HggRazor");
  assert(smhTree);

  TFile* fs = new TFile("~/Work/data/HggRazorRun2/MC/CMSSW_7_6_March15_Ntuples/T2bH-Hgg-sbm470-sbw1-chi2m230-chi2w0p1-chi1m100_CMSSW_7_6_March15_1pb_weighted.root");
  assert(fs);
  TTree* sTree = (TTree*)fs->Get("HggRazor");
  assert(sTree);
  
 
  //------------------------------
  //------------------------------
  //Non-resonant bkg histograms
  //------------------------------
  //------------------------------
    
  //---------------
  //full mgg region
  //---------------
  TString bkgCut = "1.37*weight*pileupWeight*"+lumi+"*(" + cut + ")";
  std::cout << "bkgCut: " << bkgCut << std::endl;
  bkgTree->Draw("t1Rsq:MR>>bkgH(197,150,10000, 2000,0, 10)", bkgCut, "goff");
  TH2F* bkgH = (TH2F*)gDirectory->Get("bkgH");
  //-----------------
  //mgg signal region
  //-----------------
  TString bkgCutSR = "1.37*weight*pileupWeight*"+lumi+"*(" + massWindowCut + ")";
  std::cout << "bkgCutSR: " << bkgCutSR << std::endl;
  bkgTree->Draw("t1Rsq:MR>>bkgSRH(197,150,10000, 2000,0, 10)", bkgCutSR, "goff");
  TH2F* bkgSRH = (TH2F*)gDirectory->Get("bkgSRH");
  
  //------------------------------
  //------------------------------
  //SM higgs bkg histograms
  //------------------------------
  //------------------------------
  TString smhCutSR = "weight*pileupWeight*"+lumi+"*(" + massWindowCut + ")";
  std::cout << "smhCutSR: " << smhCutSR << std::endl;
  smhTree->Draw("t1Rsq:MR>>smhSRH(197,150,10000, 2000,0, 10)", smhCutSR, "goff");
  TH2F* smhSRH = (TH2F*)gDirectory->Get("smhSRH");

  //------------------------------
  //------------------------------
  //signal model histogram
  //------------------------------
  //------------------------------
  TString sCutSR = "weight*"+lumi+"*(" + massWindowCut + ")";//SR:Signal Region
   std::cout << "sCutSR: " << sCutSR << std::endl;
  sTree->Draw("t1Rsq:MR>>sSRH(197,150, 10000, 2000, 0, 10)", sCutSR, "goff");
  TH2F* sSRH = (TH2F*)gDirectory->Get("sSRH");
  
  
  TFile* fout = new TFile("test_qnot.root", "RECREATE");
  bkgH->Write("bkgH");
  bkgSRH->Write("bkgSRH");
  smhSRH->Write("smhSRH");
  sSRH->Write("sSRH");
  

  double b_nr[100];//non-resonant bkg
  double b_smh[100];//smh bkg
  double b[100];
  double s[100];//signal is fixed
  double obs[100];
  double bFull[100];

  double bestBkgNR[2];
  double bestBkgSMH[2];
  double bestBkgFull[2];
  double bestS[2];
  
  //------------------------
  //------------------------
  //Binning Algorithm 
  //------------------------
  //------------------------
  //MR
  TH1F* sigmaMR;
  TH1F* fbkg0MR;
  TH1F* fbkg1MR;
  TH1F* bkg0MR;
  TH1F* bkg1MR;
  TH1F* smh0MR;
  TH1F* smh1MR;
  TH1F* s0MR;
  TH1F* s1MR;
  //Rsq
  TH1F* sigmaRsq;
  TH1F* fbkg0Rsq;
  TH1F* fbkg1Rsq;
  TH1F* bkg0Rsq;
  TH1F* bkg1Rsq;
  TH1F* smh0Rsq;
  TH1F* smh1Rsq;
  TH1F* s0Rsq;
  TH1F* s1Rsq;

  bool splitI     = true;
  bool splitII    = true;
  int binMRsplit  = 1;
  int binRSQsplit = 1;
  int nsplits = 2;
  float MR_L = 150.;
  float MR_H = 10000.;
  float Rsq_L = .0;
  float Rsq_H = 10.;
  int nMRbins  = 197;
  int nRSQbins = 2000;

  finalBin bE;
  bE.xl = 1;
  bE.xh = nMRbins;
  bE.yl = 1;
  bE.yh = nRSQbins;
  bE._final = false;
  std::map<int, finalBin> myMap;
  int nfbins = 0;
  myMap[nfbins] = bE;
  float maxSignificance = -99;
  for ( int k = 0; k < 60; k++ )//loop to create bins
    {
      std::cout << "[INFO]: iteration: " << k << std::endl;
           
      int partitionType = -1;//MR = 0, Rsq = 1
      int maxIbin = -1;
      int maxBin = -1;
      std::cout << "==========" << std::endl;
      for( int ibin = 0; ibin <= k; ibin++ )//loop around already created bins
	{
	  //std::cout << "[INFO]: splitting ibin: " << ibin << std::endl;
	  //std::cout << " analyzing bin " << ibin << " :"<< myMap[ibin].xl << "-" << myMap[ibin].xh << "," << myMap[ibin].yl << "-" << myMap[ibin].yh << std::endl;

	  //MR
	  TString index = Form("_%d_%d", k, ibin);
	  sigmaMR = new TH1F("sigmaMR"+index, "sigma-MR", 27, 150, 1500);
	  fbkg0MR = new TH1F("fbkg0MR"+index, "fbkg0-MR", 27, 150, 1500);
	  fbkg1MR = new TH1F("fbkg1MR"+index, "fbkg1-MR", 27, 150, 1500);
	  bkg0MR  = new TH1F("bkg0MR"+index, "bkg0-MR", 27, 150, 1500);
	  bkg1MR  = new TH1F("bkg1MR"+index, "bkg1-MR", 27, 150, 1500);
	  smh0MR  = new TH1F("smh0MR"+index, "smh0-MR", 27, 150, 1500);
	  smh1MR  = new TH1F("smh1MR"+index, "smh1-MR", 27, 150, 1500);
	  s0MR    = new TH1F("s0MR"+index, "s0-MR", 27, 150, 1500);
	  s1MR    = new TH1F("s1MR"+index, "s1-MR", 27, 150, 1500);
	  //Rsq
	  sigmaRsq = new TH1F("sigmaRsq"+index, "sigma-Rsq", 200, 0, 1);
	  fbkg0Rsq = new TH1F("fbkg0Rsq"+index, "fbkg0-Rsq", 200, 0, 1);//full bkg bin0
	  fbkg1Rsq = new TH1F("fbkg1Rsq"+index, "fbkg1-Rsq", 200, 0, 1);//full bkg bin1
	  bkg0Rsq  = new TH1F("bkg0Rsq"+index, "bkg0-Rsq", 200, 0, 1);//sr bkg bin0
	  bkg1Rsq  = new TH1F("bkg1Rsq"+index, "bkg1-Rsq", 200, 0, 1);//sr bkg bin1
	  smh0Rsq  = new TH1F("smh0Rsq"+index, "smh0-Rsq", 200, 0, 1);//sr smh bin0
	  smh1Rsq  = new TH1F("smh1Rsq"+index, "smh1-Rsq", 200, 0, 1);//sr smh bin1
	  s0Rsq    = new TH1F("s0Rsq"+index, "s0-Rsq", 200, 0, 1);//sr signal bin0
	  s1Rsq    = new TH1F("s1Rsq"+index, "s1-Rsq", 200, 0, 1);//sr signal bin1
	  
	  for (int bb = 0; bb < k+2; bb++ )
	    {
	      b[bb]     = 0.0;
	      s[bb]     = 0.0;
	      obs[bb]   = 0;
	      b_nr[bb]  = 0.0;
	      b_smh[bb] = 0.0;
	    }
	  int former_ibin = 0;
	  for ( int ifb = 0; ifb <= k; ifb++ )//add bins that are not going to be partitioned
	    {
	      if ( ifb != ibin )//exclude bin being partitioned
		{
		  b[former_ibin]   = myMap[ifb].b_nr + myMap[ifb].b_smh;
		  s[former_ibin]   = myMap[ifb].s;
		  obs[former_ibin] = b[former_ibin] + s[former_ibin];
		  //std::cout << "adding existing bin info: " << myMap[ifb].xl << "-" << myMap[ifb].xh << "," << myMap[ifb].yl << "-" << myMap[ifb].yh << std::endl;
		  //std::cout << "adding existing bin info: " << b[former_ibin] << " " << s[former_ibin] << std::endl;
		  former_ibin++;
		}
	    }

	  //for (int bb = 0; bb < k+2; bb++ ) std::cout << "Init-->b[" << bb << "] = " << b[bb] << std::endl;
	  if( myMap[ibin]._final )
	    {
	      //std::cout << ibin << " skipping final bin: " << myMap[ibin].xl << "-" << myMap[ibin].xh << "," << myMap[ibin].yl << "-" << myMap[ibin].yh << std::endl; 
	      continue;//do not split bin if already final
	    }
	  
	  float MR  = 175.;
	  for ( int i = 1; i <= 197; i++ )//sum over MR bin ins histogram
	    {
	      if ( i > myMap[ibin].xl && i < myMap[ibin].xh ) //checks that i (MR bin) is above the boundary of the current ibin
		{
		  bFull[k+0] = bkgH->Integral( myMap[ibin].xl, i-1, myMap[ibin].yl, myMap[ibin].yh );
		  bFull[k+1] = bkgH->Integral( i, myMap[ibin].xh, myMap[ibin].yl, myMap[ibin].yh );
		  b_nr[k+0]  = bkgSRH->Integral( myMap[ibin].xl, i-1, myMap[ibin].yl, myMap[ibin].yh );
		  b_nr[k+1]  = bkgSRH->Integral( i, myMap[ibin].xh, myMap[ibin].yl, myMap[ibin].yh );
		  b_smh[k+0] = smhSRH->Integral( myMap[ibin].xl, i-1, myMap[ibin].yl, myMap[ibin].yh );
		  b_smh[k+1] = smhSRH->Integral( i, myMap[ibin].xh, myMap[ibin].yl, myMap[ibin].yh );
		  s[k+0]     = sSRH->Integral( myMap[ibin].xl, i-1, myMap[ibin].yl, myMap[ibin].yh );
		  s[k+1]     = sSRH->Integral( i, myMap[ibin].xh, myMap[ibin].yl, myMap[ibin].yh );
		  b[k+0]     = b_nr[k+0] + b_smh[k+0];
		  b[k+1]     = b_nr[k+1] + b_smh[k+1];
		  obs[k+0]   = b[k+0] + s[k+0];
		  obs[k+1]   = b[k+1] + s[k+1];
		}
	      else
		{
		  continue;
		}
	      if ( b[k+0] == 0  || b[k+1] == 0 ) continue;//avoid nans
	      if ( bFull[k+0] < minNevents || bFull[k+1] <= minNevents ) continue;//avoid bins with small counts
	      
	      double mu   = GetBestFitSignalStrength( k+2, b, s, obs );
	      double qnot = GetQnotTestStatistics( k+2, b, s, obs, mu );
	      /*std::cout << "==========" << std::endl;
	      std::cout << "MR: " << MR << " best fit mu: " << mu << " ==> nsigma = " << sqrt( qnot )
			<< " --> b: " << b[k+0] << ", " << b[k+1] << "; obs: " << obs[k+0] << ", " << obs[k+1] << std::endl;
	      */
	      //-----------------------
	      //filling info histograms
	      //-----------------------
	      sigmaMR->SetBinContent(i, sqrt(qnot) );
	      fbkg0MR->SetBinContent( i, bFull[k+0] );
	      fbkg1MR->SetBinContent( i, bFull[k+1] );
	      bkg0MR->SetBinContent( i, b_nr[k+0] );
	      bkg1MR->SetBinContent( i, b_nr[k+1] );
	      smh0MR->SetBinContent( i, b_smh[k+0] );
	      smh1MR->SetBinContent( i, b_smh[k+1] );
	      s0MR->SetBinContent( i, s[k+0] );
	      s1MR->SetBinContent( i, s[k+1] );
	      MR += 50.; 
	    }

	  if ( sigmaMR->GetMaximum() > 1.001*maxSignificance && sigmaMR->GetMaximum() > 0 )
	    {
	      maxSignificance = sigmaMR->GetMaximum();
	      maxBin = sigmaMR->GetMaximumBin();
	      partitionType = 1;
	      maxIbin = ibin;
	      bestBkgNR[0]   = bkg0MR->GetBinContent(maxBin);
	      bestBkgNR[1]   = bkg1MR->GetBinContent(maxBin);
	      bestBkgSMH[0]  = smh0MR->GetBinContent(maxBin);
	      bestBkgSMH[1]  = smh1MR->GetBinContent(maxBin);
	      bestBkgFull[0] = fbkg0MR->GetBinContent(maxBin);
	      bestBkgFull[1] = fbkg1MR->GetBinContent(maxBin);
	      bestS[0]       = s0MR->GetBinContent(maxBin);
	      bestS[1]       = s1MR->GetBinContent(maxBin);
	      std::cout << "MR: Improved significance: " << maxSignificance << " " << maxBin << " maxIbin " << maxIbin << std::endl;
	      std::cout << "new bins: " << myMap[ibin].xl << "," << maxBin << " - " << maxBin << "," <<  myMap[ibin].xh << std::endl;
	      std::cout << "Bkg0: " << bestBkgFull[0] << "," << bestBkgFull[1] << std::endl;
	      //std::cout << " analyzing bin " << ibin << " :"<< myMap[ibin].xl << "-" << myMap[ibin].xh << "," << myMap[ibin].yl << "-" << myMap[ibin].yh << std::endl;
	      //for (int bb = 0; bb < k+2; bb++ ) std::cout << "b[" << bb << "] = " << b[bb] << std::endl;
	    }
	  
	  float Rsq  = 0.0025;
	  for ( int i = 1; i <= 2000; i++ )
	    {
	      if ( i > myMap[ibin].yl && i < myMap[ibin].yh ) 
		{
		  //std::cout << myMap[ibin].xl << " " <<  myMap[ibin].xh << " " << myMap[ibin].yl << " " << i << std::endl;
		  bFull[k+0] = bkgH->Integral( myMap[ibin].xl, myMap[ibin].xh, myMap[ibin].yl, i-1 );
		  bFull[k+1] = bkgH->Integral( myMap[ibin].xl, myMap[ibin].xh, i, myMap[ibin].yh );
		  b_nr[k+0]  = bkgSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, myMap[ibin].yl, i-1 );
		  b_nr[k+1]  = bkgSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, i, myMap[ibin].yh );
		  b_smh[k+0] = smhSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, myMap[ibin].yl, i-1 );
		  b_smh[k+1] = smhSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, i, myMap[ibin].yh );
		  s[k+0]     = sSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, myMap[ibin].yl, i-1 );
		  s[k+1]     = sSRH->Integral( myMap[ibin].xl, myMap[ibin].xh, i, myMap[ibin].yh );
		  b[k+0]     = b_nr[k+0] + b_smh[k+0];
		  b[k+1]     = b_nr[k+1] + b_smh[k+1];
		  obs[k+0]   = b[k+0] + s[k+0];
		  obs[k+1]   = b[k+1] + s[k+1];
		}
	      else
		{
		  continue;
		}
	       
	      if ( b[k+0] == 0 || b[k+1] == 0 ) continue;
	      if ( bFull[k+0] <= minNevents || bFull[k+1] <= minNevents ) continue;//avoid bins with small counts
	      double mu   = GetBestFitSignalStrength( k+2, b, s, obs );
	      double qnot = GetQnotTestStatistics( k+2, b, s, obs, mu );
	      
	      /*std::cout << "==========" << std::endl;
	      std::cout << " Rsq: " << Rsq << " best fit mu: " << mu << " ==> nsigma = " << sqrt( qnot )
			<< " --> bF: " << bFull[k+0] << ", " << bFull[k+1]
			<< " --> b: " << b[k+0] << ", " << b[k+1]
			<< "; s: " << s[k+0] << ", " << s[k+1]
			<< "; obs: " << obs[k+0] << ", " << obs[k+1] << std::endl;
	      */
	      
	      //-----------------------
	      //filling info histograms
	      //-----------------------
	      sigmaRsq->SetBinContent(i, sqrt(qnot) );
	      fbkg0Rsq->SetBinContent( i, bFull[k+0] );
	      fbkg1Rsq->SetBinContent( i, bFull[k+1] );
	      bkg0Rsq->SetBinContent( i, b_nr[k+0] );
	      bkg1Rsq->SetBinContent( i, b_nr[k+1] );
	      smh0Rsq->SetBinContent( i, b_smh[k+0] );
	      smh1Rsq->SetBinContent( i, b_smh[k+1] );
	      s0Rsq->SetBinContent( i, s[k+0] );
	      s1Rsq->SetBinContent( i, s[k+1] );
	      Rsq += 0.005; 
	    }

	  if ( sigmaRsq->GetMaximum() > 1.001*maxSignificance && sigmaRsq->GetMaximum() > 0)
	    {
	      maxSignificance = sigmaRsq->GetMaximum();
	      maxBin = sigmaRsq->GetMaximumBin();
	      partitionType = 2;
	      maxIbin = ibin;
	      bestBkgNR[0]   = bkg0Rsq->GetBinContent(maxBin);
	      bestBkgNR[1]   = bkg1Rsq->GetBinContent(maxBin);
	      bestBkgSMH[0]  = smh0Rsq->GetBinContent(maxBin);
	      bestBkgSMH[1]  = smh1Rsq->GetBinContent(maxBin);
	      bestBkgFull[0] = fbkg0Rsq->GetBinContent(maxBin);
	      bestBkgFull[1] = fbkg1Rsq->GetBinContent(maxBin);
	      bestS[0]       = s0Rsq->GetBinContent(maxBin);
	      bestS[1]       = s1Rsq->GetBinContent(maxBin);
	      std::cout << "Rsq: Improved significance: " << maxSignificance << " " << maxBin << " maxIbin " << maxIbin << std::endl;
	      std::cout << "new bins: " << myMap[ibin].yl << "," << maxBin << " - " << maxBin << ","<<  myMap[ibin].yh << std::endl;
	      std::cout << "Bkg0: " << bestBkgFull[0] << "," << bestBkgFull[1] << std::endl;
	      
	      //for (int bb = 0; bb < k+2; bb++ ) std::cout << "b[" << bb << "] = " << b[bb] << std::endl;
	    }
	  
	  sigmaMR->Write();
	  fbkg0MR->Write();
	  fbkg1MR->Write();
	  bkg0MR->Write();
	  bkg1MR->Write();
	  smh0MR->Write();
	  smh1MR->Write();
	  s0MR->Write();
	  s1MR->Write();
	  
	  sigmaRsq->Write();
	  fbkg0Rsq->Write();
	  fbkg1Rsq->Write();
	  bkg0Rsq->Write();
	  bkg1Rsq->Write();
	  smh0Rsq->Write();
	  smh1Rsq->Write();
	  s0Rsq->Write();
	  s1Rsq->Write();
	  //MR
	  delete sigmaMR;
	  delete fbkg0MR;
	  delete fbkg1MR;
	  delete bkg0MR;
	  delete bkg1MR;
	  delete smh0MR;
	  delete smh1MR;
	  delete s0MR;
	  delete s1MR;
	  //Rsq
	  delete sigmaRsq;
	  delete fbkg0Rsq;
	  delete fbkg1Rsq;
	  delete bkg0Rsq;
	  delete bkg1Rsq;
	  delete smh0Rsq;
	  delete smh1Rsq;
	  delete s0Rsq;
	  delete s1Rsq;
	}
      
      if ( partitionType == 1 )
	{
	  std::map<int, finalBin> tmpMap;
	  int itmp = 0;
	  for ( auto tmp: myMap )
	    {
	      if ( tmp.first != maxIbin )
		{
		  tmpMap[itmp] = tmp.second;
		  itmp++;
		}
	      else
		{
		  //std::cout << "skippen maxIBin: " << maxIbin << " when cloning map" << std::endl;
		}
	    }
	  std::cout << "splitting in MR @ " << maxBin << ", significance = " << maxSignificance << " nsigmas" << std::endl;
	  std::cout << "MR @ " << maxBin << " NBkg0: " << bestBkgFull[0]
		    << " NBkg1: " << bestBkgFull[1] << std::endl;
	  finalBin fb;
	  //filling low bin
	  fb.xl     = myMap[maxIbin].xl;
	  fb.xh     = maxBin-1;
	  fb.yl     = myMap[maxIbin].yl;
	  fb.yh     = myMap[maxIbin].yh;
	  fb.b_full = bestBkgFull[0];
	  fb.b_nr   = bestBkgNR[0];
	  fb.b_smh  = bestBkgSMH[0];
	  fb.s      = bestS[0];
	  fb._final = false;
	  tmpMap[itmp] = fb;
	  if ( (fb.b_full - minNevents) < minNevents )
	    {
	      std::cout << "MR: First partition is at minimum bkg events: " << fb.b_full << std::endl;
	      tmpMap[itmp]._final = true;
	    }
	  itmp++;
	  //filling high bin
	  fb.xl     = maxBin;
	  fb.xh     = myMap[maxIbin].xh;
	  fb.yl     = myMap[maxIbin].yl;
	  fb.yh     = myMap[maxIbin].yh;
	  fb.b_full = bestBkgFull[1];
	  fb.b_nr   = bestBkgNR[1];
	  fb.b_smh  = bestBkgSMH[1];
	  fb.s      = bestS[1];
	  fb._final = false;
	  tmpMap[itmp] = fb;
	  if ( (fb.b_full - minNevents) < minNevents )
	    {
	      std::cout << "MR: Second partition is at minimum bkg events: " << fb.b_full << std::endl;
	      tmpMap[itmp]._final = true;
	    }
	  myMap = tmpMap;
	}
      else if ( partitionType == 2 )
	{
	  std::cout << "splitting in Rsq @ " << maxBin << ", significance = " << maxSignificance << " nsigmas" << std::endl;
	  std::cout << "Rsq @ " << maxBin << " NBkg0: " << bestBkgFull[0]
		    << " NBkg1: " << bestBkgFull[1] << std::endl;
	  
	  std::map<int, finalBin> tmpMap;
	  int itmp = 0;
	  for ( auto tmp: myMap )
	    {
	      if ( tmp.first != maxIbin )
		{
		  tmpMap[itmp] = tmp.second;
		  itmp++;
		}
	      else
		{
		  //std::cout << "skippen maxIBin: " << maxIbin << " when cloning map" << std::endl;
		}
	    }
	  finalBin fb;
	  //fill low bin
	  fb.xl     = myMap[maxIbin].xl;
	  fb.xh     = myMap[maxIbin].xh;
	  fb.yl     = myMap[maxIbin].yl;
	  fb.yh     = maxBin-1;
	  fb.b_full = bestBkgFull[0];
	  fb.b_nr   = bestBkgNR[0];
	  fb.b_smh  = bestBkgSMH[0];
	  fb.s      = bestS[0];
	  fb._final = false;
	  tmpMap[itmp] = fb;
	  if ( (fb.b_full - minNevents) < minNevents )
	    {
	      std::cout << "Rsq: First partition is at minimum bkg events: " << fb.b_full << std::endl;
	      tmpMap[itmp]._final = true;
	    }
	  itmp++;
	  //filling high bin
	  fb.xl     = myMap[maxIbin].xl;
	  fb.xh     = myMap[maxIbin].xh;
	  fb.yl     = maxBin;
	  fb.yh     = myMap[maxIbin].yh;
	  fb.b_full = bestBkgFull[1];
	  fb.b_nr   = bestBkgNR[1];
	  fb.b_smh  = bestBkgSMH[1];
	  fb.s      = bestS[1];
	  fb._final = false;
	  tmpMap[itmp] = fb;
	  if ( (fb.b_full - minNevents) < minNevents )
	    {
	      std::cout << "Rsq: Second partition is at minimum bkg events: " << fb.b_full << std::endl;
	      tmpMap[itmp]._final = true;
	    }
	  myMap = tmpMap;
	}
      else
	{
	  std::cout << "[INFO]: Significance has not been improved by more than 1%, best significance is: " << maxSignificance << std::endl;
	  std::cout << "[INFO]: Binning process has terminated" << std::endl;
	  break;
	}

      for(int kk = 0; kk < 2; kk++)
	{
	  bestBkgNR[kk]   = -1;
	  bestBkgSMH[kk]  = -1;
	  bestBkgFull[kk] = -1;
	  bestS[kk]       = -1;
	}
    }
  /*
  TRandom3 r1(0);
  TRandom3 r2(0);
  TH1F* h_qnot = new TH1F("h_qnot", "h_qnot", 1000, 0, 20);
  for ( int i = 0; i < 1e7; i ++ )
    {
      for ( int j = 0; j < n; j++ )
	{
	  b[j] = b_nr[j] + b_smh[j];
	  double bNonRes = r1.PoissonD( b_nr[j] );
	  double bSMH = r2.PoissonD( b_smh[j] );
	  obs[j] = bNonRes + bSMH;
	}
      double mu   = GetBestFitSignalStrength( n, b, s, obs );
      double qnot = GetQnotTestStatistics( n, b, s, obs, mu );
      h_qnot->Fill( qnot );
    }
  */

  TH2Poly* h2p = new TH2Poly();
  TH2Poly* h2p_nr = new TH2Poly();
  TH2Poly* h2p_smh = new TH2Poly();
  TH2Poly* h2p_s = new TH2Poly();
  TH2Poly* h2p_SoverSqrtB = new TH2Poly();
 
  int ctr = 0;
  for ( auto tmp : myMap )
    {
      /*std::cout << "BIN #" << tmp.first << "--> xl: " << tmp.second.xl << ", xh: " << tmp.second.xh
		<< ", yl: " << tmp.second.yl << ", yh: " << tmp.second.yh
		<< "; b: " << tmp.second.b_nr + tmp.second.b_smh 
		<< ", s: " << tmp.second.s << std::endl;
      */
      double MR_bw  = (MR_H-MR_L)/nMRbins;
      double Rsq_bw = (Rsq_H-Rsq_L)/nRSQbins;
      
      double MR_low   = 150 + MR_bw*(tmp.second.xl - 1);
      double MR_high  = 150 + MR_bw*((tmp.second.xh+1) - 1);//had to artificially subtract 1 from the xh in the code above, restoring the correct values here
      double Rsq_low  = Rsq_bw*(tmp.second.yl - 1);
      double Rsq_high = Rsq_bw*((tmp.second.yh+1) - 1);//had to artificially subtract 1 from the xh in the code above, restoring the correct values here
      if ( Rsq_high > 1.0 ) Rsq_high  = 1.0;
      if ( MR_high > 9000.0 ) MR_high = 10000.0;

      /*std::cout << "BIN #" << tmp.first << "--> xl: " << MR_low << ", xh: " << MR_high
		<< ", yl: " << Rsq_low << ", yh: " << Rsq_high
		<< "; b: " << tmp.second.b_nr + tmp.second.b_smh 
		<< ", s: " << tmp.second.s << std::endl;
      */

     
      std::cout << "float bin_" << categoryMode << ctr << "[4] = {" << MR_low<< "," << Rsq_low << "," <<  MR_high << "," <<  Rsq_high << "};\n";
      ctr++;
      
      h2p->AddBin( MR_low, Rsq_low, MR_high, Rsq_high );
      h2p_nr->AddBin( MR_low, Rsq_low, MR_high, Rsq_high );
      h2p_smh->AddBin( MR_low, Rsq_low, MR_high, Rsq_high );
      h2p_s->AddBin( MR_low, Rsq_low, MR_high, Rsq_high );
      h2p_SoverSqrtB->AddBin( MR_low, Rsq_low, MR_high, Rsq_high );
      int bin = h2p->FindBin( MR_low+0.5*MR_bw, Rsq_low+0.5*Rsq_bw );
      h2p->SetBinContent( bin, tmp.second.b_full);
      h2p_nr->SetBinContent( bin, tmp.second.b_nr);
      h2p_smh->SetBinContent( bin, tmp.second.b_smh);
      h2p_s->SetBinContent( bin, tmp.second.s);
      h2p_SoverSqrtB->SetBinContent( bin, tmp.second.s/sqrt(tmp.second.b_nr+tmp.second.s) );
    }

  std::cout << "std::vector<float*> SetBinning_" << categoryMode << "()\n{\n";
  std::cout << "std::vector<float*> myVec;\n";
  for( int i = 0; i < myMap.size(); i++ )
    {
      std::cout << "myVec.push_back(bin_" << categoryMode << i << ");\n";
    }
  std::cout << "return myVec;\n}\n";


  std::cout << "\\begin{table*}[htbH]\n\\begin{center}\n";
  std::cout << "\\topcaption{" << categoryMode << " category binning. The non-resonant, "
	    << "SM Higgs, and signal () expected yields for an integrated luminosity correspondint to 2.3~$\\mathrm{fb}^{-1}$ \\label{tab:binning-" << categoryMode << "}}\n";
  std::cout << "\\def\\arraystretch{1.5}\n\\begin{tabular}{|c|c|c|c|}\n\\hline\n";
  std::cout << "$\\mathrm{M_{R}}\\otimes\\mathrm{R^{2}}$ & non-resonant (MC) & SM Higgs & Signal\\\\\n\\hline\n";
  for( auto tmp : myMap )
    {
      double MR_bw  = (MR_H-MR_L)/nMRbins;
      double Rsq_bw = (Rsq_H-Rsq_L)/nRSQbins;
      
      double MR_low   = 150 + MR_bw*(tmp.second.xl - 1);
      double MR_high  = 150 + MR_bw*((tmp.second.xh+1) - 1);//had to artificially subtract 1 from the xh in the code above, restoring the correct values here
      double Rsq_low  = Rsq_bw*(tmp.second.yl - 1);
      double Rsq_high = Rsq_bw*((tmp.second.yh+1) - 1);//had to artificially subtract 1 from the xh in the code above, restoring the correct values here
      if ( Rsq_high > 1.0 ) Rsq_high  = 1.0;
      if ( MR_high > 9000.0 ) MR_high = 10000.0;
      
      TString s = Form("%.0f-%.0f\\GeV $\\otimes$ %.3f-%.3f & %.2f & %.2f & %.2f \\\\", MR_low, MR_high, Rsq_low, Rsq_high, tmp.second.b_nr, tmp.second.b_smh, tmp.second.s);
      std::cout << s << std::endl;
      /*std::cout << MR_low << "-" << MR_high << "\\GeV $\\otimes$ "
		<< Rsq_low << "-"  << Rsq_high;
      std::cout << " & " << tmp.second.b_nr << " & "
		<< tmp.second.b_smh << " & "  << tmp.second.s << "\\\\" << std::endl;
      */
    }
  std::cout << "\\hline\n\\end{tabular}\n\\end{center}\n\\end{table*}" << std::endl;
  
  //h_qnot->Write();
  h2p->Write("h2p_bFull");
  h2p_nr->Write("h2pSR_bNR");
  h2p_smh->Write("h2pSR_bSMH");
  h2p_s->Write("h2pSR_Signal");
  h2p_SoverSqrtB->Write("h2pSR_SoverSqrtB");
  fout->Close();
  
  return 0;
}
