#ifndef HGG_RAZOR_CLASS_HH
#define HGG_RAZOR_CLASS_HH
//C++ INCLUDES
#include <iostream>
#include <string>
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TH2Poly.h>
//LOCAL INCLUDES
#include "HggTree.hh"
#include "HggAux.hh"

class HggRazorSystematics : public HggTree
{
public:
  HggRazorSystematics( );
  HggRazorSystematics( TTree* tree );
  HggRazorSystematics( TTree* tree, TString processName, TString boxName, bool info = false, bool debug = false );
  ~HggRazorSystematics( );

  bool InitMrRsqTH2Poly( );
  void Loop();
  void PrintBinning();
  float GetYields( float mr, float rsq, float mgg_l, float mgg_h );
  float GetYields( float mr, float rsq, float mgg_l, float mgg_h, double& err );
  bool WriteOutput( TString outName = "default" );

  float GetHighPtGB( double mr, double r2 );
  float GetHighResGB( double mr, double r2 );
  std::pair<float, float> GetFacScaleSystematic( float mr, float rsq );
  std::pair<float, float> GetRenScaleSystematic( float mr, float rsq );
  std::pair<float, float> GetFacRenScaleSystematic( float mr, float rsq );
  static float GetHggBF( );

  bool SetBinningMap( std::map<std::pair<float, float>, std::vector<float>> myMap ){ this->binningMap = myMap; return true;};
  bool SetNeventsHisto( TH1F* histo );
  bool SetFacScaleWeightsHisto( TH1F* histo );
  bool SetPdfWeightsHisto( TH1F* histo );
private:

  static constexpr float hggBF = 2.28e-3;
  static float Lumi;//lumi in pb-1
  
  bool _debug;
  bool _info;
  TString processName;
  TString boxName;

  //----------------------------
  // D e f i n e   B i n n i n g
  //----------------------------
  std::map<std::pair<float, float>, std::vector<float>> binningMap;
  
  //TH2Poly Histos
  TH2Poly *h2p;
  TH2Poly *h2p_facScaleUp;
  TH2Poly *h2p_facScaleDown;
  TH2Poly *h2p_renScaleUp;
  TH2Poly *h2p_renScaleDown;
  TH2Poly *h2p_facRenScaleUp;
  TH2Poly *h2p_facRenScaleDown;

  //Systematic Histos(user must set these)
  TH1F* NEvents;
  TH1F* SumScaleWeights;
  TH1F* SumPdfWeights;
  
  // o u t p u t  r o o t  f i l e
  //------------------------------
  TFile* fout;
};

#endif
