#ifndef DefinePdfs_HH
#define DefinePdfs_HH
//C++ INCLUDES
//ROOT INCLUDES
#include <TTree.h>
#include <TString.h>
//ROOFIT INCLUDES
#include <RooRealVar.h>
#include <RooWorkspace.h>
//LOCAL INCLUDES

TString MakeSingleGauss( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeSingleGaussNE( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeDoubleGauss( TString tag, RooRealVar& mgg, RooWorkspace& w );//mu is fixed to be the same in both Gaussians
TString MakeDoubleGaussNE( TString tag, RooRealVar& mgg, RooWorkspace& w, bool _globalScale = false );//mu is fixed to be the same in both Gaussians
TString MakeFullDoubleGauss( TString tag, RooRealVar& mgg, RooWorkspace& w );//all parameters are floated
TString MakeFullDoubleGaussNE( TString tag, RooRealVar& mgg, RooWorkspace& w, bool _globalScale = false, bool _categoryScale = false, TString category = "category" );//all parameters are floated NonExtended
TString MakeDoubleCB( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeFullTripleGauss( TString tag, RooRealVar& mgg, RooWorkspace& w );//all parameters are floated
TString MakeDoubleExp( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeDoubleExpN1N2( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeSingleExp( TString tag, RooRealVar& mgg, RooWorkspace& w );
TString MakeSingleExpNE( TString tag, RooRealVar& mgg, RooWorkspace& w );//nonExtended
TString MakeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakeDoublePowN1N2(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakeSinglePowNE(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakePoly2(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakePoly3(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString MakePoly4(TString tag, RooRealVar& mgg,RooWorkspace& w);

#endif
