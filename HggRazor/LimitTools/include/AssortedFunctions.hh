#include <TH1.h>


double GetEffSigma( TH1* hist );
double GetBestFitSignalStrength( int n, double* b, double* s, double* obs );
double GetQnotTestStatistics( int n, double* b, double* s, double* obs, double mu );
