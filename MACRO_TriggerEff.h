
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TFitResult.h"
#include "TMatrixD.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>

#include <vector>
#include "math.h"

#include <TMultiGraph.h>

#ifndef MACRO_TriggerEff_H
#define MACRO_TriggerEff_H


TF1* fitGamma(TGraphAsymmErrors* g, double q, double p);
double  fitf(double *v, double *par);
double getIntdfdp3(const int fi);

#endif
