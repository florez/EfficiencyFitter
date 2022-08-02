//////////////////////////////////////////
//
//       Andrés Flórez, July 2022       
//
//////////////////////////////////////////
#include "TFile.h"
#include "MACRO_TriggerEff.h"
using namespace std;

void setTDRStyle();

void BACKUP_MACRO_TriggerEff(){
  
  cout << "working" << endl;


  TCanvas* c1 = new TCanvas("c1","c1",200,200,800,700);
  c1->SetGrid();

  // ------------> top pad of the plot
   TPad *top = new TPad("top", "top",0,0.25,1,1);
   top->SetTickx();
   top->SetTicky();
   top->Draw();
   top->cd();
   top->Range(61.49425,-1.836137,167.8161,16.52523);
   top->SetFillColor(0);
   top->SetBorderMode(-1);
   top->SetBorderSize(5);
   top->SetLeftMargin(0.08);
   top->SetRightMargin(0.05);
   top->SetTopMargin(0.05);
   top->SetFrameBorderMode(0);
   top->SetFrameBorderMode(0);

  // Open the root file
  TFile *f1 = TFile::Open("SingleMuon_Run2016_TauPOG_Nov15.root","READ");
  if (f1 == 0) {
    printf("Error: cannot open Trig Eff root file \n");
    return;
  }
  // Read the histograms with the information for the numberator and denominator
  TH1F *Numerator = (TH1F*)f1->Get("NRecoTriggers2/TauJet1Pt");
  TH1F *Denominator = (TH1F*)f1->Get("NMuon1Tau1Combinations/TauJet1Pt");

  // Get the number of bins
  Int_t nbins   = Denominator->GetXaxis()->GetNbins();

  //  Most times we use variable binning. This array is to define the bin width for the different bins.
  int nb = 26;
  float BINS[27] = {0.0, 2.5,  5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0, 60.0, 67.5, 80.0, 100.0, 200.0};

  // Define historgams for your numerator and denominator to be filled with the variable binning defined above
  TH1F *H_bins_num = new TH1F("h_NUM", "h_NUM", nb, BINS);
  TH1F *H_bins_den = new TH1F("h_DEN", "h_DEN", nb, BINS);

  // Fill the histograms
  for (int b = 1; b <= H_bins_den->GetXaxis()->GetNbins(); b++ ){

    int num_bin = 0;
    int den_bin = 0;
    
    for (int i = 1; i < nbins; i++ ){

      if(((( Denominator->GetXaxis()->GetBinLowEdge(i)) + (Denominator->GetXaxis()->GetBinWidth(i))) <= ((H_bins_den->GetXaxis()->GetBinLowEdge(b))+ (H_bins_den->GetXaxis()->GetBinWidth(b)))) &&
	 ((( Denominator->GetXaxis()->GetBinLowEdge(i))+ (Denominator->GetXaxis()->GetBinWidth(i))) > ((H_bins_den->GetXaxis()->GetBinLowEdge(b))))){
	
	num_bin += Numerator->GetBinContent(i);
	den_bin += Denominator->GetBinContent(i);

      }else{
	continue;
      }
      }
      H_bins_num->SetBinContent(b, num_bin);
      H_bins_den->SetBinContent(b, den_bin);

   // }
      
  }

  // Define the TGraphAsummErrors for the efficiency
  TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( H_bins_num, H_bins_den, "b(1,1) mode" );

  // Perform the fit. Below, you will find the fitGamma function. It basically is an error/frequentist function. 
  TF1* gr_fit = fitGamma(gr1, 0, 10);
  TFitResultPtr err = gr1->Fit(gr_fit,"S E");
  
  // Get the covariance and correlations matrices. For now, we are only using the covariance matrix in this example.
  auto corrMatrix = err->GetCorrelationMatrix();
  auto covaMatrix = err->GetCovarianceMatrix();
  // This line are to print the matrices
  corrMatrix.Print();
  covaMatrix.Print();

  // Define the matrix to load the covariant matriz information. This is done automatically by the mnemat function.
  Double_t matrix[5][5];
  gMinuit->mnemat(&matrix[0][0],5);

  //Now, lets draw a pretty error band
  
  const int nn = 19; // made the band for 19 bins since the other effiencies are basically zero.
  // Define the arrys for the error band
  double errxl[nn]={0.0};
  double effplus[nn]={0.0};
  double effminus[nn]={0.0};
  double eff[nn]={0.0};
  double xval[nn]={0.0};
  // I performed the fit first and extracted the 4 fir paramters needed for the band, this can be automatized later. 
  double fit_mean = 40.45;
  double fit_sigma = 0.3259;
  double fit_p0 = 0.0226;
  double fit_p1 = 0.880077;

  double xmean_ratio = 0.;
  // Derivative of fitting function w.r.t px (p0, p1, mean, sigma)
  double dpdx[4] = {0.0}; 

  for (int i=0; i< (nn-1); i++){
    // Get the X axis bins   
    xval[i] = gr1->GetX()[i]; 

    ///////////////////////////////////
    // This is the argument for the Freq function. See the fitf function below 
    double arg = (TMath::Sqrt( gr1->GetX()[i] ) - TMath::Sqrt(fit_mean ))/(2*fit_sigma);
    // This is value of the fit 
    double fitval = fit_p0 + fit_p1*TMath::Freq(arg);
    // Save the fitval in an array to use it later for the band.
    eff[i] = fitval;

    // Derivatives of the fitting function to get the errors
    // The derivatives are w.r.t p0, p1, mean and sigma
    xmean_ratio = gr1->GetY()[i]/fit_mean ;
    dpdx[0] = 1.0;
    dpdx[1] = TMath::Freq(arg);
    dpdx[2] = fit_p1/(4.0*pow(fit_mean, 0.5)*fit_sigma);
    // The derivative has at the end in integral that cannot be solved analytically
    // Therefore, I performed the integral numerically, with a function "getIntdfdp3" defined below.
    dpdx[3] = (fit_p1/(2*fit_sigma))*getIntdfdp3(200);
   
    // Lets get the error 
    double error_v1[4] = {0.0};
    double err_eff = 0.0;
    // Multiply dpdx by the cov matrix
    for (int fl = 0; fl < 4; fl++){
       for (int cl = 0; cl < 4; cl++){
          // First multiply the derivatives by the covariance matrix
          error_v1[fl]+=(dpdx[cl]*matrix[cl][fl]);
       }
      // Now, multiply the result by the derivatives again and add the results. 
      // This gives you the error on the fit square
      err_eff += (error_v1[fl]*dpdx[fl]);
    }
    // take the square root    
    effplus[i]  = pow(err_eff , 0.5)  ;
    effminus[i] = pow(err_eff , 0.5) ;

    // This is to get the last step 
    if (i == (nn-2)){
      eff[nn-1] =  fitval;
      xval[nn-1] = 250.;
      
      effplus[nn-1] =  effplus[i];
      effminus[nn-1] = effminus[i];
    }
  }
  
  TGraphAsymmErrors *gr_fit_e = new TGraphAsymmErrors(nn, xval, eff, errxl, errxl, effminus, effplus );
  
  // Define a pretty graph
  gr1->SetLineColor(kBlack);
  gr1->SetTitle("Trigger Efficiency");
  gr1->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr1->GetXaxis()->SetLabelFont(42);
  gr1->GetXaxis()->SetLabelSize(0.035);
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(0.89);
  gr1->GetXaxis()->SetTitleFont(42);
  gr1->GetYaxis()->SetTitle("#epsilon");
  gr1->GetYaxis()->SetLabelFont(42);
  gr1->GetYaxis()->SetLabelSize(0.035);
  gr1->GetYaxis()->SetTitleSize(0.06);
  gr1->GetYaxis()->SetTitleOffset(0.69);
  gr1->GetYaxis()->SetTitleFont(42);

  gr_fit_e->SetMarkerColor(5);
  gr_fit_e->SetLineColor(3);
  gr_fit_e->SetFillColor(68);
  gr_fit_e->SetMarkerStyle(0);
  gr_fit_e->SetMarkerSize(0);

  gr1->Draw("ap"); 
  gr_fit_e->Draw("same3");
  gr1->Draw("psame");
  //////////////////////////////


  top->Modified();
  c1->cd();

  

// ------------>ratio pad
  TPad *bottom = new TPad("bottom", "",0,0,1,0.32);
  bottom->SetTickx();
  bottom->SetTicky();
  bottom->Draw();
  bottom->cd();
  bottom->SetFillColor(0);
  bottom->SetBorderMode(-1);
  bottom->SetBorderSize(5);
  bottom->SetLeftMargin(0.08);
  bottom->SetRightMargin(0.05);
  bottom->SetTopMargin(0);
  bottom->SetBottomMargin(0.3);
  bottom->SetFrameBorderMode(0);
  bottom->SetFrameBorderMode(0);

  double fit_data_ratio[nn] = {0.0};
  

  for (int i=0; i< nn; i++){

     fit_data_ratio[i] = gr1->GetY()[i]/eff[i];

  }


  TGraphAsymmErrors *gr_ratio_e = new TGraphAsymmErrors(nn, xval, fit_data_ratio, errxl, errxl, effminus, effplus );
  
  gr_ratio_e->Draw("ap");

  TLine *line = new TLine(0.4559321,1.0,273.2287,1.0);
  line->SetLineColor(2);
  line->SetLineStyle(5);
  line->SetLineWidth(3);
  line->Draw();


  bottom->Modified();
  c1->cd();
  c1->Modified();
  c1->cd();

  c1->SetSelected(c1);
  c1->ToggleToolBar();

 
}

double  fitf(double *v, double *par)
{
   double arg = 0;
    if (par[3] != 0) {
      arg = (TMath::Sqrt(v[0]) - TMath::Sqrt(par[2]))/(2*par[3]);
      double fitval =  par[0] + par[1]*TMath::Freq(arg);
      return fitval;
    }
}


TF1* fitGamma(TGraphAsymmErrors* g, double q, double p)
{

  TF1 *fit = new TF1("fit",fitf,0,300,4);
  fit->SetParameters(q, p, g->GetMean(1),g->GetRMS(1));
  fit->SetParNames("P0", "P1", "Mean","Sigma");
  g->Fit("fit","R", 0, 300);
  return fit;
}

double getIntdfdp3(const int fi)
{
   double int_v = 0.0;
   double arg = 0.0;
   double it = 0.0;
   for (int i = 0; i < fi; i++){
     
     arg = (it*it/2.0);
     int_v += arg*TMath::Exp(-arg); 
     it+=1.0;
   }

  return int_v;


}




