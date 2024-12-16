const Double_t cutdummy= -99999.0;

#include <TH1D.h>
#include <TRandom3.h>

double GetBinCenterLog(double low, double up){
  //double width = TMath::Power(10,up)-TMath::Power(10, low);
  return TMath::Power(10,(low+up)/2.);///1.075;//TMath::Power(10,low) + width/2.; //0.5*(TMath::Power(10,low)+TMath::Power(10,up));
}

TGraphAsymmErrors* ShiftTickUsingTGraph(TH1D* h, double Shift = 0.2) {
    // Convert TH1D to TGraphAsymmErrors
    int nbins = h->GetNbinsX();
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);

    double abs_shift = TMath::Abs(Shift);

    for (int i = 1; i <= nbins; ++i) {
        double x = h->GetXaxis()->GetBinCenter(i) + Shift; // Apply the shift
        double y = h->GetBinContent(i);
        double exl = h->GetXaxis()->GetBinCenter(i) - h->GetXaxis()->GetBinLowEdge(i) + Shift; // Symmetric x-errors for now (set to 0)
        double exh = h->GetXaxis()->GetBinUpEdge(i) - h->GetXaxis()->GetBinCenter(i) - Shift;
        double eyl = h->GetBinError(i); // y-errors from histogram
        double eyh = h->GetBinError(i);

        g->SetPoint(i - 1, x, y);
        g->SetPointError(i - 1, exl, exh, eyl, eyh);
    }

    // Copy style attributes from the histogram
    g->SetLineColor(h->GetLineColor());
    g->SetLineWidth(h->GetLineWidth());
    g->SetLineStyle(h->GetLineStyle());
    g->SetMarkerStyle(h->GetMarkerStyle());
    g->SetMarkerSize(h->GetMarkerSize());
    g->SetMarkerColor(h->GetMarkerColor());
    g->SetFillColor(h->GetFillColor());
    g->SetFillStyle(h->GetFillStyle());

    return g;
}


TH1D* SmearHistogram(TH1D* hist_or, double sigma=3) {
  TH1D *hist = (TH1D*)hist_or->Clone(TString::Format("%s_smeared", hist_or->GetName()));
    // Create a random generator
    TRandom3 randGen(0); // Seed with 0 for a random seed

    // Loop over all bins of the histogram
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        // Get bin content and bin center
        double binContent = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);

        // Get the center of the bin
        double binCenter = hist->GetBinCenter(i);

        // Apply Gaussian smearing to the bin center
        double smearedValue = randGen.Gaus(binCenter, sigma);

        // Update the histogram
        hist->SetBinContent(i, binContent);  // Optionally modify binContent
        hist->SetBinError(i, binError);     // Keep bin errors the same
    }
    return hist;
}

void DivideBinWidthLog(TH1D * h){
    double err;
    for (int j = 1; j <= h->GetNbinsX();j++){
    err = h->GetBinError(j);
    double width = TMath::Power(10,h->GetXaxis()->GetBinUpEdge(j))-TMath::Power(10, h->GetXaxis()->GetBinLowEdge(j));
    cout << j << " width: " << width << endl;
    //double width = h->GetBinCenter(j)*0.1*TMath::Log(10);
      h->SetBinContent(j,((h->GetBinContent(j))/(width)));
      h->SetBinError(j,((err)/(width)));
    }
}

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

//==============================================================================
// Example Unfolding
//==============================================================================

TH2D* SwapIJ(TH2D* NMatrix){
  TH2D * Matrix = (TH2D*)NMatrix->Clone(TString::Format("%s_swap", NMatrix->GetName()));
  int Dim = NMatrix->GetXaxis()->GetNbins();
  for (int i = 1; i <= Dim ; i++){
    for (int j = 1 ; j <= Dim ; j++) {
          Matrix->SetBinContent(i,j,NMatrix->GetBinContent(j,i));//norm by rows
          //Matrix->SetBinContent(i,j,NMatrix->GetBinContent(i,j)/Sum);//norm by columns
        }
    
  }
  return Matrix;
}


TH2D* NormalizeMatrix(TH2D* NMatrix){
  TH2D * Matrix = (TH2D*)NMatrix->Clone(TString::Format("%s_normalized", NMatrix->GetName()));
  int Dim = NMatrix->GetXaxis()->GetNbins();
  double Sum;
  for (int i = 1; i <= Dim ; i++){
    Sum = 0.0;
    for (int j = 1 ; j <= Dim ; j++) {
    Sum = Sum + NMatrix->GetBinContent(j,i);//norm by rows
    //Sum = Sum + NMatrix->GetBinContent(i,j); //norm by columns
    }
    if (Sum > 0.){
        for (int j = 1 ; j <= Dim ; j++) {
          Matrix->SetBinContent(j,i,NMatrix->GetBinContent(j,i)/Sum);//norm by rows
          //Matrix->SetBinContent(i,j,NMatrix->GetBinContent(i,j)/Sum);//norm by columns
        }
    }
  }
  return Matrix;
}

int RooUnfoldExample()
{
  
  cout << "==================================== TRAIN ====================================" << endl;
  TFile *  rfile = new TFile("astro.root", "read");
  TH2D * matrix = (TH2D*)rfile->Get("astro_matrix_short");
  TH2D * matrix_nooverflow = (TH2D*)rfile->Get("astro_matrix_short_nooverflow");
  TH1D * N = (TH1D*)rfile->Get("N"); //N_overflow
  N->SetName("N_or");
  TH1D * Ncorr = (TH1D*)rfile->Get("N");//Ncorr_overflow
  TH1D * Ntruth = (TH1D*)Ncorr->Clone("Ntruth");
  //Ntruth->Reset();

TH1D * h_frompdf = (TH1D*)N->Clone("h_frompdf"); //https://arxiv.org/pdf/2008.06486 fig.14 blue point
h_frompdf->Reset();


h_frompdf->SetBinContent(1 , 4.329421296327236e+37);
h_frompdf->SetBinContent(2 , 4.0643452499288163e+37);
h_frompdf->SetBinContent(3 , 3.84575171901762e+37);
h_frompdf->SetBinContent(4 , 3.84575171901762e+37);
h_frompdf->SetBinContent(5 , 4.329421296327236e+37);
h_frompdf->SetBinContent(6 , 4.835579728406051e+37);
h_frompdf->SetBinContent(7 , 5.358427269539197e+37);
h_frompdf->SetBinContent(8 , 5.798777047170974e+37);
h_frompdf->SetBinContent(9 , 5.798777047170974e+37);
h_frompdf->SetBinContent(10, 6.0323418750209005e+37);
h_frompdf->SetBinContent(11, 5.574255561017896e+37);
h_frompdf->SetBinContent(12, 5.574255561017896e+37);
h_frompdf->SetBinContent(13, 5.358427269539197e+37);
h_frompdf->SetBinContent(14, 3.4432023310243814e+37);
h_frompdf->SetBinContent(15, 2.894040389083963e+37);
h_frompdf->SetBinContent(16, 1.3137534799428786e+37);
h_frompdf->SetBinContent(17, 7.800815599105648e+36);
h_frompdf->SetBinContent(18, 8.309583735309199e+36);
h_frompdf->SetBinContent(19, 5.297543905247407e+36);
h_frompdf->SetBinContent(20, 8.441878217224814e+36);

h_frompdf->SetBinError(1 , 2.1250309e+36);
h_frompdf->SetBinError(2 , 2.1250309e+36);
h_frompdf->SetBinError(3 , 2.1250309e+36);
h_frompdf->SetBinError(4 , 2.1250309e+36);
h_frompdf->SetBinError(5 , 2.1250309e+36);
h_frompdf->SetBinError(6 , 2.1250309e+36);
h_frompdf->SetBinError(7 , 2.1250309e+36);
h_frompdf->SetBinError(8 , 3.1250309e+36);
h_frompdf->SetBinError(9 , 3.1250309e+36);
h_frompdf->SetBinError(10, 3.1250309e+36);
h_frompdf->SetBinError(11, 3.1250309e+36);
h_frompdf->SetBinError(12, 3.1250309e+36);
h_frompdf->SetBinError(13, 3.1250309e+36);
h_frompdf->SetBinError(14, 3.3834613e+36);
h_frompdf->SetBinError(15, 3.3813675e+36);
h_frompdf->SetBinError(16, 3.5957909e+36);
h_frompdf->SetBinError(17, 4.8279511e+36);
h_frompdf->SetBinError(18, 3.0457440e+36);
h_frompdf->SetBinError(19, 0);
h_frompdf->SetBinError(20, 0);

h_frompdf->SetMarkerColor(kBlue);
h_frompdf->SetLineColor(kBlue);
h_frompdf->SetMarkerStyle(21);
h_frompdf->SetMarkerSize(3);
h_frompdf->SetLineWidth(3);

TH1D * h_frompdf2 = (TH1D*)N->Clone("h_frompdf2"); //https://arxiv.org/pdf/2008.06486 fig.14 red point
h_frompdf2->Reset();


h_frompdf2->SetBinContent(1 , 5.191797106187108e+37);
h_frompdf2->SetBinContent(2 , 4.835579728406051e+37);
h_frompdf2->SetBinContent(3 , 4.722357123022113e+37);
h_frompdf2->SetBinContent(4 , 4.685208472155697e+37);
h_frompdf2->SetBinContent(5 , 5.070233854622015e+37);
h_frompdf2->SetBinContent(6 , 5.400913819499949e+37);
h_frompdf2->SetBinContent(7 , 6.12838104720031e+37);
h_frompdf2->SetBinContent(8 , 6.684589442794821e+37);
h_frompdf2->SetBinContent(9 , 7.349090977461116e+37);
h_frompdf2->SetBinContent(10, 7.584959252938284e+37);
h_frompdf2->SetBinContent(11, 7.06454325046431e+37);
h_frompdf2->SetBinContent(12, 7.953031353047906e+37);
h_frompdf2->SetBinContent(13, 8.079649267497752e+37);
h_frompdf2->SetBinContent(14, 1.0159207110885235e+38);
h_frompdf2->SetBinContent(15, 8.74362613595501e+37);
h_frompdf2->SetBinContent(16, 5.443737241986397e+37);
h_frompdf2->SetBinContent(17, 6.17697250033443e+37);
h_frompdf2->SetBinContent(18, 8.851533660407454e+36);
h_frompdf2->SetBinContent(19, 2.738389609352454e+37);
h_frompdf2->SetBinContent(20, 2.7819867191707236e+37);

h_frompdf2->SetBinError(1 , 5.2041600e+35);
h_frompdf2->SetBinError(2 , 5.2041600e+35);
h_frompdf2->SetBinError(3 , 5.2041600e+35);
h_frompdf2->SetBinError(4 , 5.2041600e+35);
h_frompdf2->SetBinError(5 , 5.2041600e+35);
h_frompdf2->SetBinError(6 , 5.2041600e+35);
h_frompdf2->SetBinError(7 , 5.2041600e+35);
h_frompdf2->SetBinError(8 , 5.2041600e+35);
h_frompdf2->SetBinError(9 , 5.2041600e+35);
h_frompdf2->SetBinError(10, 5.2041600e+35);
h_frompdf2->SetBinError(11, 5.2041600e+36);
h_frompdf2->SetBinError(12, 7.9059478e+36);
h_frompdf2->SetBinError(13, 1.0161552e+37);
h_frompdf2->SetBinError(14, 1.2776954e+37);
h_frompdf2->SetBinError(15, 1.4155810e+37);
h_frompdf2->SetBinError(16, 1.6768202e+37);
h_frompdf2->SetBinError(17, 2.1619922e+37);
h_frompdf2->SetBinError(18, 1.5473118e+37);
h_frompdf2->SetBinError(19, 3.1063655e+37);
h_frompdf2->SetBinError(20, 0);

h_frompdf2->SetMarkerColor(kRed);
h_frompdf2->SetLineColor(kRed);
h_frompdf2->SetMarkerStyle(21);
h_frompdf2->SetMarkerSize(3);
h_frompdf2->SetLineWidth(3);



  TH1D * h_paper = (TH1D*)N->Clone("h_paper");
  h_paper->Reset();
h_paper->SetBinContent(1 , 1.9383*1e-18) ;
h_paper->SetBinContent(2 , 9.076 *1e-19) ;
h_paper->SetBinContent(3 , 4.310 *1e-19) ;
h_paper->SetBinContent(4 , 2.164 *1e-19) ;
h_paper->SetBinContent(5 , 1.227 *1e-19) ;
h_paper->SetBinContent(6 , 6.852 *1e-20) ;
h_paper->SetBinContent(7 , 3.796 *1e-20) ;
h_paper->SetBinContent(8 , 2.055 *1e-20) ;
h_paper->SetBinContent(9 , 1.035 *1e-20) ;
h_paper->SetBinContent(10, 0.533 *1e-20) ;
h_paper->SetBinContent(11, 2.492 *1e-21) ;
h_paper->SetBinContent(12, 1.252 *1e-21) ;
h_paper->SetBinContent(13, 5.98  *1e-22) ;
h_paper->SetBinContent(14, 1.93  *1e-22) ;
h_paper->SetBinContent(15, 8.10  *1e-23) ;
h_paper->SetBinContent(16, 1.86  *1e-23) ;
h_paper->SetBinContent(17, 5.5   *1e-24) ;
h_paper->SetBinContent(18, 2.9   *1e-24) ;
h_paper->SetBinContent(19, 9.5*1e-25);
h_paper->SetBinContent(20, 7.5*1e-25);

h_paper->SetBinError(1 , 0.0067*1e-18) ;
h_paper->SetBinError(2 , 0.042 *1e-19) ;
h_paper->SetBinError(3 , 0.025 *1e-19) ;
h_paper->SetBinError(4 , 0.016 *1e-19) ;
h_paper->SetBinError(5 , 0.011 *1e-19) ;
h_paper->SetBinError(6 , 0.074 *1e-20) ;
h_paper->SetBinError(7 , 0.049 *1e-20) ;
h_paper->SetBinError(8 , 0.032 *1e-20) ;
h_paper->SetBinError(9 , 0.021 *1e-20) ;
h_paper->SetBinError(10, 0.013 *1e-20) ;
h_paper->SetBinError(11, 0.081 *1e-21) ;
h_paper->SetBinError(12, 0.052 *1e-21) ;
h_paper->SetBinError(13, 0.32  *1e-22) ;
h_paper->SetBinError(14, 0.17  *1e-22) ;
h_paper->SetBinError(15, 0.99  *1e-23) ;
h_paper->SetBinError(16, 0.46  *1e-23) ;
h_paper->SetBinError(17, 2.5   *1e-24) ;
h_paper->SetBinError(18, 1.7   *1e-24) ;
h_paper->SetBinError(19, 1.7   *1e-24);
h_paper->SetBinError(20, 1.7   *1e-24);

  TH2D * norm_matrix = matrix_nooverflow;//matrix
  TH2D * swap_matrix = SwapIJ(matrix);


  //RooUnfoldResponse response (Ncorr,Ntruth,matrix,"test","test2");

  auto unfold = RooUnfoldBayes();
  unfold.SetRegParm( 4 );
  unfold.SetNToys( 10000 );
  unfold.SetVerbose( 3 );
  unfold.SetSmoothing( 0 );
  //RooUnfoldResponse response;
  TH1D * Ntruth2 =SmearHistogram(Ntruth,10);
  Ntruth2->Reset();
  TH2D * matrixT = (TH2D*)matrix->Clone("matrixT");
  matrixT->Rebin2D(1,2);
  ////RooUnfoldResponse * response = new RooUnfoldResponse(Ncorr, Ncorr, norm_matrix, "response", "methods");
  RooUnfoldResponse * response = new RooUnfoldResponse(Ntruth2, Ntruth2, norm_matrix);
  //RooUnfoldResponse * response2 = new RooUnfoldResponse(Ntruth2, Ntruth2, matrixT);
  //RooUnfoldResponse * response2 = new RooUnfoldResponse(nullptr, nullptr, hMatrixTranspose, "response", "methods");
  unfold.SetResponse( response );
  unfold.SetMeasured( N );

  auto unfoldSvd      = RooUnfoldSvd();
  auto unfoldTUnfold  = RooUnfoldTUnfold();
  auto unfoldIds      = RooUnfoldIds();
  auto unfoldBinByBin = RooUnfoldBinByBin();

  unfoldSvd.SetResponse( response );
  unfoldSvd.SetMeasured( N );
  unfoldTUnfold.SetResponse( response ); //response2
  unfoldTUnfold.SetMeasured( N );
  unfoldIds.SetResponse( response );
  unfoldIds.SetMeasured( N );
  unfoldBinByBin.SetResponse( response );
  unfoldBinByBin.SetMeasured( N );

  unfoldSvd.SetRegParm( 16.8 );
  unfoldIds.SetRegParm( 1 );
  //RooUnfoldSvd      unfold (&response, N, 2);   // OR
  //RooUnfoldTUnfold  unnfold (&response, Ncorr);       // OR
  //RooUnfoldIds      unfold (&response, N, 1);
  //RooUnfoldBinByBin unfold (&response, Ncorr);

  //RooUnfoldBayes   unfold (&response, Ncorr, 4);    // OR

  //correction epsilon - roughly
//2.92482883875008e19	 0.9164981512690659
//3.24589409720919e19	 0.9451778899421397
//3.59846308107099e19	 0.9364664829963452
//3.97861185996286e19	 0.9509438928203018
//4.42938570613113e19	 0.9738260452440284
//4.91216676472747e19	 0.9772884996625165
//5.44770017889487e19	 0.9816205156070249
//6.04142380673262e19	 0.9847931161835063
//6.72066064599911e19	 0.9795594459325798
//7.42855039139774e19	 0.9838919712863008
//8.2383582084585e19	 0.9879341333888024
//9.16710598198342e19	 0.9925554937661049
//1.01679126345178e20	 1.0018150250247277
//1.12657261377803e20	 0.9719442876124416
//1.25730710180846e20	 0.9635217156902284
//1.38957186545905e20	 0.9635064334138485
//1.54087938286369e20	 0.9634906417282558
//1.70866245320103e20	 0.9634748500426632
//1.90104280671985e20	 0.9634585489478579
//2.10802632740404e20	 0.9631529034202585

  //for (int i =1;i<=N->GetXaxis()->GetNbins();i++){
  //  N->Set
  //}
  cout << "==================================== UNFOLD ===================================" << endl;
  //RooUnfoldBayes   unfold (&response, Ncorr, 4);    // OR
  //RooUnfoldSvd     unfold (&response, N, 2);   // OR
  //RooUnfoldTUnfold unfold (&response, Ncorr);       // OR
  //RooUnfoldIds     unfold (&response, N, 1);
  //RooUnfoldBinByBin     unfold (&response, Ncorr);

  TH1D* hUnfold= (TH1D*) unfold.Hunfold();
  TH1D* hUnfoldSvd      = (TH1D*)unfoldSvd.Hunfold();
  //TH1D* hUnfoldTUnfold  = (TH1D*)unfoldTUnfold.Hunfold();
  TH1D* hUnfoldIds      = (TH1D*)unfoldIds.Hunfold();
  TH1D* hUnfoldBinByBin = (TH1D*)unfoldBinByBin.Hunfold();

  hUnfoldSvd->SetLineColor(kGreen+3);
  hUnfoldSvd->SetLineWidth(2);
  hUnfoldIds->SetLineColor(kAzure);
  hUnfoldIds->SetLineWidth(2);
  hUnfoldBinByBin->SetLineColor(kPink-8);
  hUnfoldBinByBin->SetLineWidth(2);
  //TCanvas* c1= new TCanvas("canvas","canvas");

  //unfold.PrintTable (cout, hTrue);
  hUnfold->SetLineColor(kRed);
  hUnfold->SetLineWidth(2);
  N->SetLineWidth(2);
  Ncorr->SetLineWidth(2);
  TCanvas *can = new TCanvas("can", "can", 0,0,1600,800);
  can->Divide(2,1);
  can->cd(1);
  gPad->SetLeftMargin(0.12);
  hUnfold->SetTitle("");//ad hoc norm based on paper 
  DivideBinWidthLog(hUnfold);
  DivideBinWidthLog(hUnfoldSvd);
  DivideBinWidthLog(hUnfoldIds);
  DivideBinWidthLog(hUnfoldBinByBin);
  DivideBinWidthLog(Ncorr);
  DivideBinWidthLog(N);

  std::cout << "Maximum is = " << hUnfold->GetMaximum() << endl;
  double NormFactor = 1.0;
  hUnfold->Scale(1./(60400.*NormFactor));
  hUnfoldSvd->Scale(1./(60400.*NormFactor));
  hUnfoldIds->Scale(1./(60400.*NormFactor));
  hUnfoldBinByBin->Scale(1./(60400.*NormFactor));
  Ncorr->Scale(1./60400.);
  N->Scale(1./60400.);

  double exponent = 3;
  for (int i = 1; i<=hUnfold->GetXaxis()->GetNbins(); i++){
    //h_frompdf->SetBinContent(i,h_frompdf->GetBinContent(i)/TMath::Power(GetBinCenterLog(h_frompdf->GetXaxis()->GetBinLowEdge(i), h_frompdf->GetXaxis()->GetBinUpEdge(i)),3));
    //h_frompdf->SetBinError(i,h_frompdf->GetBinError(i)/TMath::Power(GetBinCenterLog(h_frompdf->GetXaxis()->GetBinLowEdge(i), h_frompdf->GetXaxis()->GetBinUpEdge(i)),3));
    //h_frompdf2->SetBinContent(i,h_frompdf2->GetBinContent(i)/TMath::Power(GetBinCenterLog(h_frompdf->GetXaxis()->GetBinLowEdge(i), h_frompdf->GetXaxis()->GetBinUpEdge(i)),3));
    //h_frompdf2->SetBinError(i,h_frompdf2->GetBinError(i)/TMath::Power(GetBinCenterLog(h_frompdf->GetXaxis()->GetBinLowEdge(i), h_frompdf->GetXaxis()->GetBinUpEdge(i)),3));

    hUnfold->SetBinContent(i,hUnfold->GetBinContent(i)*TMath::Power(GetBinCenterLog(hUnfold->GetXaxis()->GetBinLowEdge(i), hUnfold->GetXaxis()->GetBinUpEdge(i)),exponent));
    Ncorr->SetBinContent(i,Ncorr->GetBinContent(i)*TMath::Power(GetBinCenterLog(Ncorr->GetXaxis()->GetBinLowEdge(i), Ncorr->GetXaxis()->GetBinUpEdge(i)),exponent));
    N->SetBinContent(i,N->GetBinContent(i)*TMath::Power(GetBinCenterLog(N->GetXaxis()->GetBinLowEdge(i), N->GetXaxis()->GetBinUpEdge(i)),exponent));
    h_paper->SetBinContent(i,h_paper->GetBinContent(i)*TMath::Power(GetBinCenterLog(h_paper->GetXaxis()->GetBinLowEdge(i), h_paper->GetXaxis()->GetBinUpEdge(i)),exponent));

    hUnfold->SetBinError(i,hUnfold->GetBinError(i)*TMath::Power(GetBinCenterLog(hUnfold->GetXaxis()->GetBinLowEdge(i), hUnfold->GetXaxis()->GetBinUpEdge(i)),exponent));
    Ncorr->SetBinError(i,Ncorr->GetBinError(i)*TMath::Power(GetBinCenterLog(Ncorr->GetXaxis()->GetBinLowEdge(i), Ncorr->GetXaxis()->GetBinUpEdge(i)),exponent));
    N->SetBinError(i,N->GetBinError(i)*TMath::Power(GetBinCenterLog(N->GetXaxis()->GetBinLowEdge(i), N->GetXaxis()->GetBinUpEdge(i)),exponent));
    h_paper->SetBinError(i,h_paper->GetBinError(i)*TMath::Power(GetBinCenterLog(h_paper->GetXaxis()->GetBinLowEdge(i), h_paper->GetXaxis()->GetBinUpEdge(i)),exponent));

    hUnfoldSvd->SetBinContent(i,hUnfoldSvd->GetBinContent(i)*TMath::Power(GetBinCenterLog(hUnfoldSvd->GetXaxis()->GetBinLowEdge(i), hUnfoldSvd->GetXaxis()->GetBinUpEdge(i)),exponent));
    hUnfoldSvd->SetBinError(i,hUnfoldSvd->GetBinError(i)*TMath::Power(GetBinCenterLog(hUnfoldSvd->GetXaxis()->GetBinLowEdge(i), hUnfoldSvd->GetXaxis()->GetBinUpEdge(i)),exponent));

    hUnfoldIds->SetBinContent(i,hUnfoldIds->GetBinContent(i)*TMath::Power(GetBinCenterLog(hUnfoldIds->GetXaxis()->GetBinLowEdge(i), hUnfoldIds->GetXaxis()->GetBinUpEdge(i)),exponent));
    hUnfoldIds->SetBinError(i,hUnfoldIds->GetBinError(i)*TMath::Power(GetBinCenterLog(hUnfoldIds->GetXaxis()->GetBinLowEdge(i), hUnfoldIds->GetXaxis()->GetBinUpEdge(i)),exponent));

    hUnfoldBinByBin->SetBinContent(i,hUnfoldBinByBin->GetBinContent(i)*TMath::Power(GetBinCenterLog(hUnfoldBinByBin->GetXaxis()->GetBinLowEdge(i), hUnfoldBinByBin->GetXaxis()->GetBinUpEdge(i)),exponent));
    hUnfoldBinByBin->SetBinError(i,hUnfoldBinByBin->GetBinError(i)*TMath::Power(GetBinCenterLog(hUnfoldBinByBin->GetXaxis()->GetBinLowEdge(i), hUnfoldBinByBin->GetXaxis()->GetBinUpEdge(i)),exponent));

      // hUnfold->SetBinContent(i,hUnfold->GetBinContent(i)*TMath::Power(TMath::Power(10,hUnfold->GetXaxis()->GetBinCenter(i)),exponent));
      // Ncorr->SetBinContent(i,Ncorr->GetBinContent(i)*TMath::Power(TMath::Power(10,Ncorr->GetXaxis()->GetBinCenter(i)),exponent));
      // N->SetBinContent(i,N->GetBinContent(i)*TMath::Power(TMath::Power(10,N->GetXaxis()->GetBinCenter(i)),exponent));
      // h_paper->SetBinContent(i,h_paper->GetBinContent(i)*TMath::Power(TMath::Power(10,h_paper->GetXaxis()->GetBinCenter(i)),exponent));
// 
      // hUnfold->SetBinError(i,hUnfold->GetBinError(i)*TMath::Power(TMath::Power(10,hUnfold->GetXaxis()->GetBinCenter(i)),exponent));
      // Ncorr->SetBinError(i,Ncorr->GetBinError(i)*TMath::Power(TMath::Power(10,Ncorr->GetXaxis()->GetBinCenter(i)),exponent));
      // N->SetBinError(i,N->GetBinError(i)*TMath::Power(TMath::Power(10,N->GetXaxis()->GetBinCenter(i)),exponent));
      // h_paper->SetBinError(i,h_paper->GetBinError(i)*TMath::Power(TMath::Power(10,h_paper->GetXaxis()->GetBinCenter(i)),exponent));
// 
      // hUnfoldSvd->SetBinContent(i,hUnfoldSvd->GetBinContent(i)*TMath::Power(TMath::Power(10,hUnfoldSvd->GetXaxis()->GetBinCenter(i)),exponent));
      // hUnfoldSvd->SetBinError(i,hUnfoldSvd->GetBinError(i)*TMath::Power(TMath::Power(10,hUnfoldSvd->GetXaxis()->GetBinCenter(i)),exponent));
// 
      // hUnfoldIds->SetBinContent(i,hUnfoldIds->GetBinContent(i)*TMath::Power(TMath::Power(10,hUnfoldIds->GetXaxis()->GetBinCenter(i)),exponent));
      // hUnfoldIds->SetBinError(i,hUnfoldIds->GetBinError(i)*TMath::Power(TMath::Power(10,hUnfoldIds->GetXaxis()->GetBinCenter(i)),exponent));
// 
      // hUnfoldBinByBin->SetBinContent(i,hUnfoldBinByBin->GetBinContent(i)*TMath::Power(TMath::Power(10,hUnfoldBinByBin->GetXaxis()->GetBinCenter(i)),exponent));
      // hUnfoldBinByBin->SetBinError(i,hUnfoldBinByBin->GetBinError(i)*TMath::Power(TMath::Power(10,hUnfoldBinByBin->GetXaxis()->GetBinCenter(i)),exponent));

  }
  hUnfold->SetMinimum(2e36);
  hUnfold->SetMaximum(1.5e38);
  hUnfoldSvd->SetMinimum(2e36);
  hUnfoldIds->SetMinimum(2e36);
  hUnfoldBinByBin->SetMinimum(2e36);
  hUnfold->GetYaxis()->SetTitle("J(E)#times E^{3} [km^{-2} yr^{-1} sr^{-1} eV^{2}]");
  hUnfold->GetXaxis()->SetTitleOffset(1.15);
  hUnfold->GetXaxis()->SetTitle("log_{10}(E) [eV]");
  //Ncorr->Divide(hUnfold);
  
  N->SetLineColor(kBlack);
  
  Ncorr->SetLineColor(kBlue);
  Ncorr->SetLineStyle(2);

  for (int i=1;i<=Ncorr->GetXaxis()->GetNbins();i++){
    cout << i << " = " << N->GetBinContent(i) << endl;
  }
  
  gStyle->SetOptStat(0);
  gPad->SetLogy(1);
  TLegend * leg = new TLegend(0.15,0.15,0.5,0.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03); 
  
  //hTrue->SetLineColor(8);
  //hTrue->Draw("SAME");
  //Ncorr->Divide(N);
  //gPad->SetLogx(1);
  //Ncorr->SetMaximum(1.);
  //Ncorr->SetMinimum(0.88);
  //Ncorr->Draw("hist");

  
  //h_paper->Scale(1./1.12);
  h_paper->SetLineColor(kGreen+3);
  h_paper->SetLineWidth(2);
  h_paper->SetLineStyle(2);
  

  hUnfold->Draw("hist e1");
  hUnfoldSvd->Draw("hist e1 same");
  hUnfoldIds->Draw("hist e1 same");
  hUnfoldBinByBin->Draw("hist e1 same");
  //N->Draw("hist SAME");
  Ncorr->Draw("hist e1 same");
  h_paper->Draw("hist same e1");
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();
    TGraphAsymmErrors* graph2 = new TGraphAsymmErrors();
      for (int i = 1; i <= h_frompdf->GetNbinsX(); ++i) {
        double x = h_frompdf->GetBinCenter(i);      // Bin center
        double y = h_frompdf->GetBinContent(i);     // Bin content
        double ex = 0;//h_frompdf->GetBinWidth(i) / 2;  // Symmetric x-errors (half the bin width)
        double ey = h_frompdf->GetBinError(i);      // Symmetric y-errors

        graph->SetPoint(i - 1, x, y);               // Set point
        graph->SetPointError(i - 1, ex, ex, ey, ey); // Set symmetric errors

        double x2 = h_frompdf2->GetBinCenter(i);      // Bin center
        double y2 = h_frompdf2->GetBinContent(i);     // Bin content
        double ex2 = 0;//h_frompdf->GetBinWidth(i) / 2;  // Symmetric x-errors (half the bin width)
        double ey2 = h_frompdf2->GetBinError(i);      // Symmetric y-errors

        graph2->SetPoint(i - 1, x2, y2);               // Set point
        graph2->SetPointError(i - 1, ex2, ex2, ey2, ey2); // Set symmetric err
    }

    // Draw the graph
    //TCanvas* c = new TCanvas("c", "Graph from Histogram", 800, 600);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(2);
    graph->SetLineWidth(2);

    graph->Draw("P same"); // "AP" for points with error bars

    graph2->SetMarkerColor(kRed);
    graph2->SetLineColor(kRed);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerSize(2);
    graph2->SetLineWidth(2);

    graph2->Draw("P same"); // "AP" for points with error bars

  //h_frompdf->Draw("e1");
  
  //leg->AddEntry(N, "N");
  
  leg->AddEntry((TObject*)0, "https://arxiv.org/pdf/2008.06486", "");
  leg->AddEntry(Ncorr, "J(E)#times E^{3} derived from N Table VI");
  leg->AddEntry(h_paper, "J(E)#times E^{3} Table VI (last column)");
  leg->AddEntry(graph, "Pierre Auger Fig.14", "p");
  leg->AddEntry(graph2, "Telescope Array Fig.14", "p");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(hUnfold, "J(E)#times E^{3}(N) unfolded Bayes");
  leg->AddEntry(hUnfoldSvd, "J(E)#times E^{3}(N) unfolded Svd (par. = 16.8)");
  leg->AddEntry(hUnfoldIds, "J(E)#times E^{3}(N) unfolded Ids (par. = 1)");
  leg->AddEntry(hUnfoldBinByBin, "J(E)#times E^{3}(N) unfolded BinByBin");
  //leg->AddEntry((TObject*)0, "Bayes Unfolding", "");
  leg->Draw("same");

  TCanvas * can2 = new TCanvas("can2", "can2", 0,0,800,800);
  can2->cd();
  //TH1D * unfoldBinByBin2 = (TH1D*)hUnfoldBinByBin->Clone("unfoldBinByBin2");
  //unfoldBinByBin2->SetMaximum(2.);
  //unfoldBinByBin2->SetMinimum(0.);
  //unfoldBinByBin2->Divide(h_frompdf);
  //unfoldBinByBin2->Draw("hist same");
  norm_matrix->Draw("colz text");
  //h_paper->Divide(Ncorr);
  //h_paper->Draw("hist e1");

  //TCanvas * can3 = new TCanvas("can3", "can3", 0,0,800,800);
  can->cd(2);
  TH1D * unity = (TH1D*)h_frompdf2->Clone("unity");
  TH1D * unity_or = (TH1D*)h_frompdf->Clone("unity_or");
  TH1D * unfoldBayes2 = (TH1D*)hUnfold->Clone("unfoldBayes");
  TH1D * unfoldBinByBin2 = (TH1D*)hUnfoldBinByBin->Clone("unfoldBinByBin2");
  TH1D * unfoldSvd2 = (TH1D*)hUnfoldSvd->Clone("unfoldSvd2");
  TH1D * unfoldIds2 = (TH1D*)hUnfoldIds->Clone("unfoldIds2");
  TH1D * h_frompdf2b = (TH1D*)h_frompdf2->Clone("h_frompdf2b");
  
  unity->SetMaximum(2.);
  unity->SetMinimum(0.);
  unity->Divide(h_frompdf);
  unity->SetMarkerSize(0);
  unity->SetFillColorAlpha(kRed-10, 1.);
  unity->SetFillStyle(3008);

  unity_or->SetMaximum(2.);
  unity_or->SetMinimum(0.0001);
  unity_or->Divide(h_frompdf);
  unity_or->SetMarkerSize(0);
  unity_or->SetFillColor(kBlue -10);

  h_frompdf2b->Divide(h_frompdf);
  unfoldBinByBin2->Divide(h_frompdf);
  unfoldBayes2->Divide(h_frompdf);
  unfoldSvd2->Divide(h_frompdf);
  unfoldIds2->Divide(h_frompdf);

  TH1D * Ncorr2 = (TH1D*)Ncorr  ->Clone("Ncorr2");
  TH1D * h_paper2 = (TH1D*)h_paper->Clone("h_paper2");

  Ncorr2  ->Divide(h_frompdf);
  h_paper2->Divide(h_frompdf);

  unity_or->SetTitle("");
  unity_or->GetXaxis()->SetTitle(hUnfold->GetXaxis()->GetTitle());
  unity_or->GetYaxis()->SetTitle("Ratio = Unfolded / Pierre Auger Fig. 14");
  unity_or->SetBinError(19,2.);
  unity_or->SetBinError(20,2.);
  unity_or->GetXaxis()->SetTitleOffset(1.15);
  unity_or->Draw("e2");
  TH1D * unity2_or = (TH1D*)unity_or->Clone("unity2_or");
  unity2_or->SetFillColor(0);
  unity2_or->SetLineStyle(2);
  unity2_or->SetLineColor(kBlue-10);
  unity2_or->SetMarkerSize(2);
  unity2_or->Draw("hist same");
  unity2_or->Draw("p same"); 

  unity->SetTitle("");
  unity->GetXaxis()->SetTitle(hUnfold->GetXaxis()->GetTitle());
  unity->GetYaxis()->SetTitle("Ratio = Unfolded / Pierre Auger Fig. 14");
  unity->SetBinError(19,2.);
  unity->SetBinError(20,2.);
  unity->Draw("e2 same");
  TH1D * unity2 = (TH1D*)unity->Clone("unity2");
  unity2->SetFillColor(0);
  unity2->SetLineStyle(2);
  unity2->SetLineColor(kRed-10);
  unity2->SetMarkerSize(2);
  //unity2->Draw("hist same");
  unity2->Draw("p same"); 


  

  TGraphAsymmErrors *gunfoldBayes2 = ShiftTickUsingTGraph(unfoldBayes2, -0.02);
  TGraphAsymmErrors *gunfoldBinByBin2 = ShiftTickUsingTGraph(unfoldBinByBin2, 0.01);
  TGraphAsymmErrors *gunfoldSvd2 = ShiftTickUsingTGraph(unfoldSvd2, -0.01);
  TGraphAsymmErrors *gunfoldIds2 = ShiftTickUsingTGraph(unfoldIds2, 0.);
  gunfoldBayes2->Draw("p same");
  gunfoldBinByBin2->Draw("p same");
  gunfoldSvd2->Draw("p same");
  gunfoldIds2->Draw("p same");
  //unfoldBayes2->Draw("hist e1 same");
  //unfoldBinByBin2->Draw("hist e1 same");
  //unfoldSvd2->Draw("hist e1 same");
  //unfoldIds2->Draw("hist e1 same");

  Ncorr2  ->Draw("hist e1 same");
  h_paper2->Draw("hist e1 same");

  h_frompdf2b->SetMarkerColor(kRed);
  h_frompdf2b->SetLineColor(kRed);
  h_frompdf2b->SetMarkerStyle(21);
  h_frompdf2b->SetMarkerSize(2);
  h_frompdf2b->SetLineWidth(2);

  TLegend * leg2 = new TLegend(0.174969, 0.130085, 0.535244, 0.437134);


  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.03); 
  leg2->AddEntry((TObject*)0, "https://arxiv.org/pdf/2008.06486", "");
  leg2->AddEntry(Ncorr  , "J(E)#times E^{3} derived from N Table VI");
  leg2->AddEntry(h_paper, "J(E)#times E^{3} Table VI (last column)");
  leg2->AddEntry(graph, "Pierre Auger Fig.14", "p");
  leg2->AddEntry(graph2, "Telescope Array Fig.14", "p");
  leg2->AddEntry((TObject*)0, "", "");
  leg2->AddEntry(unfoldBayes2, "J(E)#times E^{3}(N) unfolded Bayes");
  leg2->AddEntry(unfoldSvd2, "J(E)#times E^{3}(N) unfolded Svd (par. = 16.8)");
  leg2->AddEntry(unfoldIds2, "J(E)#times E^{3}(N) unfolded Ids (par. = 1)");
  leg2->AddEntry(unfoldBinByBin2, "J(E)#times E^{3}(N) unfolded BinByBin");
  //leg->AddEntry((TObject*)0, "Bayes Unfolding", "");
  leg2->Draw("same");

  vector<TArrow*> arrows;
    TArrow * temp;
    for (int i=1; i <= unity->GetXaxis()->GetNbins(); i++){
        if (unity->GetBinContent(i) > 2.0){
          cout << i << endl;
          temp = new TArrow(unity->GetXaxis()->GetBinCenter(i), 1.8, unity->GetXaxis()->GetBinCenter(i), 2.0 , 0.01 , "|>");  
          temp->SetAngle(40);
          temp->SetLineWidth(2);
          temp->SetFillColor(kRed);
          temp->SetLineColor(kRed);
          arrows.push_back(temp);
        }
        if (unity->GetBinContent(i) < 0.0){
          temp = new TArrow(unity->GetXaxis()->GetBinCenter(i), 0.65, unity->GetXaxis()->GetBinCenter(i), 0.0 , 0.01 , "|>");  
          temp->SetAngle(40);
          temp->SetLineWidth(2);
          temp->SetFillColor(kRed);
          temp->SetLineColor(kRed);
          arrows.push_back(temp);
        }
        
    }

    for (int i=0; i < arrows.size(); i++){
      arrows[i]->Draw();
      cout << arrows[i] << endl;
    }
    gPad->RedrawAxis();
    can2->Update();
  //h_frompdf2b->Draw("hist e1 same");
  //norm_matrix->Draw("colz text");
  //h_paper->Divide(Ncorr);
  //h_paper->Draw("hist e1");

  can->SaveAs("AstroUnfold.pdf");
  can->SaveAs("AstroUnfold.png");

  return 0;
}


