///////////////////////////////////////////////////////////////////////////////
// prerequisites:
// --------------
// .L muon_capture/scripts/muon_capture_yields.C
//
// mc = new muon_capture();
// mc->fit_proton_energy_spectrum()
//
///////////////////////////////////////////////////////////////////////////////
#include "TCanvas.h"
#include "TLegend.h"

#include "ana/muon_capture.hh"

#include "murat/plot/murat_plot_functions.hh"

muon_capture* mc(nullptr);
//-----------------------------------------------------------------------------
muon_capture::muon_capture() {

  alcap = nullptr;
  twist = nullptr;

  init_data();
  
  //  init_proton_fit();

  f_prot = nullptr;
  init_proton_fit();

  f_prot_2 = nullptr;
  init_proton_fit_2();

  f_deut = nullptr;
  init_deuteron_fit();

  init_hgf_spectrum();
}

//-----------------------------------------------------------------------------
// parameterization by Ed Hungerford : https://arxiv.org/pdf/1803.08403.pdf
// Par[0] = 0: energy
//        = 1: momentum
// Par[1]    : normalization of the integral
//-----------------------------------------------------------------------------
double muon_capture::hgf_spectrum(double* X, double* Par) {

  //taken from GMC
  //
  //   Ed Hungerford  Houston University May 17 1999
  //   Rashid Djilkibaev New York University (modified) May 18 1999
  //
  //   e - proton kinetic energy (MeV)
  //   p - proton Momentum (MeV/c)
  //
  //   Generates a proton spectrum similar to that observed in
  //   mu capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

  //these numbers are in MeV!!!!
  static const double emn  = 1.4; // replacing par1 from GMC
  static const double par2 = 1.3279;
  static const double par3 = 17844.0;
  static const double par4 = .32218;
  static const double par5 = 100.;
  static const double par6 = 10.014;
  static const double par7 = 1050.;
  static const double par8 = 5.103;

  static double MP = 938. ; 
    
  double spectrumWeight;

  double e, p, w1(1.);

  if (Par[0] == 0) e = X[0];
  else {
    p  = X[0];
    e  = p*p/(2*MP);
    w1 = p/MP;
  }

  if (e >= 20) {
    spectrumWeight=par5*exp(-(e-20.)/par6);
  }
  else if (e >= 8.0 && e <= 20.0) {
    spectrumWeight=par7*exp(-(e-8.)/par8);
  }
  else if (e > emn) {
    double xw=(1.-emn/e);
    double xu=std::pow(xw,par2);
    double xv=par3*exp(-par4*e);
    spectrumWeight=xv*xu;
  }
  else {
    spectrumWeight = 0.;
  }

  return w1*Par[1]*spectrumWeight;
}

//-----------------------------------------------------------------------------
double muon_capture::fitf_prot(double* X, double* P) {
  double f;

  double e  = X[0];
  double e1 = P[3];
  //  double e2 = 10;

  if (e <= 0) {
    f = 0;
  }
  else if (e < e1) {
    f = pow(e/(1+P[1]*e),P[2])*exp(-e/P[4]);
  }
  else {
    double c4 = pow(e1/(1+P[1]*e1),P[2])*exp(-e1/P[4]);
    f = c4*exp(-(e-e1)/P[5]);
  }

  return P[0]*f;
}

//-----------------------------------------------------------------------------
void muon_capture::init_proton_fit() {

  if (f_prot == nullptr) {
    f_prot = new TF1("f_prot",fitf_prot,0,50,6);
    f_prot->SetLineColor(kBlue+2);
  }

  f_prot->SetParameter(0,0.01);
  f_prot->SetParameter(1,0.5);
  f_prot->SetParameter(2,2.);
  f_prot->SetParameter(3,7.);
  f_prot->SetParameter(4,3.);
  f_prot->SetParameter(5,6.);
}

//-----------------------------------------------------------------------------
double muon_capture::fitf_prot_2(double* X, double* P) {
  double f;

  double e  = X[0];
  double e1 = P[3];
  //  double e2 = 10;

  if (e <= 0) {
    f = 0;
  }
  else if (e < e1) {
    f = pow(e,P[6])/(1+P[1]*P[1]*(e-P[2])*(e-P[2]))*exp(-e/P[4]);
  }
  else {
    double c4 = pow(e1,P[6])/(1+P[1]*P[1]*(e1-P[2])*(e1-P[2]))*exp(-e1/P[4]);
    f = c4*exp(-(e-e1)/P[5]);
  }

  return P[0]*f;
}

//-----------------------------------------------------------------------------
void muon_capture::init_proton_fit_2() {

  if (f_prot_2 == nullptr) {
    f_prot_2 = new TF1("f_prot",fitf_prot_2,0,50,7);
    f_prot_2->SetLineColor(kBlue+2);
  }

  f_prot_2->SetParameter(0,0.005);
  f_prot_2->SetParameter(1,0.8);
  f_prot_2->SetParameter(2,3.5); 
  f_prot_2->SetParameter(3,7.0);
  f_prot_2->SetParameter(4,3.4);
  f_prot_2->SetParameter(5,5.7);
  f_prot_2->SetParameter(6,0.8);
}

//-----------------------------------------------------------------------------
// fit parameters in the deuteron energy spectrum fit based on the proton fit
void muon_capture::init_deuteron_fit() {

  if (f_deut == nullptr) {
    f_deut = new TF1("f_deut",fitf_prot ,0,50,6);
    f_deut->SetLineColor(kBlue+2);
  }

  f_deut->SetParameter(0,0.01);
  f_deut->FixParameter(1,0.5);
  f_deut->SetParameter(2,2);
  f_deut->FixParameter(3,7.755);
  f_deut->SetParameter(4,3);
  f_deut->SetParameter(5,6.);
}

//-----------------------------------------------------------------------------
void muon_capture::init_hgf_spectrum() {
					// energy spectrum
  fe_hgf = new TF1("fe_hgf",muon_capture::hgf_spectrum,0,100,2);
  fe_hgf->SetParameter(0,0);
  fe_hgf->SetParameter(1,1);
  fe_hgf->SetParameter(1,0.05/fe_hgf->Integral(0,100));
  fe_hgf->SetNpx(1000);

  fe_hgf->SetTitle("Ejected proton energy");
  fe_hgf->GetXaxis()->SetTitle("E, MeV");

					// momentum 
  fp_hgf = new TF1("fp_hgf",muon_capture::hgf_spectrum,0,1000,2);
  fp_hgf->SetParameter(0,1);
  fp_hgf->SetParameter(1,1);
  fp_hgf->SetParameter(1,0.05/fp_hgf->Integral(0,1000));
  fp_hgf->SetNpx(1000);
  
  fp_hgf->SetTitle("Ejected proton momentum");
  fp_hgf->GetXaxis()->SetTitle("P, MeV/c");
}

//-----------------------------------------------------------------------------
void muon_capture::init_data() {
  // if (! gInterpreter->IsLoaded("muon_capture/scripts/alcap_data.C")) {
  //   gInterpreter->LoadMacro("muon_capture/scripts/alcap_data.C");
  // }

  // if (! gInterpreter->IsLoaded("muon_capture/scripts/twist_data.C")) {
  //   gInterpreter->LoadMacro("muon_capture/scripts/twist_data.C");
  // }

  // gInterpreter->ProcessLine("twist = new twist_data();");
  // gInterpreter->ProcessLine("alcap = new alcap_data();");


  twist = new twist_data();
  alcap = new alcap_data();

  printf("muon_capture initialized\n");
}

//-----------------------------------------------------------------------------
// FCN=21.3977 FROM MIGRAD    STATUS=CONVERGED     413 CALLS         414 TOTAL
//                     EDM=2.60755e-07    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.0 per cent
//  EXT PARAMETER                                   STEP         FIRST   
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  p0           9.92698e-03   9.60414e-04   1.73751e-07  -2.25181e+00
//   2  p1           5.04814e-01   8.76443e-03  -2.45901e-06   8.56883e-02
//   3  p2           1.76282e+00   6.47681e-01  -4.49109e-05  -6.58915e-03
//   4  p3           7.75550e+00   4.04286e-03  -5.85623e-07   2.15198e-02
//   5  p4           3.44504e+00   3.42612e-01  -3.64853e-06  -8.02442e-03
//   6  p5           5.86909e+00   1.36226e-01  -4.86316e-05   5.45286e-03
//-----------------------------------------------------------------------------
void muon_capture::fit_proton_energy_spectrum() {

					// at this point all interpolations exist,
					// create points for the fit 
					// start from TWIST proton energy (twist.pe)
  struct dat_t {
    double x;
    double y; 
    double ey;
  };

  dat_t* ad[1000];			// have less than 100 points
  int nptot = 0;

  double x[1000], y[1000], ex[1000], ey[1000];
  for (int i=0; i<1000; i++) ex[i] = 0;
//-----------------------------------------------------------------------------
// start with TWIST
//-----------------------------------------------------------------------------
  int np = twist->pe.s->fN;

  int    npp  = 20;
  double step = (twist->pe.s->fX[np-1]-twist->pe.s->fX[0])/(npp-1);

  for (int i=0; i<npp; i++) {
    ad[nptot]     = new dat_t();

    double x      = twist->pe.s->fX[0]+step*(i-1);

    ad[nptot]->x  = x;
    ad[nptot]->y  = twist->pe.s->GetFunc()->Eval(x);

    double eu     = twist->pe.su->GetFunc()->Eval(x);
    double el     = twist->pe.sl->GetFunc()->Eval(x);

    ad[nptot]->ey = (eu-el)/2;

    nptot        = nptot+1;
  }
//-----------------------------------------------------------------------------
// add AlCap
//-----------------------------------------------------------------------------
  int npa = alcap->pe.s->fN;

  step = (alcap->pe.s->fX[npa-1]-alcap->pe.s->fX[0])/(npp-1);

  for (int i=0; i<npp; i++) {
    ad[nptot]     = new dat_t();

    double x      = alcap->pe.s->fX[0]+step*(i-1);

    ad[nptot]->x  = x;
    ad[nptot]->y  = alcap->pe.s->GetFunc()->Eval(x);

    double eu     = alcap->pe.su->GetFunc()->Eval(x);
    double el     = alcap->pe.sl->GetFunc()->Eval(x);

    ad[nptot]->ey = (eu-el)/2;

    nptot        = nptot+1;
  }

  //  printf(" nptot = %3i\n",nptot);
//-----------------------------------------------------------------------------
// sort data according to ascending X
//-----------------------------------------------------------------------------
  for (int i1=0; i1<nptot-1; i1++) {
    dat_t* ad1 = ad[i1];
    for (int i2=i1+1; i2<nptot; i2++) {
      dat_t* ad2 = ad[i2];

      if (ad1->x > ad2->x) {
	ad[i1] = ad2;
	ad[i2] = ad1;
	ad1    = ad2;
      }
    }
  }

//-----------------------------------------------------------------------------
// build a single TGraphErrors
//-----------------------------------------------------------------------------
  for (int i=0; i<nptot; i++) {
    x [i] = ad[i]->x;
    y [i] = ad[i]->y;
    ey[i] = ad[i]->ey;

    //    printf("i,x,y,ey : %3i %12.5e %12.5e %12.5e\n",i,x[i],y[i],ey[i]);
  }

  //  printf("npp = %3i\n",npp);

  gr_fit_pe = new TGraphErrors(nptot,x,y,ex,ey);
  gr_fit_pe->SetName("gr_fit_pe");
  gr_fit_pe->SetTitle("proton energy");

  //  TCanvas* c_fit = new TCanvas("c_fit","c fit",1300,900);

  //  c_fit->cd();

  gr_fit_pe->SetMarkerStyle(20);
  gr_fit_pe->SetMarkerSize(1);
  gr_fit_pe->SetMarkerColor(kRed+2);
  gr_fit_pe->SetLineColor  (kRed+2);
  gr_fit_pe->SetFillStyle(3003);
  gr_fit_pe->SetFillColor(kRed+2);

  gr_fit_pe->Draw("ap");

  gr_fit_pe->Fit(f_prot,"","");
  //  gr_fit_pe->Fit(f_prot_2,"","");
}

//-----------------------------------------------------------------------------
// AlCap doesn't have deuteron data
// FCN=2.30803 FROM HESSE     STATUS=NOT POSDEF     23 CALLS        1176 TOTAL
//                     EDM=3.4442e-07    STRATEGY= 1      ERR MATRIX NOT POS-DEF
//  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  p0           1.10022e-03   3.35040e-04   4.97657e-09  -2.56153e+01
//   2  p1           5.00000e-01     fixed    
//   3  p2           6.54336e+00   7.82961e-01   9.88501e-06  -1.31436e-02
//   4  p3           7.75500e+00     fixed    
//   5  p4           2.49239e+00   2.59680e-01   3.71896e-06  -3.47638e-02
//   6  p5           7.59651e+00   2.49517e-01   3.70865e-05  -4.15849e-04
//-----------------------------------------------------------------------------
void muon_capture::fit_deuteron_energy_spectrum() {

  dat_t* ad[1000];			// have less than 1000 points

  double x[1000], y[1000], ex[1000], ey[1000];
  for (int i=0; i<1000; i++) ex[i] = 0;
//-----------------------------------------------------------------------------
// start with TWIST
//-----------------------------------------------------------------------------
  int    npp  = 40;
  int    np   = twist->de.s->fN;

  double step = (twist->de.s->fX[np-1]-twist->de.s->fX[0])/(npp-1);

  printf(" deut: n, x[0], x[n-1] : %3i %12.5e %12.5e \n",twist->de.s->fN,twist->de.s->fX[0],twist->de.s->fX[np-1]);

  int nptot = 0;
  for (int i=0; i<npp; i++) {
    ad[nptot]     = new dat_t();

    double x      = twist->de.s->fX[0]+step*i;

    ad[nptot]->x  = x;
    ad[nptot]->y  = twist->de.s->GetFunc()->Eval(x);

    double eu     = twist->de.su->GetFunc()->Eval(x);
    double el     = twist->de.sl->GetFunc()->Eval(x);

    ad[nptot]->ey = (eu-el)/2;

    nptot         = nptot+1;
  }

//-----------------------------------------------------------------------------
// build a single TGraphErrors
//-----------------------------------------------------------------------------
  for (int i=0; i<nptot; i++) {
    x [i] = ad[i]->x;
    y [i] = ad[i]->y;
    ey[i] = ad[i]->ey;

    printf("i,x,y,ey : %3i %12.5e %12.5e %12.5e\n",i,x[i],y[i],ey[i]);
  }

  printf("npp = %3i\n",npp);

  gr_fit_de = new TGraphErrors(nptot,x,y,ex,ey);
  gr_fit_de->SetTitle("");
  gr_fit_de->GetXaxis()->SetTitle("energy, MeV");
  gr_fit_de->GetXaxis()->SetLimits(0,50);

  gr_fit_de->SetMarkerStyle(20);
  gr_fit_de->SetMarkerSize(1);
  gr_fit_de->SetMarkerColor(kMagenta+2);
  gr_fit_de->SetLineColor  (kMagenta+2);
  gr_fit_de->SetFillStyle(3001);
  gr_fit_de->SetFillColor(kMagenta+2);

  gr_fit_de->Draw("ap");
//-----------------------------------------------------------------------------
// need to redefine parameters
//-----------------------------------------------------------------------------
  gr_fit_de->Fit(f_deut,"","");
}

//-----------------------------------------------------------------------------
void muon_capture::plot_energy_spectra() {

  //  double protons_per_capture(0.05);

  TCanvas* c_e = new TCanvas("c_e","c_e",1300,800);
  
  //  plot_ejected_proton_spectrum("e",protons_per_capture);

  twist->pe.gr->Draw("same,l3");
  twist->de.gr->Draw("same,l3");
  alcap->pe.gr->Draw("same,l3");

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry("ep_fun"   ,"Mu2e protons"   ,"f");
  leg->AddEntry(twist->pe.gr,"TWIST protons"  ,"f");
  leg->AddEntry(twist->de.gr,"TWIST deuterons","f");
  leg->AddEntry(alcap->pe.gr,"AlCap protons"  ,"f");

  leg->Draw();

  c_e->Modified();
  c_e->Update();
//-----------------------------------------------------------------------------
// plot momentum distributions
//-----------------------------------------------------------------------------
  TCanvas* c_p = new TCanvas("c_p","c_p",1300,800);
  
  //  plot_ejected_proton_spectrum("p",protons_per_capture);

  twist->pp.gr->Draw("same,l3");
  twist->dp.gr->Draw("same,l3");
  alcap->pp.gr->Draw("same,l3");

  TLegend* leg2 = new TLegend(0.6,0.6,0.9,0.9);
  leg2->AddEntry("ep_fun"   ,"Mu2e protons"   ,"f");
  leg2->AddEntry(twist->pp.gr,"TWIST protons"  ,"f");
  leg2->AddEntry(twist->dp.gr,"TWIST deuterons","f");
  leg2->AddEntry(alcap->pp.gr,"AlCap protons"  ,"f");

  leg2->Draw();

  c_p->Modified();
  c_p->Update();
}


//-----------------------------------------------------------------------------
void muon_capture::plot(int Figure, int Print) {

  TString canvas_name = Form("Figure_%04i",Figure);
  TCanvas* c = new TCanvas(canvas_name,canvas_name,1300,900);

  if (Figure == 1) {// fig_0001: proton momentum, linear
    // plot name: fig_001_proton_energy_log.eps
    c->SetLogy(0);
    twist->pp.gr->GetXaxis()->SetTitle("proton momentum, MeV/c");
    twist->pp.gr->GetXaxis()->SetLimits(0,350);
    twist->pp.gr->Draw("alp");
    twist->pp.gr->Draw("p3,same");
    alcap->pp.gr->Draw("p3,same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pp.gr,"TWIST","f");
    leg->AddEntry(alcap->pp.gr,"AlCap","f");
    leg->Draw();
  }

  if (Figure == 2) { // fig_0002: proton momentum, log
    // plot name: fig_002_proton_momentum_log.eps
    c->SetLogy(1);
    twist->pp.gr->GetXaxis()->SetTitle("proton momentum, MeV/c");
    twist->pp.gr->GetXaxis()->SetLimits(0,350);
    twist->pp.gr->Draw("alp");
    twist->pp.gr->Draw("p3,same");
    alcap->pp.gr->Draw("p3,same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pp.gr,"TWIST","f");
    leg->AddEntry(alcap->pp.gr,"AlCap","f");
    leg->Draw();
  }

  if (Figure == 3) { // fig_0003: proton momenum, overlay with Mu2e model
    // plot name: fig_003_proton_momentum_twist_alcap_mu2e.eps
    //    c->SetLogy(1);

    twist->pp.gr->GetXaxis()->SetTitle("proton momentum, MeV");
    twist->pp.gr->GetXaxis()->SetLimits(0,350);
    twist->pp.gr->GetYaxis()->SetRangeUser(0,0.0009);
    twist->pp.gr->Draw("alp3");
    //    twist->pe.gr->Draw("p3,same");
    
    fp_hgf->Draw("same");

    alcap->pp.gr->Draw("p3,same");

    TGraph* gr = (TGraph*) twist->ppg.gr->Clone("gr");
    int n = gr->GetN();
    double sf = 0.05/gr->Integral();

    for (int i=0; i<n; i++) {
      gr->GetY()[i] = gr->GetY()[i]*sf;
    }
    
    gr->Draw("same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(fe_hgf,"Mu2e","f");
    leg->AddEntry(twist->pp.gr,"TWIST","f");
    leg->AddEntry(alcap->pp.gr,"AlCap","f");
    leg->AddEntry(gr,"G4 precompound, I=0.05","l");
    leg->Draw();
  }

  if (Figure == 11) {// fig_0012: proton energy, linear
    // plot name: fig_011_proton_energy.eps
    c->SetLogy(0);
    twist->pe.gr->GetXaxis()->SetTitle("proton energy, MeV");
    twist->pe.gr->GetXaxis()->SetLimits(0,50);
    twist->pe.gr->Draw("alp");
    twist->pe.gr->Draw("p3,same");
    alcap->pe.gr->Draw("p3,same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pe.gr,"TWIST","f");
    leg->AddEntry(alcap->pe.gr,"AlCap","f");
    leg->Draw();
  }

  if (Figure == 12) { // fig_0012: proton energy, log 
    // plot name: fig_012_proton_energy_log.eps
    c->SetLogy(1);
    twist->pe.gr->GetXaxis()->SetTitle("proton energy, MeV");
    twist->pe.gr->GetXaxis()->SetLimits(0,50);
    twist->pe.gr->Draw("alp");
    twist->pe.gr->Draw("p3,same");
    alcap->pe.gr->Draw("p3,same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pe.gr,"TWIST","f");
    leg->AddEntry(alcap->pe.gr,"AlCap","f");
    leg->Draw();
  }

  if (Figure == 13) { // fig_0013: proton energy, overlay with Mu2e model
    // name: fig_013_proton_energy_twist_alcap_mu2e.eps
    //    c->SetLogy(1);

    twist->pe.gr->GetXaxis()->SetTitle("proton energy, MeV");
    twist->pe.gr->GetXaxis()->SetLimits(0,50);
    twist->pe.gr->GetYaxis()->SetLimits(0,0.001);
    twist->pe.gr->Draw("alp3");
    //    twist->pe.gr->Draw("p3,same");

    fe_hgf->Draw("same");

    alcap->pe.gr->Draw("p3,same");

    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(fe_hgf,"Mu2e","f");
    leg->AddEntry(twist->pe.gr,"TWIST","f");
    leg->AddEntry(alcap->pe.gr,"AlCap","f");
    leg->Draw();
  }

  if (Figure == 14) { // fig_0014: twist proton and deutron energy spectra linear 
    // name: fig_014_twist_proton_deutron_energy.eps
    //    c->SetLogy(1);

    twist->pe.gr->SetTitle("");
    twist->pe.gr->GetXaxis()->SetTitle("energy, MeV");
    twist->pe.gr->GetXaxis()->SetLimits(0,50);
    twist->pe.gr->GetYaxis()->SetLimits(0,0.001);
    twist->pe.gr->Draw("alp3");

    alcap->pe.gr->Draw("lp3,same");

    twist->de.gr->GetXaxis()->SetTitle("energy, MeV");
    twist->de.gr->GetXaxis()->SetLimits(0,50);
    twist->de.gr->GetYaxis()->SetLimits(0,0.001);
    twist->de.gr->Draw("lp3,same");


    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pe.gr,"TWIST protons","f");
    leg->AddEntry(alcap->pe.gr,"AlCap","f");
    leg->AddEntry(twist->de.gr,"TWIST deuterons","f");
    leg->Draw();
  }

  if (Figure == 15) { // fig_0015: twist proton and deutron energy spectra log
    // name: fig_015_twist_proton_deutron_energy_log.eps

    c->SetLogy(1);

    twist->pe.gr->SetTitle("");
    twist->pe.gr->GetXaxis()->SetTitle("energy, MeV");
    twist->pe.gr->GetXaxis()->SetLimits(0,50);
    twist->pe.gr->GetYaxis()->SetLimits(0,0.001);
    twist->pe.gr->Draw("alp3");

    alcap->pe.gr->Draw("lp3,same");

    twist->de.gr->GetXaxis()->SetTitle("energy, MeV");
    twist->de.gr->GetXaxis()->SetLimits(0,50);
    twist->de.gr->GetYaxis()->SetLimits(0,0.001);
    twist->de.gr->Draw("lp3,same");


    TLegend* leg = new TLegend(0.6,0.7,0.8,0.8);
    leg->SetLineWidth(0);
    leg->AddEntry(twist->pe.gr,"TWIST protons","f");
    leg->AddEntry(alcap->pe.gr,"AlCap","f");
    leg->AddEntry(twist->de.gr,"TWIST deuterons","f");
    leg->Draw();
  }

  if (Figure == 21) { // fit proton energy spectrum
    // plot name: fig_021_fit_proton_energy_spectrum.eps
    fit_proton_energy_spectrum();
  }

  if (Figure == 22) { // fit deuteron energy spectrum
    // plot name: fig_022_fit_deuteron_energy_spectrum.eps
    fit_deuteron_energy_spectrum();
  }




}
