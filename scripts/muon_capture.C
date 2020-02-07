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
#include "TGraphErrors.h"
#include "TInterpreter.h"

#include "twist_data.hh"
#include "alcap_data.hh"
#include "murat/plot/murat_plot_functions.hh"


class muon_capture {
public:

  struct dat_t {
    double x;
    double y; 
    double ey;
  };

  twist_data*     twist;
  alcap_data*     alcap;

  TGraphErrors*   gr_fit_pe;
  TGraphErrors*   gr_fit_de;

  TF1*            f_prot;
  TF1*            f_deut;

  muon_capture();

  static double fitf(double* X, double* P);

  void   init_data();
  void   init_proton_fit ();
  void   init_deuteron_fit ();

  void   fit_proton_energy_spectrum();

  void   fit_deuteron_energy_spectrum();

  void   plot_energy_spectra();

  void   plot(int Figure, int Print);
};

muon_capture*  mc;
twist_data*    twist;
alcap_data*    alcap;

//-----------------------------------------------------------------------------
muon_capture::muon_capture() {

  alcap = nullptr;
  twist = nullptr;

  init_data();
  
  f_prot = new TF1("f_prot",fitf,0,50,7);
  f_deut = new TF1("f_deut",fitf,0,50,7);

  f_prot->SetLineColor(kBlue+2);
  init_proton_fit();

  f_deut->SetLineColor(kBlue+2);
  init_deuteron_fit();
}

//-----------------------------------------------------------------------------
double muon_capture::fitf(double* X, double* P) {
  double f;

  double e  = X[0];
  double e0 =  P[0];
  double e1 =  P[6];
  double e2 = 10;

  if (e <= e0) {
    f = 0;
  }
  else if (e < e1) {
    f = pow(1.-e0/e,P[2])*exp(-e/P[3]);
  }
  else if (e < e2) {
    double c4 = pow(1.-e0/e1,P[2])*exp(-e1/P[3]);
    f = c4*exp(-(e-e1)/P[4]);
  }
  else {
    double c6 = pow(1.-e0/e1,P[2])*exp(-e1/P[3])*exp(-(e2-e1)/P[4]);
    f = c6*exp(-(e-e2)/P[5]);
  }

  return P[1]*f;
}

//-----------------------------------------------------------------------------
void muon_capture::init_proton_fit() {

  f_prot->FixParameter(0,1.4);
  f_prot->SetParameter(1,0.025);
  f_prot->SetParameter(2,1.1);
  f_prot->SetParameter(3,4.5);
  f_prot->SetParameter(4,5.);
  f_prot->SetParameter(5,6.);
  f_prot->SetParameter(6,3.);
}

//-----------------------------------------------------------------------------
void muon_capture::init_deuteron_fit() {

  f_deut->FixParameter(0,1.4);
  f_deut->SetParameter(1,0.01);
  f_deut->SetParameter(2,1.3);
  f_deut->SetParameter(3,4);
  f_deut->SetParameter(4,5.);
  f_deut->SetParameter(5,8.);
  f_deut->SetParameter(6,5.);
}

//-----------------------------------------------------------------------------
void muon_capture::init_data() {
  if (! gInterpreter->IsLoaded("muon_capture/scripts/alcap_data.C")) {
    gInterpreter->LoadMacro("muon_capture/scripts/alcap_data.C");
  }

  if (! gInterpreter->IsLoaded("muon_capture/scripts/twist_data.C")) {
    gInterpreter->LoadMacro("muon_capture/scripts/twist_data.C");
  }

  gInterpreter->ProcessLine("twist = new twist_data();");
  gInterpreter->ProcessLine("alcap = new alcap_data();");
  printf("initialized\n");
}

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

  printf(" nptot = %3i\n",nptot);
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

    printf("i,x,y,ey : %3i %12.5e %12.5e %12.5e\n",i,x[i],y[i],ey[i]);
  }

  printf("npp = %3i\n",npp);

  gr_fit_pe = new TGraphErrors(nptot,x,y,ex,ey);
  gr_fit_pe->SetName("gr_fit_pe");
  gr_fit_pe->SetTitle("proton energy");

  TCanvas* c_fit = new TCanvas("c_fit","c fit",1300,900);

  c_fit->cd();

  gr_fit_pe->SetMarkerStyle(20);
  gr_fit_pe->SetMarkerSize(1);
  gr_fit_pe->SetMarkerColor(kRed+2);
  gr_fit_pe->SetLineColor  (kRed+2);
  gr_fit_pe->SetFillStyle(3003);
  gr_fit_pe->SetFillColor(kRed+2);

  gr_fit_pe->Draw("ap");

  gr_fit_pe->Fit(f_prot,"","");
}

//-----------------------------------------------------------------------------
// AlCap doesn't have deuteron data
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

  TCanvas* c_fit_de = new TCanvas("c_fit_de","c fit",1300,900);

  c_fit_de->cd();

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

  double protons_per_capture(0.05);

  TCanvas* c_e = new TCanvas("c_e","c_e",1300,800);
  
  plot_ejected_proton_spectrum("e",protons_per_capture);
  twist->pe.gr->Draw("same,l3");
  twist->de.gr->Draw("same,l3");
  alcap->pe.gr->Draw("same,l3");

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry("ep_fun"   ,"Mu2e protons"   ,"f");
  leg->AddEntry(twist->pe.gr,"TWIST protons"  ,"f");
  leg->AddEntry(twist->de.gr,"TWIST deuterons","f");
  leg->AddEntry(alcap->pe.gr,"AlCap protons"  ,"f");

  leg->Draw();
//-----------------------------------------------------------------------------
// plot momentum distributions
//-----------------------------------------------------------------------------
  TCanvas* c_p = new TCanvas("c_p","c_p",1300,800);
  
  plot_ejected_proton_spectrum("p",protons_per_capture);
  twist->pp.gr->Draw("same,l3");
  twist->dp.gr->Draw("same,l3");
  alcap->pp.gr->Draw("same,l3");

  TLegend* leg2 = new TLegend(0.6,0.6,0.9,0.9);
  leg2->AddEntry("ep_fun"   ,"Mu2e protons"   ,"f");
  leg2->AddEntry(twist->pp.gr,"TWIST protons"  ,"f");
  leg2->AddEntry(twist->dp.gr,"TWIST deuterons","f");
  leg2->AddEntry(alcap->pp.gr,"AlCap protons"  ,"f");

  leg2->Draw();
}


//-----------------------------------------------------------------------------
void muon_capture::plot(int Figure, int Print) {

  TString canvas_name = Form("Figure_%04i",Figure);
  TCanvas* c = new TCanvas(canvas_name,canvas_name,1300,900);

  if (Figure == 1) {
    twist->pp.gr->Draw("alp");
    twist->pp.gr->Draw("p3,same");
    alcap->pp.gr->Draw("p3,same");
  }

}
