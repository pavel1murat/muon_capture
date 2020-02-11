///////////////////////////////////////////////////////////////////////////////
// prerequisites:
// --------------
// .L muon_capture/scripts/muon_capture_yields.C
//
// mc = new muon_capture();
// mc->fit_proton_energy_spectrum()
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __ana_muon_capture__
#define __ana_muon_capture__

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

  TGraphErrors*   gr_fit_pe;
  TGraphErrors*   gr_fit_de;

  TF1*            f_prot;
  TF1*            f_prot_2;
  TF1*            f_deut;

  TF1*            fp_hgf;  // https://arxiv.org/pdf/1803.08403.pdf
  TF1*            fe_hgf;  // 

  muon_capture();

  static double fitf_prot  (double* X, double* P);
  static double fitf_prot_2(double* X, double* P);

  void   init_data        ();

  void   init_proton_fit  ();
  void   init_proton_fit_2();
  void   init_deuteron_fit();
  void   init_hgf_spectrum();

  void   fit_proton_energy_spectrum();

  void   fit_deuteron_energy_spectrum();

  static double hgf_spectrum(double* X, double* Par);

  void   plot_energy_spectra();

  void   plot(int Figure, int Print = 0);
};

extern muon_capture*  mc;

#endif
