#ifndef __ana_alcap_data_hh__
#define __ana_alcap_data_hh__
///////////////////////////////////////////////////////////////////////////////
// alcap data - in energy, keV, plotted dN/d(0.5 MeV)
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/alg/smooth.hh"
#include "TGraph.h"
#include "TGraphErrors.h"

struct alcap_data {
  double MP;

  struct plot_data {
    TGraphErrors* gr;
    TGraph*       gru;
    TGraph*       grl;
    smooth*       s;
    smooth*       sl;
    smooth*       su;
  };

  plot_data       pp;
  plot_data       pe;

  alcap_data();
  ~alcap_data();

  void  init();
};

extern alcap_data* alcap;

#endif
