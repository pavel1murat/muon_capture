#ifndef __twist_data_hh__
#define __twist_data_hh__
///////////////////////////////////////////////////////////////////////////////
// TWIST data - momentum in MeV/c
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/alg/smooth.hh"
#include "TGraph.h"
#include "TGraphErrors.h"

struct twist_data {

  struct plot_data_1 {
    TGraph*       gr;
    smooth*       s;
  };

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
  
  plot_data       dp;
  plot_data       de;

  TGraph*         wh;    // supposedly, Ed's function
  smooth*         swh;

  plot_data_1     ppg;   // proton momentum   (G4 precompound)
  plot_data_1     peg;   // proton energy     (G4 precompound)

  plot_data_1     dpg;   // deuteron momentum (G4 precompound)
  plot_data_1     deg;   // deuteron energy   (G4 precompound)

  TGraph*         pp_dog;   // proton   momentum (P), data-over-geant
  TGraph*         dp_dog;   // deuteron momentum (P), data-over-geant
  
} twist;

void init_twist_data         ();

#endif
