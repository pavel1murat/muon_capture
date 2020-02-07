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

  double          MP;
  double          MD;

  twist_data();
  ~twist_data();

  void init();

  void init_deuteron_data            ();
  void init_proton_data              ();
  void init_wh_data                  ();
  void init_proton_precompound_data  ();
  void init_deuteron_precompound_data();
  
  void init_proton_dog               ();   // data-over-geant
  void init_deuteron_dog             ();   // data-over-geant

  
};

#endif
