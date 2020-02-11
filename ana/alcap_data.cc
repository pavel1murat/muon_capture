///////////////////////////////////////////////////////////////////////////////
// alcap data - in energy, keV, plotted dN/d(0.5 MeV)
///////////////////////////////////////////////////////////////////////////////
#include "ana/alcap_data.hh"

alcap_data* alcap(nullptr);

//-----------------------------------------------------------------------------
alcap_data::alcap_data() {
  MP = 938.;
  init();
}

//-----------------------------------------------------------------------------
alcap_data::~alcap_data() {
}

//-----------------------------------------------------------------------------
void alcap_data::init() {

  float data [] = {
		   2247.79,0.00280366,
		   2769.91,0.00277421,
		   3265.49,0.00262696,
		   3769.91,0.00272120,
		   4274.34,0.00260340,
		   4761.06,0.00221466,
		   5274.34,0.00194372,
		   5761.06,0.00182003,
		   6265.49,0.00183770,
		   6761.06,0.00149018,
		   7247.79,0.00128403,
		   7752.21,0.00112500,
		   8247.79,0.00112500,
		   8761.06,0.000912958,
		   9256.64,0.000842277,
		   9761.06,0.000848168,
		   10256.6,0.000783377,
		   10752.2,0.000689136,
		   11247.8,0.000618456,
		   11761.1,0.000636126,
		   -1
  };

  float data_upper[] = {
			2156.64,0.00354581,
			2256.64,0.00354581,
			2699.11,0.00338089,
			3274.34,0.00294503,
			3752.21,0.00289202,
			4274.34,0.00272120,
			4814.16,0.00225589,
			5274.34,0.00200262,
			5769.91,0.00189071,
			6256.64,0.00190249,
			6752.21,0.00157853,
			7238.94,0.00134293,
			7752.21,0.00118979,
			8247.79,0.00117212,
			8761.06,0.000960079,
			9238.94,0.000889398,
			9743.36,0.000895288,
			10265.5,0.000830497,
			10743.4,0.000742147,
			11230.1,0.000659686,
			11761.1,0.000689136,
			12266.1,0.000719136,
			-1
  };

  float data_lower[] = {
			2147.79,0.00206741,
			2247.79,0.00206741,
			2769.91,0.00218521,
			3265.49,0.00228534,
			3769.91,0.00253272,
			4247.79,0.00247382,
			4761.06,0.00213809,
			5265.49,0.00187893,
			5778.76,0.00176702,
			6265.49,0.00175524,
			6778.76,0.00141950,
			7265.49,0.00121335,
			7761.06,0.00107199,
			8256.64,0.00106021,
			8761.06,0.000848168,
			9256.64,0.000783377,
			9752.21,0.000777487,
			10247.8,0.000712696,
			10752.2,0.000636126,
			11256.6,0.000571335,
			11761.1,0.000583115,
			12266.1,0.000598115,
			-1
  };

  float x [1000], y [1000], ex [1000], ey [1000]; 
  float xe[1000], ye[1000], exe[1000], eye[1000]; 

  int npt = 0;
  for(int i=0; data_upper[2*i] >= 0; i++) {
    float e = data_upper[2*i  ]/1000;
    float p = sqrt(2*MP*e);

    x[npt]  = p;
    y[npt]  = data_upper[2*i+1]*2*p/MP;   // plotted - dN/dE

    xe[npt] = e;
    ye[npt] = data_upper[2*i+1]*2;
    npt++;
  }

  pp.gru = new TGraph(npt,x,y);
  pp.gru->SetName("alcap_pp_gru");
  pp.gru->SetLineColor(kBlue+2);
  pp.su = new smooth(pp.gru,x[0],x[npt-1]);

  pe.gru = new TGraph(npt,xe,ye);
  pe.gru->SetName("alcap_pe_gru");
  pe.gru->SetLineColor(kBlue+2);
  pe.su = new smooth(pe.gru,x[0],x[npt-1]);

  npt = 0;
  for(int i=0; data_lower[2*i] >= 0; i++) {
    float e = data_lower[2*i  ]/1000;
    float p = sqrt(2*MP*e);
    
    x[npt]  = p;
    y[npt]  = data_lower[2*i+1]*2*p/MP;
    
    xe[npt] = e;
    ye[npt] = data_lower[2*i+1]*2;
    
    npt++;
  }

  pp.grl = new TGraph(npt,x,y);
  pp.grl->SetName("alcap_pp_grl");
  pp.grl->SetLineColor(kBlue+2);
  pp.sl = new smooth(pp.grl,x[0],x[npt-1]);

  pe.grl = new TGraph(npt,xe,ye);
  pe.grl->SetName("alcap_pe_grl");
  pe.grl->SetLineColor(kBlue+2);
  pe.sl = new smooth(pe.grl,x[0],x[npt-1]);

  npt = 0;

  //  printf(" ----- data\n");
  for(int i=0; data[2*i] >= 0; i++) {
    float e = data[2*i  ]/1000;
    float p = sqrt(2*MP*e);
    
    x[npt]    = p;
    y[npt]    = data[2*i+1]*2*p/MP;

    float eu  = pp.su->GetFunc()->Eval(p);
    float el  = pp.sl->GetFunc()->Eval(p);

    float err = (eu-el)/2;
    ex[npt]   = 0;
    ey[npt]   = err;
    
    xe[npt]    = e;
    ye[npt]    = data[2*i+1]*2;

    eu        = pe.su->GetFunc()->Eval(e);
    el        = pe.sl->GetFunc()->Eval(e);

    err       = (eu-el)/2;
    exe[npt]  = 0;
    eye[npt]  = err;

    //    printf("%5i %12.5e %12.5e %12.5e %12.5e %12.5e\n",npt,x[npt],y[npt],eu,el,ey[npt]);
    npt++;
  }

  pp.gr = new TGraphErrors(npt,x,y,ex,ey);
  pp.gr->SetName("alcap_pp_gr");
  pp.gr->SetFillColor(kBlue+2);
  pp.gr->SetFillStyle(3002);
  pp.gr->SetLineColor(kBlue+2);

  pp.gr->SetMarkerColor(kBlue+2);
  pp.gr->SetMarkerStyle(20);
  pp.gr->SetMarkerSize(1.0);

  pp.s = new smooth(pp.gr,x[0],x[npt-1]);

  pe.gr = new TGraphErrors(npt,xe,ye,exe,eye);
  pe.gr->SetName("alcap_pe_gr");
  pe.gr->SetFillColor(kBlue+2);
  pe.gr->SetFillStyle(3002);
  pe.gr->SetLineColor(kBlue+2);

  pe.gr->SetMarkerColor(kBlue+2);
  pe.gr->SetMarkerStyle(20);
  pe.gr->SetMarkerSize(1.0);

  pe.s = new smooth(pe.gr,xe[0],xe[npt-1]);

  //  alcap.gr->Draw("A3");
}		
