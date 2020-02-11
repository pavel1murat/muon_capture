//

//-----------------------------------------------------------------------------
double fitf(double* X, double* P) {
  double f;


  double e  = X[0];
  double e0 = 1.4;  // P[0];
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


TF1* f1;
//-----------------------------------------------------------------------------
void fit() {

  f1 = new TF1("f1",fitf,0,30,7);

  f1->FixParameter(0,1.4);
  f1->SetParameter(1,0.025);
  f1->SetParameter(2,1.1);
  f1->SetParameter(3,4.5);
  //  f1->SetParameter(4,0.1);
  f1->SetParameter(4,5.);
  //  f1->SetParameter(5,0.035);
  f1->SetParameter(5,6.);

  f1->SetParameter(6,3.);
  
  f1->SetLineColor(kBlue+2);
}
