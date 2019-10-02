#include <iostream>
//#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "ViewMedium.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {
  TApplication app("app", &argc, argv);

  cout << "Running" << endl;
  
  MediumMagboltz *gas = new MediumMagboltz();
  gas->LoadGasFile("ar_70_co2_30_2019_10_01.gas");
  gas->PrintGas();
  
  ViewMedium *view = new ViewMedium();
  view->SetMedium(gas);
  
  TCanvas *cV = new TCanvas("cV", "", 600, 600);
  view->SetCanvas(cV);
  view->PlotElectronVelocity('e', 400., 0., 0.);
  cV->SaveAs("./electron_velocity.eps");  

  //app.Run(kTRUE);
}
