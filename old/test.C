// X-ray conversion
// -------------------------------------------------------------------
// Simulate the Fe55 spectrum
//
// Thanks to Dorothea Pfeiffer and Heinrich Schindler from CERN for their help.
// -------------------------------------------------------------------
// Lucian Scharenberg
// scharenberg@physik.uni-bonn.de
// 05 APR 2018

#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

void primaries() {
	TFile *f1 = new TFile("../data/Xray-mca/30_04_2019/1_2mm_coll_25um_Cu.root");
	gDirectory->ls();
	TH1F *spectrum = (TH1F*)f1->Get("spectrum");
	//spectrum->Draw();

	TH1F *histo = new TH1F("histo", "histogram", spectrum->GetNbinsX(), spectrum->GetBinCenter(0), spectrum->GetBinCenter(spectrum->GetNbinsX()));
	for (int i=0; i<1000000; i++) {
		float c = spectrum->GetRandom();
		histo->Fill(c);
	}
	histo->Draw();
}


