// X-ray conversion
// -------------------------------------------------------------------
// Simulate the Fe55 spectrum
//
// Thanks to Dorothea Pfeiffer and Heinrich Schindler from CERN for their help.
// -------------------------------------------------------------------
// Lucian Scharenberg
// scharenberg@physik.uni-bonn.de
// 05 APR 2018

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1F.h>

#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

	if (argc < 5) {
		cout << "Usage: primaries [temperature (K)] [argon fraction (%)] [co2 fraction (%)] [rate] [histogram lower limit] [histogram upper limit]" << endl;
		return -1;
	}

	//TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();

	const float temperature = atof(argv[1]);
	const char *gas1 = argv[2];
	const char *gas2 = argv[3];
	const float gas1Fraction = atof(argv[4]);
	const float gas2Fraction = atof(argv[5]);
	const float rate = atof(argv[6]);
	const float bin_low = atof(argv[7]);
	const float bin_high = atof(argv[8]);
	cout << "Temperature (K): " << temperature << endl;
	cout << gas1 << " (%): " << gas1Fraction << endl;
	cout << gas2 << " (%): " << gas2Fraction << endl;
	cout << "Rate: " << rate << endl;

	// Make a gas medium.
	MediumMagboltz* gas = new MediumMagboltz();
	gas->SetComposition(gas1, gas1Fraction, gas2, gas2Fraction);
	gas->SetTemperature(temperature);
	gas->SetPressure(AtmosphericPressure);

	// Model TPG drift gap as a parallelepiped.
	// Sides [cm]
	const double length = 20;
	const double width = 20;
	const double depth = 4;
	SolidBox box(length/2., width/2., depth/2., length/2., width/2., depth/2.);

	// Combine gas and box to a simple geometry.
	GeometrySimple geo;
	geo.AddSolid(&box, gas);

	// Make a component with constant electric field.
	// 400 V/cm along down z axis
	ComponentConstant field;
	field.SetGeometry(&geo);
	field.SetElectricField(0., 0., -400.);

	// Make a sensor.
	Sensor sensor;
	sensor.AddComponent(&field);

	// Get the source spectrum
	TFile *f1 = new TFile("data/1_2mm_coll_25um_Cu.root");
	gDirectory->ls();
	TH1F *spectrum = (TH1F*)f1->Get("spectrum");

	// Use Heed for simulating the photon absorption.
	TrackHeed track;
	track.SetSensor(&sensor);
	track.EnableElectricField();

	// Histogram
	TH1::StatOverflows(true);
	const int nBins = 50;
	TH1F hTotalCharge("hTotalCharge", "Total primary charge", nBins, bin_low, bin_high);
	//const float rate = 79124.433333;
	//const float rate = 74722.2;
	const int acquisitionTime = 10; // time over which the current is measured (seconds)
	const int nSamples = 500;
	const int nEvents =  (int) (rate * acquisitionTime); // 100000
	cout << nEvents << " total events" << endl;
	long timeStamp = time(0);
	for (int j=0; j < nSamples; j++) {
	        float speed = 1. / (time(0)-timeStamp);
		float remainingTime = (nSamples - j - 1)/speed;
		cout << j << "/" << nSamples << ", " << remainingTime << "s remaining" << "\n";
		timeStamp = time(0);
		int nPrimary = 0;
		cout << "\n";
		for (int i = 0; i < nEvents; ++i) {
		  if (i % 100 == 0) cout << "\r" << i << "/" << nEvents;
			// Initial coordinates of the photon.
			const double x0 = length/2.;
			const double y0 = 0.;
			const double z0 = depth/2.;
			const double t0 = 0.;
			double x, y, z, t, e, dx, dy, dz;
			int ne = 0;
			do {
				const double egamma = spectrum->GetRandom()*1e3;
				track.TransportPhoton(x0, y0, z0, t0, egamma, 0., 1., 0., ne);
				if (ne > 0) track.GetElectron(0, x, y, z, t, e, dx, dy, dz);
				//if (x<width/2.-0.5 || x>width/2.+0.5)
			} while (ne==0 || y<width/2.-0.5 || y>width/2.+0.5);
			nPrimary += ne;
		}
		cout << "\n";
		//for (int k=1; k < hElectrons.GetNbinsX(); k++) nPrimary += hElectrons.GetBinContent(k);
		hTotalCharge.Fill(nPrimary);
		//cout << nPrimary << endl;
	}
	cout << "Mean " << hTotalCharge.GetMean() << ", RMS " << hTotalCharge.GetRMS() << endl;

	TCanvas c("c", "", 600, 600);
	c.cd();
	hTotalCharge.SetFillColor(kGreen + 2);
	hTotalCharge.SetLineColor(kGreen + 2);
	hTotalCharge.GetXaxis()->SetTitle("Total primary electrons");
	hTotalCharge.Draw();

	char pathstr[100];
	sprintf(pathstr, "./images/primaries_%drate_%d%s_%d%s.png", (int) rate, (int) gas1Fraction, gas1, (int) gas2Fraction, gas2);
	c.SaveAs(pathstr);
	sprintf(pathstr, "./data/primaries_%drate_%d%s_%d%s.root", (int) rate, (int) gas1Fraction, gas1, (int) gas2Fraction, gas2);
	hTotalCharge.SaveAs(pathstr);

	//app.Run(false);
}
