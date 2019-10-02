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
#include <TFile.h>
#include <TH1F.h>

#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

	//TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();

	// Make a gas medium.
	MediumMagboltz* gas = new MediumMagboltz();
	gas->SetComposition("Ar", 70., "CO2", 30.);
	gas->SetTemperature(293.15);
	gas->SetPressure(AtmosphericPressure);

	// Model TPG drift gap as a parallelepiped.
	// Sides [cm]
	const double length = 10;
	const double width = 10;
	const double depth = 4;
	SolidBox box(length/2., width/2., depth/2., length/2., width/2., depth/2.);

	// Combine gas and box to a simple geometry.
	GeometrySimple geo;
	geo.AddSolid(&box, gas);
	// View geometry
	//ViewGeometry* viewGeometry = new ViewGeometry();
	//viewGeometry->SetGeometry(&geo);
	//viewGeometry->Plot();

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
	const int nBins = 1500;
	TH1::StatOverflows(true);
	TH1F hElectrons("hElectrons", "Number of electrons", nBins, -0.5, nBins - 0.5);
	const int nEvents = 1000000;
	for (unsigned int i = 0; i < nEvents; ++i) {
		if (i % 1000 == 0) std::cout << i << "/" << nEvents << "\n";
		// Initial coordinates of the photon.
		const double x0 = length/2.;
		const double y0 = 0.;
		const double z0 = depth/2.;
		const double t0 = 0.;
		const double egamma = spectrum->GetRandom()*1e3;
		int ne = 0;
		track.TransportPhoton(x0, y0, z0, t0, egamma, 0., 1., 0., ne);
		if (ne > 0) hElectrons.Fill(ne);
	}

	cout << "Mean " << hElectrons.GetMean() << ", RMS " << hElectrons.GetRMS() << endl;
	TCanvas c("c", "", 600, 600);
	c.cd();
	hElectrons.SetFillColor(kBlue + 2);
	hElectrons.SetLineColor(kBlue + 2);
	hElectrons.GetXaxis()->SetTitle("Electron number");
	hElectrons.SetTitle("Distribution of primary electrons");
	hElectrons.Draw();

	char pathstr[100];
	sprintf(pathstr, "./images/primaries_single_70Ar_30CO2_%devents.png", nEvents);
	c.SaveAs(pathstr);
	sprintf(pathstr, "./data/primaries_single_70Ar_30CO2_%devents.root", nEvents);
	hElectrons.SaveAs(pathstr);

	//app.Run(true);
}
