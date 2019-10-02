// Drift time
// --------------------------------
// Simulate the distribution of the drift times
// of the primary electrons in a parallelepipedal box
// (TPG drift gap)

#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1F.h>

#include "Garfield/AvalancheMC.hh"
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

	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();

	// Make a gas medium.
	MediumMagboltz* gas = new MediumMagboltz();
	gas->SetComposition("Ar", 70., "CO2", 30.);
	gas->SetTemperature(293.15);
	gas->SetPressure(AtmosphericPressure);

	// Model TPG drift gap as a parallelepiped.
	// Sides [cm]
	const double length = 6;
	const double width = 1;
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
	// View electric field
	ViewField* view = new ViewField();
	view->SetComponent(&field);
	view->PlotSurface("ez");

	// Histogram
	const int nBins = 1500;
	TH1::StatOverflows(true);
	TH1F hDriftTime("hDriftTime", "Electron drift times", nBins, 1e8, 2e8);
	const int nEvents = 1;
	for (unsigned int i = 0; i < nEvents; ++i) {
		if (i % 1000 == 0) std::cout << i << "/" << nEvents << "\n";
		// Initial coordinates of the electron.
		double x0 = length/2.;
		double y0 = width/2.;
		double z0 = depth/2.;
		double t0 = 0.;
		double x1, y1, z1, t1; // final position and time
		int status = 0;

		Sensor sensor;
		sensor.AddComponent(&field);

		AvalancheMC *avalanche = new AvalancheMC();
		//avalanche->SetSensor(&sensor);
		avalanche->DriftElectron(x0, y0, z0, t0);

		avalanche->GetElectronEndpoint(0,
			x0, y0, z0, t0,
			x1, y1, z1, t1,
			status);
		hDriftTime.Fill(t1);
	}

	cout << "Mean " << hDriftTime.GetMean() << ", RMS " << hDriftTime.GetRMS() << endl;
	TCanvas c("c", "Drift time distribution", 600, 600);
	c.cd();
	hDriftTime.SetFillColor(kBlue + 2);
	hDriftTime.SetLineColor(kBlue + 2);
	hDriftTime.GetXaxis()->SetTitle("Drift time");
	hDriftTime.SetTitle("Distribution of electron drift times");
	//hDriftTime.Draw();

	char pathstr[21];
	sprintf(pathstr, "./results/drift_fluctuations%d.root", nEvents);
	c.SaveAs(pathstr);
	sprintf(pathstr, "./results/drift_fluctuations%d.png", nEvents);
	c.SaveAs(pathstr);

	cout << "Done." << endl;

	app.Run(true);
}