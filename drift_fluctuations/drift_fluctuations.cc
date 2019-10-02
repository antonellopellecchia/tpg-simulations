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
#include "Garfield/ViewDrift.hh"
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
	gas->LoadGasFile("ar_70_co2_30.gas");
	gas->PrintGas();

	// Model TPG drift gap as a parallelepiped.
	// Sides [cm]
	const double length = 12;
	const double width = 1.2;
	const double depth = 4;
	SolidBox* box = new SolidBox(length/2., width/2., depth/2., length/2., width/2., depth/2.);

	// Combine gas and box to a simple geometry.
	GeometrySimple *geo = new GeometrySimple();
	geo->AddSolid(box, gas);
	// View geometry
	ViewGeometry *viewGeometry = new ViewGeometry();
	viewGeometry->SetGeometry(geo);
	viewGeometry->Plot();

	// Make a component with constant electric field.
	// 400 V/cm along down z axis
	ComponentConstant *field = new ComponentConstant();
	field->SetGeometry(geo);
	field->SetElectricField(0., 0., -400.);
	// View electric field
	ViewField *view = new ViewField();
	view->SetComponent(field);
	//view->PlotSurface("v");

	Sensor *sensor = new Sensor();
	sensor->AddComponent(field);

	float w0 = 23.42; // gaussian beam waist
	TF1 *gaussianIntensity = new TF1("gauss_intensity", "TMath::Exp(-x**2/[0]**2)", -100, 100); // gaussian beam intensity distribution
	gaussianIntensity->SetParameter(0, w0);

	const int nBins = 150;
	TH1::StatOverflows(true);
	TH1F *hDriftTime = new TH1F("hDriftTime", "Electron drift times", nBins, 2000, 2300);
	const int nEvents = 1000;
	for (unsigned int i = 0; i < nEvents; ++i) {
		// Initial coordinates of the electron.
		double x0 = length/2.;
		double y0 = width/2.;
		double z0 = depth/2.+gaussianIntensity->GetRandom()*1e-4;
		double t0 = .1;
		double x1, y1, z1, t1; // final position and time
		int status = 0;

		AvalancheMC *drift = new AvalancheMC();
		drift->SetSensor(sensor);
		drift->DriftElectron(x0, y0, z0, t0);

		//drift->GetEndPoint(x1, y1, z1, t1, status);
		drift->GetElectronEndpoint(0, x0, y0, z0, t0, x1, y1, z1, t1, status);
		//cout << "Time spread: " << drift->GetArrivalTimeSpread() << endl;
		hDriftTime->Fill(t1); // time in ns
		if (i % 10 == 0) cout << i << "/" << nEvents << ", " << x1 << " " << y1 << " " << z1 << " " << t1 << endl;
	}

	AvalancheMC *drift = new AvalancheMC();
	ViewDrift *viewDrift = new ViewDrift();
	drift->SetSensor(sensor);
	drift->EnablePlotting(viewDrift);
	drift->DriftElectron(length/2., width/2., depth/2., 0);
	viewDrift->Plot();

	cout << "Mean " << hDriftTime->GetMean() << ", RMS " << hDriftTime->GetRMS() << endl;
	TCanvas *c = new TCanvas("c", "Drift time distribution", 600, 600);
	c->cd();
	hDriftTime->SetFillColor(kBlue + 2);
	hDriftTime->SetLineColor(kBlue + 2);
	hDriftTime->GetXaxis()->SetTitle("Drift time");
	hDriftTime->SetTitle("Distribution of electron drift times");
	hDriftTime->Draw();
	hDriftTime->Draw();

	char pathstr[200];
	sprintf(pathstr, "./results/drift_fluctuations%d.root", nEvents);
	c->SaveAs(pathstr);
	sprintf(pathstr, "./results/drift_fluctuations%d.png", nEvents);
	c->SaveAs(pathstr);

	cout << "Done." << endl;

	app.Run(true);
}
