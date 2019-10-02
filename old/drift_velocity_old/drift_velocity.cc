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
#include <TDatime.h>
#include <TGraphErrors.h>

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
	//ViewGeometry *viewGeometry = new ViewGeometry();
	//viewGeometry->SetGeometry(geo);
	//viewGeometry->Plot();

	// Vary the uniform electric field.
	// 100 V/cm along down z axis up to 600 V/cm
	const int fieldNumber = 5;
	TGraphErrors *velocityGraph = new TGraphErrors(fieldNumber);
	const int nEvents = 500;
	for (int fieldIndex = 0; fieldIndex < fieldNumber; fieldIndex++) {
		float electricField = 100. + fieldIndex*1000.;

		ComponentConstant *field = new ComponentConstant();
		field->SetGeometry(geo);
		field->SetElectricField(0., 0., electricField);
		// View electric field

		Sensor *sensor = new Sensor();
		sensor->AddComponent(field);

		float w0 = 23.42; // gaussian beam waist
		TF1 *gaussianIntensity = new TF1("gauss_intensity", "TMath::Exp(-x**2/2/[0]**2)", -100, 100); // gaussian beam intensity distribution
		gaussianIntensity->SetParameter(0, w0);

		const int nBins = 150;
		TH1::StatOverflows(true);
		TH1F *hDriftTime = new TH1F("hDriftTime", "Electron drift times", nBins, 2000, 2300);
		TH1F *hEndPositionZ = new TH1F("hEndPositionZ", "Electron end z", nBins, -0.005, 0.005);
		TH1F *hVelocity = new TH1F("hEndPositionZ", "Electron end z", nBins, -0.005, 0.005);

		cout << "Field: " << electricField << " V/cm" << endl;
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

			hDriftTime->Fill(t1-t0); // time in ns
			hEndPositionZ->Fill(z1-z0); // end position in cm
			hVelocity->Fill((z0-z1)/(t1-t0)); // average electron velocity
			if (i % 100 == 0) cout << i << "/" << nEvents << ", " << x1 << " " << y1 << " " << z1 << " " << t1 << "\n";
		}
		cout << endl;
		cout << "Times mean " << hDriftTime->GetMean() << ", RMS " << hDriftTime->GetRMS() << endl;
		cout << "Velocity mean " << hVelocity->GetMean() << ", RMS " << hVelocity->GetRMS() << endl << endl;

		velocityGraph->SetPoint(fieldIndex, electricField, hVelocity->GetMean());
		velocityGraph->SetPointError(fieldIndex, 0, hVelocity->GetRMS());
	}

	TCanvas *c = new TCanvas("c", "Drift velocity", 800, 600);
	c->cd();
	velocityGraph->GetXaxis()->SetTitle("Electric field (V/cm)");
	velocityGraph->GetYaxis()->SetTitle("Electron drift velocity (cm/ns)");
	velocityGraph->SetTitle("Distribution of electron drift times");

	velocityGraph->SetMarkerStyle(21);
	velocityGraph->SetMarkerColor(9);

	velocityGraph->Draw("AP");

	char pathstr[200];
	TDatime *datetime = new TDatime();
	sprintf(pathstr, "./results/velocity_%devents_%s.root", nEvents, datetime->AsString());
	c->SaveAs(pathstr);
	sprintf(pathstr, "./results/velocity_%devents_%s.pdf", nEvents, datetime->AsString());
	c->SaveAs(pathstr);

	cout << "Done." << endl;

	app.Run(true);
}
