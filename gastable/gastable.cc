#include <iostream>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {
  const double pressure = AtmosphericPressure;
  const double temperature = 293.15;
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("Ar", 70., "CO2", 30.);
 
  // Set the field range to be covered by the gas table. 
  const int nFields = 40;
  const double emin = 100.;
  const double emax = 600.;
  // Flag to request logarithmic spacing.
  const bool useLog = false;
  gas->SetFieldGrid(emin, emax, nFields, useLog); 

  cout << "Running gas table calculation..." << endl;
  const int ncoll = 10;
  gas->GenerateGasTable(ncoll);
  cout << "Saving gas table on disk..." << endl;
  gas->WriteGasFile("ar_70_co2_30_2019_10_01.gas");
  cout << "Done." << endl;
}
