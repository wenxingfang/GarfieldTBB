#include <iostream>

#include <TApplication.h>

#include "Garfield/SolidWire.hh" 
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include <TCanvas.h>
#include <TROOT.h>

using namespace Garfield;

int main(int argc, char * argv[]) {

  //TApplication app("app", &argc, argv);

  MediumMagboltz gas;
  gas.SetComposition("he", 50.,"isobutane", 50.);  // cepc gas
  gas.SetTemperature(293.15);
  gas.SetPressure(760.0);
  gas.SetMaxElectronEnergy(200.);
  gas.EnablePenningTransfer(0.44, 0.0, "He");
  gas.LoadGasFile("He_50_isobutane_50.gas");
  gas.LoadIonMobility("/workfs/higgs/fangwx/fork_Tracker_20201124/CEPCSW/garfieldpp/Data/IonMobility_He+_He.txt");

  MediumConductor metal;

  // Geometry.
  GeometrySimple geo;
  // Field Wire radius [cm]
  const double rFWire = 110.e-4;
  // Signa Wire radius [cm]
  const double rSWire = 25.e-4;
  const double halflength = 200;
  // Cell radius [cm]
  const double m_rCell = 0.5;
  // Voltages
  const double vSWire = 2000.;
  const double vFWire = 0.;
  SolidWire* wire1 = new SolidWire(0, 0, 0, rSWire, halflength, 0, 0, 1);
  wire1->SetBoundaryPotential(vSWire);
  wire1->SetLabel("s1");
  SolidWire* wire2 = new SolidWire(-m_rCell, -m_rCell, 0, rFWire, halflength, 0, 0, 1);
  wire2->SetBoundaryPotential(vFWire);
  SolidWire* wire3 = new SolidWire(0       , -m_rCell, 0, rFWire, halflength, 0, 0, 1);
  wire3->SetBoundaryPotential(vFWire);
  SolidWire* wire4 = new SolidWire(m_rCell , -m_rCell, 0, rFWire, halflength, 0, 0, 1);
  wire4->SetBoundaryPotential(vFWire);
  SolidWire* wire5 = new SolidWire(-m_rCell, 0       , 0, rFWire, halflength, 0, 0, 1);
  wire5->SetBoundaryPotential(vFWire);
  SolidWire* wire6 = new SolidWire( m_rCell, 0       , 0, rFWire, halflength, 0, 0, 1);
  wire6->SetBoundaryPotential(vFWire);
  SolidWire* wire7 = new SolidWire(-m_rCell, m_rCell , 0, rFWire, halflength, 0, 0, 1);
  wire7->SetBoundaryPotential(vFWire);
  SolidWire* wire8 = new SolidWire( m_rCell, m_rCell , 0, rFWire, halflength, 0, 0, 1);
  wire8->SetBoundaryPotential(vFWire);
  SolidWire* wire9 = new SolidWire( 0      , m_rCell , 0, rFWire, halflength, 0, 0, 1);
  wire9->SetBoundaryPotential(vFWire);

  geo.AddSolid(wire1, &metal);
  geo.AddSolid(wire2, &metal);
  geo.AddSolid(wire3, &metal);
  geo.AddSolid(wire4, &metal);
  geo.AddSolid(wire5, &metal);
  geo.AddSolid(wire6, &metal);
  geo.AddSolid(wire7, &metal);
  geo.AddSolid(wire8, &metal);
  geo.AddSolid(wire9, &metal);
  geo.SetMedium(&gas);

  
  
  //TCanvas* cD = nullptr;
  //cD = new TCanvas("DCH_geomView2d", "", 600, 600);
  ViewGeometry geomView2d;
  //geomView2d.SetCanvas(cD);
  geomView2d.SetGeometry(&geo);
  geomView2d.SetArea(-1, -1, -1, 1, 1, 2);
  geomView2d.SetPlane(0, 0, 1, 0, 0, 0);
  geomView2d.Plot2d();
  //cD->SaveAs("DCH_geomView2d.png");
  //delete cD;
  /*
  //cD = new TCanvas("DCH_geomView3d", "", 600, 600);
  ViewGeometry geomView3d;
  //geomView3d.SetCanvas(cD);
  geomView3d.SetGeometry(&geo);
  geomView3d.Plot();
  //cD->SaveAs("DCH_geomView3d.png");
  //delete cD;
  */
   
  ComponentNeBem3d nebem;

  nebem.SetGeometry(&geo);
  //nebem.SetTargetElementSize(0.01);
  nebem.UseSVDInversion();
  nebem.Initialise();
  nebem.SetMagneticField(0., 0., -3.);
  /*
  cD = new TCanvas(str1.c_str(), "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(cD);
  fieldView.SetComponent(&m_nebem);
  fieldView.SetArea(-1, -1, -1, 1, 1, 1);
  fieldView.SetPlane(0, 0, 1, 0, 0, 0  );
  fieldView.PlotContour("e");
  str0 = "/fieldView_ContourE_";
  cD->SaveAs((m_plot_path+str0+str1+str2).c_str());
  delete cD;
  */
  ///
  /// Make a sensor.
  ///
  Sensor sensor;
  sensor.AddComponent(&nebem);
  sensor.AddElectrode(&nebem, "s1");

  // Set the signal time window. [ns]
  const double tstep = 0.5;
  const double tmin = -0.5 * 0.5;
  const unsigned int nbins = 1000;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  sensor.ClearSignal();

  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  aval.EnableSignalCalculation();
  aval.EnableInducedChargeCalculation();
  aval.EnableMagneticField();
  //aval.EnableDebugging();
  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);
  drift.EnableSignalCalculation();
  drift.EnableInducedChargeCalculation();
  drift.EnableMagneticField();
  /*    
  TCanvas* cD = nullptr;
  ViewDrift driftView;
  if (m_plotDrift) {
    cD = new TCanvas(str1.c_str(), "", 600, 600);
    driftView.SetCanvas(cD);
    aval.EnablePlotting(&driftView);
  }
  */
  int N_primary = 50;     
  for(int n=0; n<N_primary; n++){
      double xe = 0.1, ye = 0.1, ze = 0., te = 0., ee = 0.;
      xe = 0.5*m_rCell;
      ye = -m_rCell+ (n+1)*(1.5*m_rCell/N_primary);
      double dx = 0., dy = 0., dz = 0.;
      aval.AvalancheElectron(xe, ye, ze, te, ee, dx, dy, dz);
      
      int ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      const int np = aval.GetNumberOfElectronEndpoints();
      std::cout << "ne="<<ne<<",ni="<<ni<<" np=" << np << std::endl;

      double xe1, ye1, ze1, te1, e1;
      double xe2, ye2, ze2, te2, e2;
      int status;
      for (int k = 0; k < np; k++)
      {
      //    info() << "k="<<k<<endmsg;
          aval.GetElectronEndpoint(k, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
          drift.DriftIon(xe1, ye1, ze1, te1);
      }
  }


  //app.Run(true);
}


