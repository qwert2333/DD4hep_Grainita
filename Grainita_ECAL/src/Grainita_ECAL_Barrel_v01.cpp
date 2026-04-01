//**************************************************************************
// \file    Grainita_ECAL_Barrel_v01.cpp
// \brief:  Grainita ECAL full detector geometry (first ver.)
// \author: Fangyi Guo (fangyi.guo@cern.ch) start from DD4HEP tutorial in 
//          DRD Calo (https://github.com/DRDCalo/DD4hepTutorials)
// \date:   March 2026
//**************************************************************************

// Includers from DD4hep
#include "DDRec/Vector3D.h"
#include <DD4hep/DetFactoryHelper.h>

using namespace dd4hep;

// Build geometry
//
static Ref_t create_detector(Detector &description, xml_h e,
                             SensitiveDetector sens) {
  std::cout << "--> Grainita_ECAL_Barrel_v01::create_detector() start" << std::endl;

  // Get info from the xml file
  //
  sens.setType("calorimeter");
  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();
  std::cout << "--> Going to create " << det_name << ", with ID: " << x_det.id() << std::endl;

  //xml_dim_t x_dim = x_det.dimensions();
  //const double CaloX = x_dim.x();
  //const double CaloY = x_dim.y();
  //const double CaloZ = x_dim.z();

  //Define readout
  //Readout readout = sens.readout();
  //Segmentation seg = readout.segmentation();


  //Define material
  Material  air(description.material("Air"));
  Material  MatCarbonfiber(description.material("CarbonFiber"));
  Material  MatCrystal(description.material( x_det.attr<std::string>(_U(material)) ));
  Material  MatWLSfiber(description.material("WLSFiber"));
  Material  MatPMMAcladding(description.material("PMMA"));


  // Detector global size
  auto innerR = description.constant<double>("InnerRadius");
  auto crystal_thick = description.constant<double>("CrystalThickness");
  auto back_space = description.constant<double>("BackSpace");
  auto halfZ = description.constant<double>("HalfZ");
  auto fiber_r = description.constant<double>("FiberRadius");
  auto cladding_thick = description.constant<double>("CladdingThickness");
  auto nSec_phi = description.constant<int>("NumberOfPhiSectors");
  auto nModule_z = description.constant<int>("NumberOfZModules");
  auto dPhi_fiber = description.constant<double>("FiberSeg_DPhi"); //Unit in degree
  auto dz_fiber = description.constant<double>("FiberSeg_Dz");
  auto Cframe_thick = description.constant<double>("FrameThickness");
  auto Cseg_thick = description.constant<double>("SegThickness");
  double outerR = innerR + crystal_thick + back_space + 2*Cframe_thick; 
  std::cout<<"  -- Dimension: inner R "<<innerR<<", outer R "<<outerR<<", half Z "<<halfZ<<std::endl;
  std::cout<<"     Z segmentation: "<<nModule_z<<", Phi segmentation: "<<nSec_phi<<std::endl;

  // Create the geometry
  DetElement ECAL(det_name, x_det.id());
  Volume worldVol = description.pickMotherVolume(ECAL);

  TGeoTube* EcalBarrel = new TGeoTube(innerR, outerR, halfZ);
  Volume EcalBarrelVol("EcalBarrel", EcalBarrel, air);

  // ====== Carbon fiber frame as supporting, at inner and outer of barrel. ====== //
  TGeoTube* innerFrame = new TGeoTube(innerR, innerR+Cframe_thick, halfZ);
  Volume innerFrame_vol("InnerCarbonFrame", innerFrame, MatCarbonfiber);
  TGeoTube* outerFrame = new TGeoTube(outerR-Cframe_thick, outerR, halfZ);
  Volume outerFrame_vol("outerCarbonFrame", outerFrame, MatCarbonfiber);

  //innerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  //outerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  //EcalBarrelVol.placeVolume(innerFrame_vol);
  //EcalBarrelVol.placeVolume(outerFrame_vol);

  
  // ======== Define a sector ======== //
  double dphi_sec = 360./nSec_phi; // In degree
  double half_dz_sec = halfZ/nModule_z; 
  double dphi_seg = Cseg_thick / innerR * (180./M_PI);
  TGeoTubeSeg* Sector = new TGeoTubeSeg(innerR+Cframe_thick, outerR-Cframe_thick, half_dz_sec, 0., dphi_sec);
  Volume Sector_vol("sector_vol", Sector, air);
  std::cout<<"  ---- For one sector: dPhi = "<<dphi_sec<<", dZ = "<<half_dz_sec<<std::endl;

  //    ------- Carbon fiber at 4 sides --------    //
  Assembly carbon_frame("carbon_frame");

  TGeoTubeSeg* frame_phi = new TGeoTubeSeg( innerR+Cframe_thick, 
                                            outerR-Cframe_thick, 
                                            half_dz_sec, 
                                            0, 
                                            dphi_seg); //Frame in phi side
  TGeoTubeSeg* frame_z = new TGeoTubeSeg( innerR+Cframe_thick, 
                                          outerR-Cframe_thick, 
                                          Cseg_thick/2., 
                                          dphi_seg, 
                                          dphi_sec-dphi_seg); //Frame in z side

  Volume frame_phi_vol("CarbonFramePhiMin", frame_phi, MatCarbonfiber);
  Volume frame_z_vol("CarbonFramePhiMin", frame_z, MatCarbonfiber);

  carbon_frame.placeVolume(frame_phi_vol);
  Transform3D trans(RotationZYX((dphi_sec-dphi_seg)*M_PI/180., 0., 0.), Position(0., 0., 0.));
  carbon_frame.placeVolume(frame_phi_vol, trans);

  trans = Transform3D(Translation3D(0., 0., half_dz_sec+Cseg_thick/2.));
  carbon_frame.placeVolume(frame_z_vol, trans); 
  trans = Transform3D(Translation3D(0., 0., -1*(half_dz_sec+Cseg_thick/2.)));
  carbon_frame.placeVolume(frame_z_vol, trans); 

  //carbon_frame.setVisAttributes(description, "CarbonFrameVis");
  //Sector_vol.placeVolume(carbon_frame);

  //    ------- Sensitive material: Grains subtract the holes -------  //
  TGeoTubeSeg* Crystal_module = new TGeoTubeSeg( innerR+Cframe_thick, 
                                          innerR+Cframe_thick+crystal_thick, 
                                          half_dz_sec-Cseg_thick, 
                                          dphi_seg, 
                                          dphi_sec-dphi_seg);
  Tube hole(0., (fiber_r + cladding_thick)+0.002*mm, crystal_thick/2.0);

  Solid crystal_solid(Crystal_module);
  Solid crystal_with_hole = crystal_solid;


  int Nsig_phi = floor((dphi_sec-2*dphi_seg)/dPhi_fiber);
  int Nsig_z = floor(2*(half_dz_sec-Cseg_thick)/dz_fiber); 
  std::cout<<"  Fiber segmentation in crystals: Phi "<<Nsig_phi<<", Z "<<Nsig_z<<std::endl;

  Solid holesLayer; 

  bool firstHole = true;

  //Subtract the holes from crystals
  for(int iphi=0; iphi<Nsig_phi; iphi++) {
    double phi = (dphi_seg + (iphi+0.5)*dPhi_fiber) * M_PI/180.0;


    double rmid = innerR + Cframe_thick + crystal_thick/2.0;
    double x = rmid * cos(phi);
    double y = rmid * sin(phi);
//std::cout<<"Hole position: x "<<x<<", y "<<y<<", Phi "<<phi<<std::endl;

    Rotation3D rot = RotationZ(phi) * RotationY(M_PI/2.0);
    Transform3D tr(rot, Position(x,y,0)); 
//std::cout << "Transform matrix:\n";
//std::cout << tr << std::endl;

    if(firstHole) {
      //holesLayer = Transform3D(tr) * hole;
      holesLayer = UnionSolid("holes_union_first", hole, hole, tr);
      firstHole = false;
    }
    else {
      holesLayer = UnionSolid("holes_union", holesLayer, hole, tr);
    }
  }
  for(int iz=0; iz<Nsig_z; iz++) {
    double zpos = -half_dz_sec + Cseg_thick/2. + (iz+0.5)*dz_fiber;
    Transform3D trz(Position(0,0,zpos));

    crystal_with_hole = SubtractionSolid("crystal_sub", crystal_with_hole, holesLayer, trz);
  }

  Volume crystal_vol("crystal_vol", crystal_with_hole, MatCrystal);
  crystal_vol.setSensitiveDetector(sens);
  crystal_vol.setVisAttributes(description, "SensitiveVis");



  //    ------- WLS fibers and PMMA cladding --------    //
  TGeoTube* Fiber = new TGeoTube(0., fiber_r, (crystal_thick+back_space)/2.0);
  TGeoTube* Cladding = new TGeoTube(fiber_r+0.001*mm, fiber_r+0.001*mm+cladding_thick, (crystal_thick+back_space)/2.0);

  Volume fiberCoreVol("fiberCore", Fiber, MatWLSfiber);
  Volume fiberCladVol("fiberClad", Cladding, MatPMMAcladding);  

  fiberCoreVol.setVisAttributes(description, "WLSFiberVis");
  fiberCladVol.setVisAttributes(description, "FiberCladdingVis");

  // Place fibers
  for(int iphi=0; iphi<Nsig_phi; iphi++) {

    double phi = (dphi_seg + (iphi+0.5)*dPhi_fiber) * M_PI/180.0;

    for(int iz=0; iz<Nsig_z; iz++) {

      double zpos = -half_dz_sec + Cseg_thick/2. + (iz+0.5)*dz_fiber;
      double rmid = innerR + Cframe_thick + (crystal_thick+back_space)/2.0;

      double x = rmid * cos(phi);
      double y = rmid * sin(phi);

      Rotation3D rot = RotationZ(phi) * RotationY(M_PI/2.0);

      Transform3D tr(rot, Position(x,y,zpos));

      Sector_vol.placeVolume(fiberCoreVol, tr);
      Sector_vol.placeVolume(fiberCladVol, tr);
    }
  }

  Sector_vol.placeVolume(crystal_vol);


  // ========  Loop to place sectors in the whole cylender  ========  //
  for(int isec_z=0; isec_z<nModule_z; isec_z++){
    for(int isec_phi=0; isec_phi<nSec_phi; isec_phi++) {

      double phi = isec_phi * dphi_sec * M_PI/180.0;
      double zpos = -halfZ + (isec_z+0.5)*half_dz_sec*2;

      Transform3D tr(RotationZ(phi), Position(0,0,zpos));

      PlacedVolume pv = EcalBarrelVol.placeVolume(Sector_vol, tr);
      pv.addPhysVolID("sector", isec_phi).addPhysVolID("module",isec_z);
    }
  }


  // Finalize geometry

  //Check overlap
  //description.manager().CheckOverlaps(0.05);
  PlacedVolume EcalBarrel_plv = worldVol.placeVolume(EcalBarrelVol);
  EcalBarrel_plv.addPhysVolID("system", x_det.id());
  ECAL.setPlacement(EcalBarrel_plv);

  std::cout << "--> Grainita_ECAL_Barrel_v01::create_detector() end" << std::endl;
  return ECAL;
}

DECLARE_DETELEMENT(Grainita_ECAL_Barrel_v01, create_detector)

//**************************************************************************
