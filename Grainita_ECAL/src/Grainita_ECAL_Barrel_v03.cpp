//**************************************************************************
// \file    Grainita_ECAL_Barrel_v03.cpp
// \brief:  Grainita ECAL full detector geometry, with pointing structure
//          and individual volume for cells and fibers. 
// \author: Fangyi Guo (fangyi.guo@cern.ch) start from DD4HEP tutorial in 
//          DRD Calo (https://github.com/DRDCalo/DD4hepTutorials)
// \date:   April 2026
//**************************************************************************

// Includers from DD4hep
#include "detectorSegmentations/DRparamBarrel_k4geo.h"

#include "DDRec/Vector3D.h"
#include <DD4hep/DetFactoryHelper.h>

#include <algorithm>
#include <cmath>

using namespace dd4hep;

// Build geometry
//
static Ref_t create_detector(Detector &description, xml_h e,
                             SensitiveDetector sens) {
  std::cout << "--> Grainita_ECAL_Barrel_v03::create_detector() start" << std::endl;

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
  auto nSec_phi = description.constant<int>("NumberOfPhiSectors");
  auto nModule_z = description.constant<int>("NumberOfZModules");

  auto Cframe_thick = description.constant<double>("FrameThickness");
  auto Cseg_thick = description.constant<double>("SegThickness");
  auto module_tilt = description.constant<double>("ModuleTiltAngle"); //Unit in degree
  double module_tilt_rad = module_tilt*M_PI/180.;

  auto fiber_r = description.constant<double>("FiberRadius");
  auto cladding_thick = description.constant<double>("CladdingThickness");
  auto fiber_pitch = description.constant<double>("FiberPitchSize"); 


  double dphi_sec = 2*M_PI/nSec_phi;
  double outerR = (innerR + crystal_thick + back_space + 2*Cframe_thick)/cos(dphi_sec/2.); 
  std::cout<<"  -- Dimension: inner R "<<innerR<<", outer R "<<outerR<<", inner half Z "<<halfZ<<std::endl;
  std::cout<<"     Z segmentation: "<<nModule_z<<", Phi segmentation: "<<nSec_phi<<std::endl;

  //Calculate the shift from phi direction tilt. 
  // Note: The tilt works on the outer vertex.
  // The section of sector in phi-plane is defined with 4 points (in R-theta coordinate):
  //    Inner: A(Rin, -dphi_sec/2.),  D(Rin, dphi_sec/2.)
  //    Original outer withouth tilt: (Rout, -dphi_sec/2.), (Rout, dphi_sec/2.)
  //    Outer with tilt: B(Rout, -dphi_sec/2. + tilt), C(Rout, dphi_sec/2. + tilt)
  // BUT this gives non-parallel inner and outer surface. We force the outer surface be parallel with inner surface
  // in the plane of x=Rout*cos(dphi_sec/2.). 
  //    (in x-y coordinate)
  //    B' (Rout*cos(dphi_sec/2.), y_B'), C'(Rout*cos(dphi_sec/2.), y_C')
  // where B' and C' are from extending AB and DC to x=Rout*cos(dphi_sec/2.) plane. 
  // To get y_B' and y_C', we calculated the function of AB and DC. 
  //
  //        \
  //     \   \ /|C'
  //      \  // |
  //       \//  |
  //      D |   |
  //        |   |
  //       A --- B'

  //Slope and intercept of line AB
  double tmp_slope1 = -(outerR*sin(dphi_sec/2.-module_tilt_rad) - innerR*tan(dphi_sec/2.)) / 
                      (outerR*cos(dphi_sec/2.-module_tilt_rad) - innerR);
  double tmp_intercept1 = -innerR*tan(dphi_sec/2) - innerR*tmp_slope1;

  //Slope and intercept of line DC
  double tmp_slope2 = (outerR*sin(dphi_sec/2+module_tilt_rad) - innerR*tan(dphi_sec/2)) /
                      (outerR*cos(dphi_sec/2+module_tilt_rad) - innerR);
  double tmp_intercept2 = innerR*tan(dphi_sec/2) - innerR*tmp_slope2;

  //Calculate the shift regarding to the non-tilt point (Rout, +-dphi_sec/2.).  
  double shift_phi_neg = fabs( fabs(tmp_slope1*outerR*cos(dphi_sec/2.)+tmp_intercept1 ) - outerR*sin(dphi_sec/2.) );
  double shift_phi_pos = fabs( fabs(tmp_slope2*outerR*cos(dphi_sec/2.)+tmp_intercept2 ) - outerR*sin(dphi_sec/2.) );

  std::cout<<"-- Consider a tilt angle in phi to avoid pointing cracks. Tilt angle (in deg): "<<module_tilt<<", in rad: "<<module_tilt_rad<<std::endl;
  std::cout<<"Input vars: Rin "<<innerR<<", Rout "<<outerR<<", phi "<<dphi_sec/2.<<", delta "<<module_tilt_rad<<std::endl;
  std::cout<<"  Calculate the new boundary lines: "<<std::endl;
  std::cout<<"    in -y side: k = "<<tmp_slope1<<", b = "<<tmp_intercept1<<std::endl;
  std::cout<<"    in +y side: k = "<<tmp_slope2<<", b = "<<tmp_intercept2<<std::endl;
  std::cout<<"  Calculated +y direction shift: "<<shift_phi_pos<<std::endl;
  std::cout<<"  Calculated -y direction shift: "<<shift_phi_neg<<std::endl;
  double halfZ_out = outerR*cos(dphi_sec/2.)/tan(atan(innerR / halfZ) - module_tilt_rad);
  double tilt_rad = atan( (outerR*sin(dphi_sec/2.)+shift_phi_pos)/(outerR*cos(dphi_sec/2.)) );
  std::cout<<"  halfZ out = "<<outerR*cos(dphi_sec/2.)<<" / tan("<<atan(innerR / halfZ) - module_tilt_rad<<") "<<std::endl;
  outerR = sqrt( pow(outerR*cos(dphi_sec/2.), 2) + pow((outerR*sin(dphi_sec/2.) + shift_phi_pos) , 2) );
  std::cout<<"  Outer radius and halfZ considering this: "<<outerR<<", "<<halfZ_out<<std::endl;

  // Create the geometry
  DetElement ECAL(det_name, x_det.id());
  Volume worldVol = description.pickMotherVolume(ECAL);

  TGeoTube* EcalBarrel = new TGeoTube(innerR, outerR, halfZ_out);
  Volume EcalBarrelVol("EcalBarrel", EcalBarrel, air);
  //EcalBarrelVol.setVisAttributes(description, "CaloVis");

  const double cladding_outer_r = fiber_r + cladding_thick;
  const double fiber_geom_length = std::max(0., crystal_thick - 105.0 * mm);
  TGeoTube* FiberCore = new TGeoTube(0., fiber_r, fiber_geom_length / 2.0);
  TGeoTube* FiberCladding = new TGeoTube(fiber_r, cladding_outer_r, fiber_geom_length / 2.0);
  Volume fiberCoreVol("fiberCore", FiberCore, MatWLSfiber);
  Volume fiberCladVol("fiberClad", FiberCladding, MatPMMAcladding);
  fiberCoreVol.setVisAttributes(description, "WLSFiberVis");
  fiberCladVol.setVisAttributes(description, "FiberCladdingVis");

  // ====== Carbon fiber frame as supporting, at inner and outer of barrel. ====== //
  PolyhedraRegular innerFrame(nSec_phi, -dphi_sec/2., innerR, innerR + Cframe_thick, 2*halfZ);
  Volume innerFrame_vol("InnerCarbonFrame", innerFrame, MatCarbonfiber);

  PolyhedraRegular outerFrame(nSec_phi, tilt_rad, 
                              outerR*cos(dphi_sec/2.)-Cframe_thick, 
                              outerR*cos(dphi_sec/2.), 
                              2*halfZ_out);
  Volume outerFrame_vol("outerCarbonFrame", outerFrame, MatCarbonfiber);

  innerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  outerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  EcalBarrelVol.placeVolume(innerFrame_vol);
  //EcalBarrelVol.placeVolume(outerFrame_vol);

  
  // ======== Define a sector (backspace is included in sector)======== //
  double innerR_sector = innerR + Cframe_thick;
  double outerR_sector = (innerR_sector + crystal_thick + back_space)/cos(dphi_sec/2.);;
  halfZ_out = outerR_sector*cos(dphi_sec/2.)/tan(atan(innerR_sector / halfZ) - module_tilt_rad);
std::cout<<"  halfZ out = "<<outerR_sector*cos(dphi_sec/2.)<<" / tan("<<atan(innerR_sector / halfZ) - module_tilt_rad<<") = "<<halfZ_out<<std::endl;

  // Sector inner and outer width
  double inner_width_neg = fabs(tmp_slope1*innerR_sector + tmp_intercept1);
  double inner_width_pos = fabs(tmp_slope2*innerR_sector + tmp_intercept2);
  double outer_width_neg = fabs(tmp_slope1*outerR_sector*cos(dphi_sec/2.) + tmp_intercept1);
  double outer_width_pos = fabs(tmp_slope2*outerR_sector*cos(dphi_sec/2.) + tmp_intercept2);

  std::cout<<std::endl;
  std::cout<<"  Calculate the new boundary lines in Sector: "<<std::endl;
  std::cout<<"    in -y side: k = "<<tmp_slope1<<", b = "<<tmp_intercept1<<std::endl;
  std::cout<<"    in +y side: k = "<<tmp_slope2<<", b = "<<tmp_intercept2<<std::endl;
  std::cout<<"  Calculated +y direction shift: "<<shift_phi_pos<<std::endl;
  std::cout<<"  Calculated -y direction shift: "<<shift_phi_neg<<std::endl;


  double vertices_sector[16] = {
    halfZ,  inner_width_pos,
    halfZ,  -inner_width_neg,
    -halfZ, -inner_width_neg,
    -halfZ, inner_width_pos,

    halfZ_out,  outer_width_pos,
    halfZ_out,  -outer_width_neg,
    -halfZ_out, -outer_width_neg,
    -halfZ_out, outer_width_pos,
  };

  EightPointSolid Sector("Sector", (crystal_thick+back_space)/2., vertices_sector );
  Volume Sector_vol("sector_vol", Sector, air);
  Sector_vol.setVisAttributes(description, "CaloVis");


  /// ------ Define module along Z axis ------  ///
  double theta_min = atan(innerR_sector/halfZ);
  double theta_max = M_PI - theta_min; 
  double theta_module = (theta_max-theta_min)/nModule_z; 
  for(int iz=0; iz<nModule_z; iz++){
    if(iz != nModule_z / 2 - 1 ) continue; 

    double theta_min_module = theta_min + iz*theta_module; 
    double theta_max_module = theta_min_module + theta_module; 
    double tilt_min = theta_min_module<M_PI/2. ? -module_tilt*M_PI/180. : module_tilt*M_PI/180.; 
    double tilt_max = theta_max_module<M_PI/2. ? -module_tilt*M_PI/180. : module_tilt*M_PI/180.; 
    double tilt_min_local = atan( (outerR_sector*cos(dphi_sec/2.)-innerR_sector)/ 
                                    (outerR_sector*cos(dphi_sec/2.)/tan(theta_min_module+tilt_min) - 
                                    innerR_sector/tan(theta_min_module)) );
    double tilt_max_local = atan( (outerR_sector*cos(dphi_sec/2.)-innerR_sector)/
                                    (outerR_sector*cos(dphi_sec/2.)/tan(theta_max_module+tilt_max) -
                                    innerR_sector/tan(theta_max_module)) );

    std::cout<<"Module #"<<iz<<" Geometry parameters: "<<std::endl;
    std::cout<<"Sector outerR - innerR = " << (outerR_sector*cos(dphi_sec/2.)-innerR_sector);
    std::cout<<", width (min) = "<<fabs(outerR_sector*cos(dphi_sec/2.)/tan(theta_min_module+tilt_min) -
                                    innerR_sector/tan(theta_min_module));
    std::cout<<", width (max) = "<<fabs(outerR_sector*cos(dphi_sec/2.)/tan(theta_max_module+tilt_max) -
                                    innerR_sector/tan(theta_max_module))<<std::endl;

    // Module height: keeping the total length = crystal thickness, so h = thick*sin(theta)
    // Exception in the most central module: h = thick. 
    double height_module = (theta_max_module<M_PI/2.) ? crystal_thick*fabs(sin(tilt_max_local)) : crystal_thick*fabs(sin(tilt_min_local));
    if( (theta_max_module-M_PI/2.)*(theta_min_module-M_PI/2.)<0 ) height_module = crystal_thick;

    //Module width: 
    double outer_width_module_neg = fabs(tmp_slope1*(innerR_sector+height_module) + tmp_intercept1);
    double outer_width_module_pos = fabs(tmp_slope2*(innerR_sector+height_module) + tmp_intercept2);


    std::cout<<"  Theta range: "<<theta_min_module<<", "<<theta_max_module<<", tilt angle "<<tilt_min<<", "<<tilt_max<<", local tilt angle "<<tilt_min_local<<", "<<tilt_max_local<<std::endl;
    std::cout<<"  height: "<<height_module<<", inner width "<<inner_width_neg<<", "<<inner_width_pos;
    std::cout<<", outer width "<<outer_width_module_neg<<", "<<outer_width_module_pos<<std::endl;
    
    double vertices_frame[16] = {
      innerR_sector/tan(theta_min_module), inner_width_pos,
      innerR_sector/tan(theta_min_module), -inner_width_neg,
      innerR_sector/tan(theta_max_module), -inner_width_neg,
      innerR_sector/tan(theta_max_module), inner_width_pos,

      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local), outer_width_module_pos,
      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local), -outer_width_module_neg,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local), -outer_width_module_neg,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local), outer_width_module_pos,
    };

    std::cout<<"  Inner surface vertices: "<<std::endl;
    for(int i=0; i<8; i+=2) std::cout<<"    ["<<vertices_frame[i]<<", "<<vertices_frame[i+1]<<"], "<<std::endl;
    std::cout<<"  Outer surface vertices: "<<std::endl;
    for(int i=8; i<16; i+=2) std::cout<<"    ["<<vertices_frame[i]<<", "<<vertices_frame[i+1]<<"], "<<std::endl;


    double vertices_crystal[16] = {
      innerR_sector/tan(theta_min_module)-Cseg_thick, inner_width_pos-Cseg_thick,
      innerR_sector/tan(theta_min_module)-Cseg_thick, -inner_width_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+Cseg_thick, -inner_width_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+Cseg_thick, inner_width_pos-Cseg_thick,

      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local)-Cseg_thick, outer_width_module_pos-Cseg_thick,
      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local)-Cseg_thick, -outer_width_module_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local)+Cseg_thick, -outer_width_module_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local)+Cseg_thick, outer_width_module_pos-Cseg_thick,
    };

    EightPointSolid ModuleFrame(Form("Module_%d", iz), height_module/2., vertices_frame );
    Volume ModuleFrame_vol(Form("module_vol%d", iz), ModuleFrame, MatCarbonfiber);
    ModuleFrame_vol.setVisAttributes(description, "CarbonFrameVis");

    EightPointSolid ModuleCrystal(Form("ModuleCrystal_%d", iz), height_module/2., vertices_crystal );
    Volume ModuleCrystal_vol(Form("moduleCrystal_vol%d", iz), ModuleCrystal, air);
    ModuleCrystal_vol.setVisAttributes(description, "SeeThroughVis");
/*
    // Define cell geometry: use a trapezoid like done in DRCalo 
    double theta_max_module_tilt = atan(height_module/(vertices_crystal[12]-vertices_crystal[4]));
    double theta_min_module_tilt = atan(height_module/(vertices_crystal[8]-vertices_crystal[0]));
    double phi_pos_module_tilt = atan( (vertices_crystal[9]-vertices_crystal[1]) / height_module );
    double phi_neg_module_tilt = atan( (vertices_crystal[11]-vertices_crystal[3]) / height_module );

    int n_cell_x = fabs(vertices_crystal[0]-vertices_crystal[4]) / (fiber_pitch*sin(theta_max_module_tilt));
    int n_cell_y = fabs(vertices_crystal[1]-vertices_crystal[3]) / (fiber_pitch*sin(phi_neg_module_tilt));
    double cell_depth = height_module;
    double cell_x1 = fabs(vertices_crystal[0]-vertices_crystal[4])/n_cell_x;  
    double cell_x2 = fabs(vertices_crystal[8]-vertices_crystal[12])/n_cell_x;
    double cell_y1 = fabs(vertices_crystal[1]-vertices_crystal[3])/n_cell_y;
    double cell_y2 = fabs(vertices_crystal[9]-vertices_crystal[11])/n_cell_y;

    Trap CellCrystal(Form("CrystalCell_%d_template", iz), cell_depth / 2., template_trap_theta, template_trap_phi,
                     template_inner_dy, template_inner_dx, template_inner_dx, 0.,
                     template_outer_dy, template_outer_dx, template_outer_dx, 0.);
    Volume CellCrystal_vol(Form("crystalCell_vol%d_template", iz), CellCrystal, MatCrystal);
    CellCrystal_vol.setVisAttributes(description, "SensitiveVis");
    CellCrystal_vol.setSensitiveDetector(sens);

    // Place fiber inside the cell. 

    // Loop to place cells in module, with rotation to make them pointing to IP.
    for (int ix_cell = 0; ix_cell < n_cell_x; ++ix_cell) {
      const double x_fraction_center = (static_cast<double>(ix_cell) + 0.5) / static_cast<double>(n_cell_x);
      const double inner_x_center = inner_x_min_edge + x_fraction_center * inner_x_span;
      const double outer_x_center = outer_x_min_edge + x_fraction_center * (outer_x_max_edge - outer_x_min_edge);

      for (int iy_cell = 0; iy_cell < n_cell_y; ++iy_cell) {
        const double y_fraction_center = (static_cast<double>(iy_cell) + 0.5) / static_cast<double>(n_cell_y);
        const double inner_y_center = inner_y_min_edge + y_fraction_center * inner_y_span;
        const double outer_y_center = outer_y_min_edge + y_fraction_center * (outer_y_max_edge - outer_y_min_edge);
        const double target_axis_x = outer_x_center - inner_x_center;
        const double target_axis_y = outer_y_center - inner_y_center;
        const double cell_rot_z = std::atan2(target_axis_y, target_axis_x) - template_trap_phi;
        const double cos_rot = std::cos(cell_rot_z);
        const double sin_rot = std::sin(cell_rot_z);
        const double rotated_axis_x = cos_rot * template_axis_dx - sin_rot * template_axis_dy;
        const double rotated_axis_y = sin_rot * template_axis_dx + cos_rot * template_axis_dy;
        const double cell_origin_x = inner_x_center + 0.5 * rotated_axis_x;
        const double cell_origin_y = inner_y_center + 0.5 * rotated_axis_y;

        Transform3D cell_tr(RotationZ(cell_rot_z), Translation3D(cell_origin_x, cell_origin_y, 0.));
        PlacedVolume cell_pv = ModuleCrystal_vol.placeVolume(CellCrystal_vol, cell_tr);
        cell_pv.addPhysVolID("theta", ix_cell);
        cell_pv.addPhysVolID("phi", iy_cell);
        cell_pv.addPhysVolID("rho", 0);
      }
    }
*/
    PlacedVolume pv = ModuleFrame_vol.placeVolume(ModuleCrystal_vol);
    pv.addPhysVolID("stave", iz);

    Transform3D tr(Translation3D(0, 0, -(crystal_thick+back_space-height_module)/2.));
    Sector_vol.placeVolume(ModuleFrame_vol, tr);
  }



  //========  Loop to place sectors  ========  //
  //for(int isec_phi=0; isec_phi<nSec_phi; isec_phi++){
  //  double phi = isec_phi * dphi_sec;
  //  Transform3D tr(RotationY(M_PI/2.0), Translation3D(innerR_sector+(crystal_thick+back_space)/2., 0., 0.));
  //  PlacedVolume pv = EcalBarrelVol.placeVolume(Sector_vol, RotationZ(phi) * tr);
  //  pv.addPhysVolID("sector", isec_phi);
  //  break;
  //}
  */

  // DRparamBarrel_k4geo tower-array demo. This mirrors the FiberDualReadoutCalo
  // pattern: set parameters, init(), build Trap from computed dimensions, and
  // place each tower with the computed transform.
  dd4hep::DDSegmentation::DRparamBarrel_k4geo drParam;
  const int demoThetaTowers = std::max(1, nModule_z);
  const int demoPhiTowers = std::max(1, nSec_phi);
  const double thetaCoverage = 2. * std::atan(halfZ / innerR_sector);
  const double demoDeltaTheta = thetaCoverage / static_cast<double>(demoThetaTowers);
  const double demoThetaStart = -0.5 * thetaCoverage;

  drParam.SetIsRHS(true);
  drParam.SetInnerX(innerR_sector);
  drParam.SetTowerH(crystal_thick);
  drParam.SetNumZRot(demoPhiTowers);
  drParam.SetSipmHeight(back_space);

  Assembly towerArrayAssembly("DRparamTowerArrayAssembly");

  for (int iTheta = 0; iTheta < demoThetaTowers; ++iTheta) {
    const double thetaCenter = demoThetaStart + (static_cast<double>(iTheta) + 0.5) * demoDeltaTheta;

    drParam.SetDeltaTheta(demoDeltaTheta);
    drParam.SetThetaOfCenter(thetaCenter);
    drParam.SetCurrentTowerNum(iTheta);
    drParam.init();

    Trap towerTrap("DRparamTower_" + std::to_string(iTheta),
                   crystal_thick / 2., 0., 0.,
                   drParam.GetH1(), drParam.GetBl1(), drParam.GetTl1(), 0.,
                   drParam.GetH2(), drParam.GetBl2(), drParam.GetTl2(), 0.);
    Volume towerVol("DRparamTowerVol_" + std::to_string(iTheta), towerTrap, MatCrystal);
    towerVol.setVisAttributes(description, "SensitiveVis");
    towerVol.setSensitiveDetector(sens);

    //Trap sipmTrap("DRparamTowerSipm_" + std::to_string(iTheta),
    //              back_space / 2., 0., 0.,
    //              drParam.GetH2(), drParam.GetBl2(), drParam.GetTl2(), 0.,
    //              drParam.GetH2sipm(), drParam.GetBl2sipm(), drParam.GetTl2sipm(), 0.);
    //Volume sipmVol("DRparamTowerSipmVol_" + std::to_string(iTheta), sipmTrap, air);
    //sipmVol.setVisAttributes(description, "CaloLayerVis");

    for (int iPhi = 0; iPhi < demoPhiTowers; ++iPhi) {
      PlacedVolume towerPV = towerArrayAssembly.placeVolume(towerVol, drParam.GetTransform3D(iPhi));
      towerPV.addPhysVolID("stave", iTheta);
      towerPV.addPhysVolID("sector", iPhi);

      //PlacedVolume sipmPV = towerArrayAssembly.placeVolume(sipmVol, drParam.GetSipmTransform3D(iPhi));
      //sipmPV.addPhysVolID("stave", iTheta);
      //sipmPV.addPhysVolID("sector", iPhi);
    }
  }

  drParam.filled();
  drParam.SetTotTowerNum(demoThetaTowers);
  drParam.finalized();

  std::cout << "--> DRparam tower-array demo:"
            << " thetaTowers=" << demoThetaTowers
            << ", phiTowers=" << demoPhiTowers
            << ", deltaTheta=" << demoDeltaTheta
            << std::endl;

  // Finalize geometry
  PlacedVolume EcalBarrel_plv = worldVol.placeVolume(towerArrayAssembly);
  EcalBarrel_plv.addPhysVolID("system", x_det.id());
  ECAL.setPlacement(EcalBarrel_plv);

  std::cout << "--> Grainita_ECAL_Barrel_v03::create_detector() end" << std::endl;
  return ECAL;
}

DECLARE_DETELEMENT(Grainita_ECAL_Barrel_v03, create_detector)

//**************************************************************************
