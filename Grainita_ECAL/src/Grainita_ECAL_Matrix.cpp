//**************************************************************************
// \file    Grainita_ECAL_Matrix.cpp
// \brief   Standalone Grainita ECAL matrix module geometry.
//          Port of Grainita_Module/src/construction.cc to DD4hep.
// \author  Fangyi Guo
//**************************************************************************

#include "DDRec/Vector3D.h"
#include <DD4hep/DetFactoryHelper.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;

static Ref_t create_detector(Detector &description, xml_h e,
                             SensitiveDetector sens) {
  std::cout << "--> Grainita_ECAL_Matrix::create_detector() start" << std::endl;

  sens.setType("calorimeter");
  xml_det_t x_det = e;
  const std::string det_name = x_det.nameStr();
  std::cout << "--> Going to create " << det_name << ", with ID: " << x_det.id()
            << std::endl;

  Material air(description.material("Air"));
  Material matCarbonFiber(description.material("CarbonFiber"));
  Material matCrystal(description.material(x_det.attr<std::string>(_U(material))));
  Material matWLSFiber(description.material("WLSFiber"));
  Material matPMMA(description.material("PMMA"));
  Material matFP(description.material("FP"));

  const double module_half_x = description.constant<double>("MatrixModuleHalfX");
  const double module_half_y = description.constant<double>("MatrixModuleHalfY");
  const double crystal_half_z = description.constant<double>("MatrixCrystalHalfZ");
  const int nbox_x = description.constant<int>("MatrixNBoxX");
  const int nbox_y = description.constant<int>("MatrixNBoxY");
  const int nfiber_x = description.constant<int>("MatrixNFiberX");
  const int nfiber_y = description.constant<int>("MatrixNFiberY");
  const int nseg_z = description.constant<int>("MatrixNZSegments");
  const double fiber_radius = description.constant<double>("MatrixFiberRadius");
  const double cladding_l1_thick = description.constant<double>("MatrixCladdingL1Thickness");
  const double cladding_l2_thick = description.constant<double>("MatrixCladdingL2Thickness");
  const double carbonframe_thick = description.constant<double>("MatrixCarbonFrameThickness");
  const double fiber_length = description.constant<double>("MatrixFiberLength");

  const double box_half_x = module_half_x / nbox_x;
  const double box_half_y = module_half_y / nbox_y;
  const double box_half_z = crystal_half_z;
  const double inner_box_half_x = box_half_x - carbonframe_thick;
  const double inner_box_half_y = box_half_y - carbonframe_thick;
  const double cell_half_x = inner_box_half_x / nfiber_x;
  const double cell_half_y = inner_box_half_y / nfiber_y;
  const double cell_half_z = crystal_half_z / nseg_z;
  const double fiber_outer_radius = fiber_radius + cladding_l1_thick + cladding_l2_thick;

  if (inner_box_half_x <= 0. || inner_box_half_y <= 0. || fiber_outer_radius <= 0.) {
    throw std::runtime_error("Grainita_ECAL_Matrix: invalid matrix dimensions; frame/fiber dimensions do not fit the box.");
  }

  std::cout << "  -- Matrix module size: " << 2. * module_half_x << " x "
            << 2. * module_half_y << " x " << 2. * crystal_half_z << std::endl;
  std::cout << "  -- Boxes: " << nbox_x << " x " << nbox_y
            << ", cells per box: " << nfiber_x << " x " << nfiber_y << " x "
            << nseg_z << std::endl;

  DetElement ecal(det_name, x_det.id());
  Volume motherVol = description.pickMotherVolume(ecal);

  Box matrixSolid(module_half_x, module_half_y,
                  std::max(crystal_half_z, fiber_length / 2. ));
  Volume matrixVol(det_name + "_envelope", matrixSolid, air);
  matrixVol.setVisAttributes(description, "SeeThrough");

  Box boxOuterSolid(box_half_x, box_half_y, box_half_z);
  Volume boxVol("matrix_box_vol", boxOuterSolid, air);
  boxVol.setVisAttributes(description, "SeeThrough");

  Box frameOuterSolid(box_half_x, box_half_y, box_half_z);
  Box frameInnerSolid(inner_box_half_x, inner_box_half_y, box_half_z + 0.1 * mm);
  SubtractionSolid frameSolid("matrix_carbon_frame_solid", frameOuterSolid, frameInnerSolid);
  Volume frameVol("matrix_carbon_frame_vol", frameSolid, matCarbonFiber);
  frameVol.setSensitiveDetector(sens);
  frameVol.setVisAttributes(description, "CarbonFrameVis");
  PlacedVolume framePV = boxVol.placeVolume(frameVol);
  framePV.addPhysVolID("component", 4);

  Box cellSolid(cell_half_x, cell_half_y, cell_half_z);
  Tube fiberHoleSolid(0., fiber_outer_radius, cell_half_z + 0.1 * mm);
  SubtractionSolid crystalCellSolid("matrix_cell_with_hole_solid", cellSolid, fiberHoleSolid);
  Volume crystalCellVol("matrix_cell_with_hole_vol", crystalCellSolid, matCrystal);
  crystalCellVol.setSensitiveDetector(sens);
  crystalCellVol.setVisAttributes(description, "SensitiveVis");

  Tube fiberCoreSolid(0., fiber_radius, fiber_length / 2.);
  Tube fiberCladL1Solid(fiber_radius, fiber_radius + cladding_l1_thick, fiber_length / 2.);
  Tube fiberCladL2Solid(fiber_radius + cladding_l1_thick, fiber_outer_radius,
                        fiber_length / 2.);

  Volume fiberCoreVol("matrix_fiber_core_vol", fiberCoreSolid, matWLSFiber);
  Volume fiberCladL1Vol("matrix_fiber_clad_l1_vol", fiberCladL1Solid, matPMMA);
  Volume fiberCladL2Vol("matrix_fiber_clad_l2_vol", fiberCladL2Solid, matFP);
  fiberCoreVol.setSensitiveDetector(sens);
  fiberCladL1Vol.setSensitiveDetector(sens);
  fiberCladL2Vol.setSensitiveDetector(sens);
  fiberCoreVol.setVisAttributes(description, "WLSFiberVis");
  fiberCladL1Vol.setVisAttributes(description, "FiberCladdingVis");
  fiberCladL2Vol.setVisAttributes(description, "FiberCladdingVis");

  std::vector<Position> fiberPositions;
  fiberPositions.reserve(nfiber_x * nfiber_y);

  const double x0_cell = -0.5 * (nfiber_x - 1) * 2. * cell_half_x;
  const double y0_cell = -0.5 * (nfiber_y - 1) * 2. * cell_half_y;
  const double z0_cell = -0.5 * (nseg_z - 1) * 2. * cell_half_z;

  for (int iz = 0; iz < nseg_z; ++iz) {
    const double z = z0_cell + iz * 2. * cell_half_z;
    for (int iy = 0; iy < nfiber_y; ++iy) {
      const double y = y0_cell + iy * 2. * cell_half_y;
      for (int ix = 0; ix < nfiber_x; ++ix) {
        const double x = x0_cell + ix * 2. * cell_half_x;
        PlacedVolume cellPV = boxVol.placeVolume(crystalCellVol, Position(x, y, z));
        cellPV.addPhysVolID("layer", iz);
        cellPV.addPhysVolID("row", iy);
        cellPV.addPhysVolID("column", ix);
        if (iz == 0) {
          fiberPositions.emplace_back(x, y, 0.);
        }
      }
    }
  }

  for (size_t i = 0; i < fiberPositions.size(); ++i) {
    const int iy = static_cast<int>(i) / nfiber_x;
    const int ix = static_cast<int>(i) % nfiber_x;
    PlacedVolume clad2PV = boxVol.placeVolume(fiberCladL2Vol, fiberPositions[i]);
    clad2PV.addPhysVolID("component", 3);
    clad2PV.addPhysVolID("row", iy);
    clad2PV.addPhysVolID("column", ix);

    PlacedVolume clad1PV = boxVol.placeVolume(fiberCladL1Vol, fiberPositions[i]);
    clad1PV.addPhysVolID("component", 2);
    clad1PV.addPhysVolID("row", iy);
    clad1PV.addPhysVolID("column", ix);

    PlacedVolume corePV = boxVol.placeVolume(fiberCoreVol, fiberPositions[i]);
    corePV.addPhysVolID("component", 1);
    corePV.addPhysVolID("row", iy);
    corePV.addPhysVolID("column", ix);
  }

  const double x0_box = -0.5 * (nbox_x - 1) * 2. * box_half_x;
  const double y0_box = -0.5 * (nbox_y - 1) * 2. * box_half_y;
  int boxCopyNo = 0;
  for (int iy = 0; iy < nbox_y; ++iy) {
    const double y = y0_box + iy * 2. * box_half_y;
    for (int ix = 0; ix < nbox_x; ++ix) {
      const double x = x0_box + ix * 2. * box_half_x;
      PlacedVolume boxPV = matrixVol.placeVolume(boxVol, Position(x, y, 0.));
      boxPV.addPhysVolID("box", ++boxCopyNo);
      boxPV.addPhysVolID("box_x", ix);
      boxPV.addPhysVolID("box_y", iy);
    }
  }

  PlacedVolume matrixPV = motherVol.placeVolume(matrixVol);
  matrixPV.addPhysVolID("system", x_det.id());
  ecal.setPlacement(matrixPV);

  std::cout << "--> Grainita_ECAL_Matrix::create_detector() end" << std::endl;
  return ecal;
}

DECLARE_DETELEMENT(Grainita_ECAL_Matrix, create_detector)

//**************************************************************************
