//**************************************************************************
// \file    Grainita_ECAL_Barrel_v1.cpp
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
  std::cout << "--> Grainita_ECAL_Barrel_v1::create_detector() start" << std::endl;

  // Get info from the xml file
  //
  sens.setType("calorimeter");
  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();
  std::cout << "--> Going to create " << det_name << ", with ID: " << x_det.id() << std::endl;



  xml_dim_t x_dim = x_det.dimensions();

  const double CaloX = x_dim.x();
  const double CaloY = x_dim.y();
  const double CaloZ = x_dim.z();
  std::cout << "--> calo dimensions from XML description: x " << CaloX / m
            << " m, y " << CaloY / m << " m, z " << CaloZ / m << " m"
            << std::endl;

  // Retrieve number of layers to populate the calorimeter container with
  //
  auto NumberOfLayers = description.constant<int>("LayersNumber");

  // Info for subdetectors
  //
  xml_det_t x_calo = x_det.child(_Unicode(calo));
  xml_det_t x_calolayer = x_det.child(_Unicode(caloLayer));
  xml_det_t x_abslayer = x_det.child(_Unicode(absLayer));
  xml_det_t x_senslayer = x_det.child(_Unicode(sensLayer));

  auto islayersens = x_senslayer.isSensitive();

  const double CaloLayerX = x_calolayer.x();
  const double CaloLayerY = x_calolayer.y();
  const double CaloLayerZ = x_calolayer.z();

  const double AbsLayerX = x_abslayer.x();
  const double AbsLayerY = x_abslayer.y();
  const double AbsLayerZ = x_abslayer.z();

  const double SensLayerX = x_senslayer.x();
  const double SensLayerY = x_senslayer.y();
  const double SensLayerZ = x_senslayer.z();

  // Create the geometry
  //

  // Create a container for the calorimeter
  //
  Box Calo(CaloX / 2., CaloY / 2., CaloZ / 2.);
  Volume CaloVol("CaloVol", Calo,
                 description.material(x_calo.attr<std::string>(_U(material))));
  CaloVol.setVisAttributes(description, x_calo.visStr());

  

  // Finalize geometry
  //
  DetElement subdet(det_name, x_det.id());
  Volume motherVolume = description.pickMotherVolume(subdet);
  // Place the calo container inside the mother volume
  PlacedVolume CaloPlaced = motherVolume.placeVolume(CaloVol);
  subdet.setPlacement(CaloPlaced);

  std::cout << "--> Grainita_ECAL_Barrel_v1::create_detector() end" << std::endl;
  return subdet;
}

DECLARE_DETELEMENT(Grainita_ECAL_Barrel_v1, create_detector)

//**************************************************************************
