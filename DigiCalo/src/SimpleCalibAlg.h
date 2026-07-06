#ifndef SIMPLE_CALIB_ALG_H
#define SIMPLE_CALIB_ALG_H

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/SmartIF.h"
#include "k4FWCore/DataHandle.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Segmentations.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "k4Interface/IGeoSvc.h"

#include <map>
#include <string>
#include <utility>

class SimpleCalibAlg : public Gaudi::Algorithm {
public:
  SimpleCalibAlg(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  double rhoPitchCorrection(uint64_t cellID) const;
  double meanWindowResponse(double pitchPhi, double pitchTheta) const;
  double responseFunction(double distance) const;

  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
  int m_rhoIndex = -1;
  int m_thetaIndex = -1;

  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_inputHits{
      "GrainitaEcalBarrelDigiHit", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_outputHits{
      "GrainitaEcalBarrelCalibHit", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<std::string> m_readoutName{this, "ReadOutName", "GrainitaEcalBarrelRO"};
  Gaudi::Property<std::string> m_inputFormat{this, "InputFormat", "Energy"};
  Gaudi::Property<double> m_calibrationConstant{this, "CalibrationConstant", 1.0};
  Gaudi::Property<double> m_lightYield{this, "LightYield", 10.0};

  // Light response model, parameters unit in cm. 
  Gaudi::Property<bool> m_applyRhoPitchCorrection{this, "ApplyRhoPitchCorrection", true};
  Gaudi::Property<double> m_attLength{this, "AttLength", 0.34};
  Gaudi::Property<double> m_responseX0{this, "ResponseX0", 0.0856};
  Gaudi::Property<double> m_responseSlope{this, "ResponseSlope", 9.276};
  Gaudi::Property<double> m_responseIntercept{this, "ResponseIntercept", 0.206};
  Gaudi::Property<int> m_neighborCellSize{this, "NeighborCellSize", 5};
  Gaudi::Property<int> m_nScan{this, "NScan", 80};
  Gaudi::Property<int> m_referenceRhoID{this, "ReferenceRhoID", 0};

  mutable std::map<std::pair<int, int>, double> m_correctionCache;
};

#endif
