#include "SimpleCalibAlg.h"

#include <algorithm>
#include <cmath>
#include <exception>

DECLARE_COMPONENT(SimpleCalibAlg)

SimpleCalibAlg::SimpleCalibAlg(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("InputCaloHitCollection", m_inputHits, "Input CalorimeterHit collection");
  declareProperty("OutputCaloHitCollection", m_outputHits, "Output calibrated CalorimeterHit collection");
}

StatusCode SimpleCalibAlg::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) {
    return sc;
  }

  const auto& format = m_inputFormat.value();
  if (format != "Energy" && format != "Npe" && format != "ADC") {
    error() << "InputFormat must be Energy, Npe, or ADC, got " << format << endmsg;
    return StatusCode::FAILURE;
  }
  if (format == "Npe" && m_lightYield.value() <= 0.) {
    error() << "LightYield must be positive for InputFormat=Npe" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_applyRhoPitchCorrection.value()) {
    if (m_attLength.value() <= 0. || m_neighborCellSize.value() <= 0 || m_nScan.value() <= 0) {
      error() << "AttLength, NeighborCellSize, and NScan must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_responseX0.value() <= 0.) {
      error() << "ResponseX0 must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_referenceRhoID.value() < 0) {
      error() << "ReferenceRhoID must be non-negative" << endmsg;
      return StatusCode::FAILURE;
    }

    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate GeoSvc" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_geoSvc->getDetector()->readouts().find(m_readoutName) ==
        m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout " << m_readoutName << " does not exist" << endmsg;
      return StatusCode::FAILURE;
    }

    m_segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation();
    if (!m_segmentation) {
      error() << "Readout " << m_readoutName << " does not have a segmentation" << endmsg;
      return StatusCode::FAILURE;
    }

    m_decoder = m_segmentation->decoder();
    if (!m_decoder) {
      error() << "Readout " << m_readoutName << " does not have a decoder" << endmsg;
      return StatusCode::FAILURE;
    }

    try {
      m_rhoIndex = m_decoder->index("rho");
      m_thetaIndex = m_decoder->index("theta");
    } catch (const std::exception& e) {
      error() << "Readout " << m_readoutName << " must contain rho and theta fields: " << e.what()
              << endmsg;
      return StatusCode::FAILURE;
    }
  }

  info() << "SimpleCalibAlg initialized with InputFormat=" << format
         << ", CalibrationConstant=" << m_calibrationConstant.value()
         << ", ApplyRhoPitchCorrection=" << m_applyRhoPitchCorrection.value() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode SimpleCalibAlg::execute(const EventContext&) const {
  const auto* inputHits = m_inputHits.get();
  auto* outputHits = m_outputHits.createAndPut();

  for (const auto& hit : *inputHits) {
    double calibratedE = hit.getEnergy() * m_calibrationConstant.value();
    if (m_inputFormat.value() == "Npe") {
      calibratedE /= m_lightYield.value();
    }

    double rhoCorrection = rhoPitchCorrection(hit.getCellID());
    // std::cout<<"Hit energy: "<<calibratedE<<", rhoCorrection: "<<rhoCorrection<<std::endl;
    calibratedE *= rhoCorrection;

    auto outHit = outputHits->create();
    outHit.setCellID(hit.getCellID());
    outHit.setEnergy(calibratedE);
    outHit.setEnergyError(hit.getEnergyError());
    outHit.setTime(hit.getTime());
    outHit.setPosition(hit.getPosition());
    outHit.setType(hit.getType());
  }

  return StatusCode::SUCCESS;
}

StatusCode SimpleCalibAlg::finalize() {
  info() << "Cached " << m_correctionCache.size() << " rho pitch correction factors" << endmsg;
  return Gaudi::Algorithm::finalize();
}

double SimpleCalibAlg::rhoPitchCorrection(uint64_t cellID) const {
  if (!m_applyRhoPitchCorrection.value()) {
    return 1.;
  }

  const int rhoID = static_cast<int>(m_decoder->get(cellID, m_rhoIndex));
  const int thetaID = static_cast<int>(m_decoder->get(cellID, m_thetaIndex));
  const auto key = std::make_pair(rhoID, thetaID);
  const auto cached = m_correctionCache.find(key);
  if (cached != m_correctionCache.end()) {
    return cached->second;
  }

  uint64_t refCellID = cellID;
  m_decoder->set(refCellID, m_rhoIndex, m_referenceRhoID.value());

  const auto position = m_segmentation->position(cellID);
  const auto refPosition = m_segmentation->position(refCellID);
  const auto dimensions = m_segmentation->cellDimensions(cellID);
  // std::cout<<" Hit position: ("<<position.X<<", "<<position.Y<<", "<<position.Z<<") "<<std::endl;

  if (dimensions.size() < 2) {
    warning() << "Segmentation cellDimensions has size " << dimensions.size()
              << "; disabling rho pitch correction for this hit" << endmsg;
    m_correctionCache[key] = 1.;
    return 1.;
  }

  const double transverse = std::hypot(position.X, position.Y);
  const double rho = std::sqrt(position.X * position.X + position.Y * position.Y + position.Z * position.Z);
  const double refTransverse = std::hypot(refPosition.X, refPosition.Y);
  const double refRho =
  std::sqrt(refPosition.X * refPosition.X + refPosition.Y * refPosition.Y + refPosition.Z * refPosition.Z);

  const double pitchPhi = transverse * dimensions[0];
  const double pitchTheta = rho * dimensions[1];
  const double refPitchPhi = refTransverse * dimensions[0];
  const double refPitchTheta = refRho * dimensions[1];

  // std::cout<<"pitchPhi: "<<pitchPhi<<", pitchTheta: "<<pitchTheta<<", refPitchPhi: "<<refPitchPhi<<", refPitchTheta: "<<refPitchTheta<<std::endl;

  const double response = meanWindowResponse(pitchPhi, pitchTheta);
  const double refResponse = meanWindowResponse(refPitchPhi, refPitchTheta);
  const double correction = response > 0. ? refResponse / response : 1.;
  // std::cout<<"response: "<<response<<", refResponse: "<<refResponse<<", correction: "<<correction<<std::endl;

  m_correctionCache[key] = correction;
  return correction;
}

double SimpleCalibAlg::meanWindowResponse(double pitchPhi, double pitchTheta) const {
  const int neighborRadius = std::max(1, m_neighborCellSize.value()) / 2;
  double sum = 0.;
  long n = 0;

  for (int iu = 0; iu < m_nScan.value(); ++iu) {
    const double u = (-0.5 + (iu + 0.5) / m_nScan.value()) * pitchPhi;
    for (int iv = 0; iv < m_nScan.value(); ++iv) {
      const double v = (-0.5 + (iv + 0.5) / m_nScan.value()) * pitchTheta;

      double responseSum = 0.;
      for (int iphi = -neighborRadius; iphi <= neighborRadius; ++iphi) {
        for (int itheta = -neighborRadius; itheta <= neighborRadius; ++itheta) {
          const double du = u - iphi * pitchPhi;
          const double dv = v - itheta * pitchTheta;
          responseSum += responseFunction(std::sqrt(du * du + dv * dv));
        }
      }
      sum += responseSum;
      ++n;
    }
  }

  return n > 0 ? sum / n : 1.;
}

double SimpleCalibAlg::responseFunction(double distance) const {
  const double x0 = m_responseX0.value();
  const double slope = m_responseSlope.value();
  const double intercept = m_responseIntercept.value();
  const double norm = 1. / std::exp(-x0 / m_attLength.value());
  if (distance < x0) {
    return slope * distance + intercept;
  }
  return norm * std::exp(-distance / m_attLength.value());
}
