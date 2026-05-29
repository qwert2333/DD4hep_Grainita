#include "detectorSegmentations/FCCSWModularGridRhoPhiTheta_k4geo.h"

#include <cmath>

namespace dd4hep {
namespace DDSegmentation {

  FCCSWModularGridRhoPhiTheta_k4geo::FCCSWModularGridRhoPhiTheta_k4geo(const std::string& cellEncoding)
      : FCCSWGridRhoPhiTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWModularGridRhoPhiTheta_k4geo::FCCSWModularGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder)
      : FCCSWGridRhoPhiTheta_k4geo(decoder) {
    commonSetup();
  }

  void FCCSWModularGridRhoPhiTheta_k4geo::commonSetup() {
    _type = "FCCSWModularGridRhoPhiTheta_k4geo";
    _description = "Rho-Phi-theta segmentation with modular tilted virtual fibers";

    registerParameter("enable_tilted_fiber", "Enable modular tilted virtual fiber positions", m_enableTiltedFiber,
                      true);
    registerParameter("inner_R", "Detector inner radius", m_innerR, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("half_Z", "Detector inner half length", m_halfZ, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("crystal_thickness", "Crystal thickness", m_crystalThickness, 0.,
                      SegmentationParameter::LengthUnit, true);
    registerParameter("back_space", "Back space behind crystal", m_backSpace, 0., SegmentationParameter::LengthUnit,
                      true);
    registerParameter("frame_thickness", "Carbon frame thickness", m_frameThickness, 0.,
                      SegmentationParameter::LengthUnit, true);
    registerParameter("seg_thickness", "Inter-segment frame thickness", m_segThickness, 0.,
                      SegmentationParameter::LengthUnit, true);
    registerParameter("module_tilt_angle", "Module tilt angle in degrees", m_moduleTiltAngle, 0.,
                      SegmentationParameter::NoUnit, true);
    registerParameter("n_phi_sectors", "Number of phi sectors", m_nPhiSectors, 1);
    registerParameter("n_z_modules", "Number of z/theta modules", m_nZModules, 1);
    registerIdentifier("identifier_sector", "Cell ID identifier for phi sector", m_sectorID, "sector");
    registerIdentifier("identifier_stave", "Cell ID identifier for z/theta module", m_staveID, "stave");

    m_sectorIndex = decoder()->index(m_sectorID);
    m_staveIndex = decoder()->index(m_staveID);
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::position(const CellID& cID) const {
    if (!m_enableTiltedFiber) {
      return FCCSWGridRhoPhiTheta_k4geo::position(cID);
    }
    return fiberPosition(cID);
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::fiberPosition(const CellID& cID) const {
    return projectToFiber(FCCSWGridRhoPhiTheta_k4geo::position(cID), cID);
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::fiberDirection(const CellID& cID) const {
    if (!m_enableTiltedFiber) {
      return normalized(FCCSWGridRhoPhiTheta_k4geo::position(cID));
    }
    Vector3D inner = innerFiberPoint(cID);
    Vector3D outer = outerFiberPoint(cID);
    return normalized(Vector3D(outer.X - inner.X, outer.Y - inner.Y, outer.Z - inner.Z));
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::transverseDistanceToFiber(const Vector3D& globalPosition,
                                                                      const CellID& cID) const {
    Vector3D projected = projectToFiber(globalPosition, cID);
    Vector3D delta(globalPosition.X - projected.X, globalPosition.Y - projected.Y, globalPosition.Z - projected.Z);
    return std::sqrt(dot(delta, delta));
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::projectToFiber(const Vector3D& globalPosition, const CellID& cID) const {
    if (!m_enableTiltedFiber) {
      Vector3D direction = normalized(FCCSWGridRhoPhiTheta_k4geo::position(cID));
      const double distance = dot(globalPosition, direction);
      return Vector3D(direction.X * distance, direction.Y * distance, direction.Z * distance);
    }

    Vector3D start = innerFiberPoint(cID);
    Vector3D direction = fiberDirection(cID);
    Vector3D delta(globalPosition.X - start.X, globalPosition.Y - start.Y, globalPosition.Z - start.Z);
    const double distance = dot(delta, direction);
    return Vector3D(start.X + direction.X * distance, start.Y + direction.Y * distance, start.Z + direction.Z * distance);
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::innerFiberPoint(const CellID& cID) const {
    const double cylR = m_innerR + m_frameThickness;
    const double dphiSec = m_nPhiSectors > 0 ? 2. * M_PI / static_cast<double>(m_nPhiSectors) : 0.;
    const int sector = m_nPhiSectors > 0 ? static_cast<int>(decoder()->get(cID, m_sectorID)) : 0;
    const double sectorPhi = static_cast<double>(sector) * dphiSec;
    return pointFromCylRThetaPhi(cylR, theta(cID), sectorPhi + localPhiInSector(cID));
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::outerFiberPoint(const CellID& cID) const {
    const double tilt = m_moduleTiltAngle * M_PI / 180.;
    const double thetaCenter = theta(cID);
    const double innerCylR = m_innerR + m_frameThickness;
    const double dphiSec = m_nPhiSectors > 0 ? 2. * M_PI / static_cast<double>(m_nPhiSectors) : 0.;
    const double outerRSector = std::abs(std::cos(dphiSec / 2.)) > 1.e-12
                                    ? (innerCylR + m_crystalThickness + m_backSpace) / std::cos(dphiSec / 2.)
                                    : innerCylR + m_crystalThickness + m_backSpace;
    const double radialSpan = outerRSector * std::cos(dphiSec / 2.) - innerCylR;

    double thetaMinModule = thetaCenter;
    double thetaMaxModule = thetaCenter;
    if (m_nZModules > 0 && m_halfZ > 0.) {
      const int stave = static_cast<int>(decoder()->get(cID, m_staveID));
      const double thetaMin = std::atan(innerCylR / m_halfZ);
      const double thetaMax = M_PI - thetaMin;
      const double thetaModule = (thetaMax - thetaMin) / static_cast<double>(m_nZModules);
      const int clampedStave = stave < 0 ? 0 : (stave >= m_nZModules ? m_nZModules - 1 : stave);
      thetaMinModule = thetaMin + static_cast<double>(clampedStave) * thetaModule;
      thetaMaxModule = thetaMinModule + thetaModule;
    }

    const double tiltMin = thetaMinModule < M_PI / 2. ? -tilt : tilt;
    const double tiltMax = thetaMaxModule < M_PI / 2. ? -tilt : tilt;
    const double tiltMinLocal = std::atan(radialSpan / (outerRSector * std::cos(dphiSec / 2.) /
                                                        std::tan(thetaMinModule + tiltMin) -
                                                        innerCylR / std::tan(thetaMinModule)));
    const double tiltMaxLocal = std::atan(radialSpan / (outerRSector * std::cos(dphiSec / 2.) /
                                                        std::tan(thetaMaxModule + tiltMax) -
                                                        innerCylR / std::tan(thetaMaxModule)));

    double heightModule = thetaMaxModule < M_PI / 2. ? m_crystalThickness * std::abs(std::sin(tiltMaxLocal))
                                                     : m_crystalThickness * std::abs(std::sin(tiltMinLocal));
    if ((thetaMaxModule - M_PI / 2.) * (thetaMinModule - M_PI / 2.) < 0.) {
      heightModule = m_crystalThickness;
    }

    const double thetaTilt = thetaCenter < M_PI / 2. ? -tilt : tilt;
    const double tiltLocal = std::atan(radialSpan / (outerRSector * std::cos(dphiSec / 2.) /
                                                     std::tan(thetaCenter + thetaTilt) -
                                                     innerCylR / std::tan(thetaCenter)));
    const double outerCylR = innerCylR + heightModule;
    const double outerZ = innerCylR / std::tan(thetaCenter) + heightModule / std::tan(tiltLocal);
    const double dphiSecForPhi = m_nPhiSectors > 0 ? 2. * M_PI / static_cast<double>(m_nPhiSectors) : 0.;
    const int sector = m_nPhiSectors > 0 ? static_cast<int>(decoder()->get(cID, m_sectorID)) : 0;
    const double sectorPhi = static_cast<double>(sector) * dphiSecForPhi;
    const double outerLocalY = outerLocalYFromInnerPhi(localPhiInSector(cID), innerCylR, outerCylR);
    const double outerPhi = sectorPhi + std::atan2(outerLocalY, outerCylR);
    return Vector3D(outerCylR * std::cos(outerPhi), outerCylR * std::sin(outerPhi), outerZ);
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::pointFromCylRThetaPhi(double cylR, double thetaVal,
                                                                    double phiVal) const {
    return positionFromRThetaPhi(cylR, thetaVal, phiVal);
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::localPhiInSector(const CellID& cID) const {
    if (m_nPhiSectors <= 0) {
      return phi(cID);
    }

    const double dphiSec = 2. * M_PI / static_cast<double>(m_nPhiSectors);
    const int sector = static_cast<int>(decoder()->get(cID, m_sectorID));
    return normalizePhi(phi(cID) - static_cast<double>(sector) * dphiSec);
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::phiBoundaryY(double cylR, bool positiveSide) const {
    const double dphiSec = m_nPhiSectors > 0 ? 2. * M_PI / static_cast<double>(m_nPhiSectors) : 0.;
    const double cosHalfPhi = std::cos(dphiSec / 2.);
    const double safeCosHalfPhi = std::abs(cosHalfPhi) > 1.e-12 ? cosHalfPhi : 1.;
    const double tilt = m_moduleTiltAngle * M_PI / 180.;
    const double outerR = (m_innerR + m_crystalThickness + m_backSpace + 2. * m_frameThickness) / safeCosHalfPhi;

    if (positiveSide) {
      const double slope = (outerR * std::sin(dphiSec / 2. + tilt) - m_innerR * std::tan(dphiSec / 2.)) /
                           (outerR * std::cos(dphiSec / 2. + tilt) - m_innerR);
      const double intercept = m_innerR * std::tan(dphiSec / 2.) - m_innerR * slope;
      return slope * cylR + intercept;
    }

    const double slope = -(outerR * std::sin(dphiSec / 2. - tilt) - m_innerR * std::tan(dphiSec / 2.)) /
                         (outerR * std::cos(dphiSec / 2. - tilt) - m_innerR);
    const double intercept = -m_innerR * std::tan(dphiSec / 2.) - m_innerR * slope;
    return slope * cylR + intercept;
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::outerLocalYFromInnerPhi(double innerPhi, double innerCylR,
                                                                    double outerCylR) const {
    const double innerY = innerCylR * std::tan(innerPhi);
    const double innerYNeg = phiBoundaryY(innerCylR, false) + m_segThickness;
    const double innerYPos = phiBoundaryY(innerCylR, true) - m_segThickness;
    const double outerYNeg = phiBoundaryY(outerCylR, false) + m_segThickness;
    const double outerYPos = phiBoundaryY(outerCylR, true) - m_segThickness;
    const double innerWidth = innerYPos - innerYNeg;
    if (std::abs(innerWidth) < 1.e-12) {
      return outerCylR * std::tan(innerPhi);
    }

    const double fraction = (innerY - innerYNeg) / innerWidth;
    return outerYNeg + fraction * (outerYPos - outerYNeg);
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::normalizePhi(double phiVal) const {
    while (phiVal <= -M_PI) {
      phiVal += 2. * M_PI;
    }
    while (phiVal > M_PI) {
      phiVal -= 2. * M_PI;
    }
    return phiVal;
  }

  Vector3D FCCSWModularGridRhoPhiTheta_k4geo::normalized(const Vector3D& vec) const {
    const double mag2 = dot(vec, vec);
    if (mag2 <= 0.) {
      return Vector3D(0., 0., 1.);
    }
    const double invMag = 1. / std::sqrt(mag2);
    return Vector3D(vec.X * invMag, vec.Y * invMag, vec.Z * invMag);
  }

  double FCCSWModularGridRhoPhiTheta_k4geo::dot(const Vector3D& lhs, const Vector3D& rhs) const {
    return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z;
  }

} // namespace DDSegmentation
} // namespace dd4hep
