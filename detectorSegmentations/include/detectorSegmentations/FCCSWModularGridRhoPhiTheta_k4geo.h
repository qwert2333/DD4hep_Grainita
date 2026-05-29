#ifndef DETECTORSEGMENTATIONS_FCCSWMODULARGRIDRHOPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_FCCSWMODULARGRIDRHOPHITHETA_K4GEO_H

#include "detectorSegmentations/FCCSWGridRhoPhiTheta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

  class FCCSWModularGridRhoPhiTheta_k4geo : public FCCSWGridRhoPhiTheta_k4geo {
  public:
    FCCSWModularGridRhoPhiTheta_k4geo(const std::string& cellEncoding);
    FCCSWModularGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder);
    virtual ~FCCSWModularGridRhoPhiTheta_k4geo() = default;

    virtual Vector3D position(const CellID& cID) const override;

    Vector3D fiberPosition(const CellID& cID) const;
    Vector3D fiberDirection(const CellID& cID) const;
    double transverseDistanceToFiber(const Vector3D& globalPosition, const CellID& cID) const;
    Vector3D projectToFiber(const Vector3D& globalPosition, const CellID& cID) const;
    bool tiltedFiberEnabled() const { return m_enableTiltedFiber; }

  private:
    void commonSetup();
    Vector3D innerFiberPoint(const CellID& cID) const;
    Vector3D outerFiberPoint(const CellID& cID) const;
    Vector3D pointFromCylRThetaPhi(double cylR, double theta, double phi) const;
    double localPhiInSector(const CellID& cID) const;
    double phiBoundaryY(double cylR, bool positiveSide) const;
    double outerLocalYFromInnerPhi(double innerPhi, double innerCylR, double outerCylR) const;
    double normalizePhi(double phi) const;
    Vector3D normalized(const Vector3D& vec) const;
    double dot(const Vector3D& lhs, const Vector3D& rhs) const;

    bool m_enableTiltedFiber = true;
    double m_innerR = 0.;
    double m_halfZ = 0.;
    double m_crystalThickness = 0.;
    double m_backSpace = 0.;
    double m_frameThickness = 0.;
    double m_segThickness = 0.;
    double m_moduleTiltAngle = 0.;
    int m_nPhiSectors = 1;
    int m_nZModules = 1;
    std::string m_sectorID = "sector";
    std::string m_staveID = "stave";
    int m_sectorIndex = -1;
    int m_staveIndex = -1;
  };

} // namespace DDSegmentation
} // namespace dd4hep

#endif
