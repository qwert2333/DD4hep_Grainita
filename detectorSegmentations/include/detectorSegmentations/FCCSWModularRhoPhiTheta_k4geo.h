#ifndef DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

/** FCCSWModularRhoPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWModularRhoPhiTheta_k4geo.h
 * FCCSWModularRhoPhiTheta_k4geo.h
 *
 *  Segmentation in rho, theta and phi.
 *  Based on FCCSWGridPhiTheta_k4geo, addition of Rho segmentation
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWModularRhoPhiTheta_k4geo : public FCCSWGridPhiTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWModularRhoPhiTheta_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWModularRhoPhiTheta_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWModularRhoPhiTheta_k4geo() = default;

    /**  Determine the global position based on the cell ID.
     *   @warning This segmentation has no knowledge of radius, so radius = 1 is taken into calculations.
     *   @param[in] aCellId ID of a cell.
     *   return Position (radius = 1).
     */
    virtual Vector3D position(const CellID& aCellID) const override;
    /**  Determine the cell ID based on the position.
     *   @param[in] aLocalPosition (not used).
     *   @param[in] aGlobalPosition position in the global coordinates.
     *   @param[in] aVolumeId ID of a volume.
     *   return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;
 
  private:
    /// Get rho from cellID
    double rho(const CellID cID) const;
    /// the grid size in rho
    double m_grid_rho;

    /// the coordinate offset in rho / R
    double m_offsetR;   // In cylinder case: offset depends on R not rho. 
    double m_offsetrho; 

    /// the field name used for rho
    std::string m_rhoID;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for theta
    int m_thetaIndex = -1;
    /// the field index used for phi
    int m_phiIndex = -1;
    /// the field index used for rho
    int m_rhoIndex = -1;



  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_GRIDRHOPHITHETA_H */
