#include "detectorSegmentations/FCCSWModularRhoPhiTheta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWModularRhoPhiTheta_k4geo::FCCSWModularRhoPhiTheta_k4geo(const std::string& cellEncoding) : FCCSWGridPhiTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWModularRhoPhiTheta_k4geo::FCCSWModularRhoPhiTheta_k4geo(const BitFieldCoder* decoder) : FCCSWGridPhiTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWModularRhoPhiTheta_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWModularRhoPhiTheta_k4geo";
    _description = "Rho-Phi-theta segmentation in the global coordinates";

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("grid_rho", "Grid size in rho", m_grid_rho, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("offset_R", "Offset in R", m_offsetR, 0., SegmentationParameter::LengthUnit, true);
    registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");

    m_thetaIndex = decoder()->index(fieldNameTheta());
    m_phiIndex = decoder()->index(fieldNamePhi());
    m_rhoIndex = decoder()->index(m_rhoID);
  }


  /// determine the local based on the cell ID
  Vector3D FCCSWModularRhoPhiTheta_k4geo::position(const CellID& cID) const {
    return positionFromRThetaPhi(rho(cID)*std::sin(theta(cID)) , theta(cID), phi(cID));
  }


  /// determine the cell ID based on the position
  CellID FCCSWModularRhoPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                         const VolumeID& vID) const {
    CellID cID = vID;
    double lTheta = thetaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    double lrho = sqrt( radiusFromXYZ(globalPosition)*radiusFromXYZ(globalPosition) + globalPosition.Z*globalPosition.Z );

    //dd4hep::Segmentation::positionToBin: int(floor((position + 0.5 * cellSize - offset) / cellSize));
    decoder()->set(cID, m_thetaIndex, positionToBin(lTheta, gridSizeTheta(), offsetTheta()));
    decoder()->set(cID, m_phiIndex, positionToBin(lPhi, gridSizePhi(), offsetPhi() ));
    decoder()->set(cID, m_rhoIndex, positionToBin(lrho, m_grid_rho, m_offsetR/std::sin(lTheta)));
    return cID;
  }


  /// determine the distance to IP rho based on the cell ID
  double FCCSWModularRhoPhiTheta_k4geo::rho(const CellID cID) const {
    CellID rhoValue = decoder()->get(cID, m_rhoIndex);
    double m_theta = theta(cID);
    // dd4hep::Segmentation::binToPosition: bin * cellSize + offset;
    return binToPosition(rhoValue, m_grid_rho, m_offsetR/std::sin(m_theta) );
  }
  


} // namespace DDSegmentation
} // namespace dd4hep
