#ifndef GRAINITACALO_MATRIX_SDACTION_H
#define GRAINITACALO_MATRIX_SDACTION_H 1

#include "DDG4/Geant4SensDetAction.h"

#include <cmath>
#include <string>
#include <vector>

namespace dd4hep {
namespace sim {

  class GrainitaCaloMatrixSDData {
  public:
    GrainitaCaloMatrixSDData() = default;
    ~GrainitaCaloMatrixSDData() = default;

    Geant4Sensitive* sensitive{};
    int rawCollectionID = 0;
    std::string rawCollectionName = "GrainitaCalorimeterHitsRaw";
    int neighborCellSize = 5;
    double fiberAttenuationLength = 1000.; // mm
    double outerRadius = 2645.; // mm

     // Light response function parameters
    bool useLightResponseFunction = true;
    // Light response function:
    //    y = slope*x + intersect (x<x0)
    //    y = norm*exp(-x/AttLength) (x>=x0)
    double x0 = 0.856;
    double AttLength = 3.4; // mm
    double slope = 0.93;
    double intersect = 0.206;
    double norm = 1. / std::exp(-1. * x0 / AttLength);


    
    double lightResponse(double distance) const {
      if (distance < x0) {
        return slope * distance + intersect;
      }
      return norm * std::exp(-distance / AttLength);
    }

  };

  using GrainitaCaloMatrixSDAction = Geant4SensitiveAction<GrainitaCaloMatrixSDData>;

} // namespace sim
} // namespace dd4hep

#endif
