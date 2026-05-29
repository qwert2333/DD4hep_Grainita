#ifndef GRAINITACALOSDACTION_H
#define GRAINITACALOSDACTION_H 1

#include "DDG4/Geant4SensDetAction.h"

#include <vector>

namespace dd4hep {
namespace sim {

  class GrainitaCaloSDData {
  public:
    GrainitaCaloSDData() = default;
    ~GrainitaCaloSDData() = default;

    Geant4Sensitive* sensitive{};

    std::vector<double> responseInnerCell = std::vector<double>(49, 0.);
    double responseMapHalfWidth = 4.;
    int responseMapBins = 7;
    double responseIntercellNeighbor = 0.2;
    double responseIntercellDiag = 0.1;
  };

  using GrainitaCaloSDAction = Geant4SensitiveAction<GrainitaCaloSDData>;

} // namespace sim
} // namespace dd4hep

#endif
