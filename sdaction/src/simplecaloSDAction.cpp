/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "DD4hep/Segmentations.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include <cmath>

#define DEBUG

namespace dd4hep {
namespace sim {
  class simplecaloSDData {
    // Constructor and destructor
    //
  public:
    simplecaloSDData() = default;
    ~simplecaloSDData() = default;

  public:
    Geant4Sensitive* sensitive{};
  };
} // namespace sim
} // namespace dd4hep

namespace dd4hep {
namespace sim {

  // Function template specialization of Geant4SensitiveAction class.
  // Define actions
  template <>
  void Geant4SensitiveAction<simplecaloSDData>::initialize() {
    m_userData.sensitive = this;
    m_hitCreationMode = HitCreationFlags::DETAILED_MODE;
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Define collections created by this sensitivie action object
  template <>
  void Geant4SensitiveAction<simplecaloSDData>::defineCollections() {
    std::string ROname = m_sensitive.readout().name();
    std::cout<<"  Readout name: "<<ROname<<std::endl;
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>(ROname);
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Method that accesses the G4Step object at each track step.
  template <>
  bool Geant4SensitiveAction<simplecaloSDData>::process(const G4Step* aStep, G4TouchableHistory* /*history*/) {

#ifdef DEBUG
    std::cout << "-------------------------------" << std::endl;
    std::cout << "--> simplecalo: track info: " << std::endl;
    std::cout << "----> Track #: " << aStep->GetTrack()->GetTrackID() << " "
              << "Step #: " << aStep->GetTrack()->GetCurrentStepNumber() << " "
              << "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " "
              << std::endl;
    std::cout << "--> simplecalo:: position info(mm): " << std::endl;
    std::cout << "----> x: " << aStep->GetPreStepPoint()->GetPosition().x()
              << " y: " << aStep->GetPreStepPoint()->GetPosition().y()
              << " z: " << aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
    std::cout << "--> simplecalo: particle info: " << std::endl;
    std::cout << "----> Particle " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " "
              << "Dep(MeV) " << aStep->GetTotalEnergyDeposit() << " "
              << "Mat " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " "
              << "Vol " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " " << std::endl;
#endif

    auto decoder = m_sensitive.readout().idSpec().decoder();
    auto VolID = volumeID(aStep);

#ifdef DEBUG
    auto SystemID = decoder->get(VolID, "system");
    auto StaveID = decoder->get(VolID, "stave");
    auto SectorID = decoder->get(VolID, "sector");
    std::cout<< "--> Volume ID: "<<SystemID<<"  "<<StaveID<<"  "<<SectorID<<std::endl;
#endif

    G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector global = (aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition() )/2.;
    G4ThreeVector local = theTouchable->GetHistory()->GetTopTransform().TransformPoint(global);
    dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter,
                         local.y() * dd4hep::millimeter / CLHEP::millimeter,
                         local.z() * dd4hep::millimeter / CLHEP::millimeter);
    dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                          global.y() * dd4hep::millimeter / CLHEP::millimeter,
                          global.z() * dd4hep::millimeter / CLHEP::millimeter);

    auto cellID = m_segmentation->cellID(loc, glob, VolID);
    auto phiID = m_segmentation->decoder()->get(cellID, "phi");
    auto thetaID = m_segmentation->decoder()->get(cellID, "theta");

#ifdef DEBUG
    std::cout<<"--> Step global position: ("<<global.x()<<", "<<global.y()<<", "<<global.z()<<") ";
    std::cout<<" (theta, phi) = "<<"("<<global.theta()<<", "<<global.phi()<<") "<<std::endl;
    std::cout<<"    local position: ("<<local.x()<<", "<<local.y()<<", "<<local.z()<<") "<<std::endl;
    std::cout<<"  phiID: "<<phiID<<", thetaID "<<thetaID<<std::endl;
#endif

    // Create the hits and accumulate contributions from multiple steps
    //
    Geant4HitCollection* coll = collection(m_collectionID);
    Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(cellID); // the hit

    if (!hit) { // if the hit does not exist yet, create it
      hit = new Geant4Calorimeter::Hit();
      hit->cellID = cellID;
      Position HitCellPos(global.x(), global.y(), global.z());
      hit->position = HitCellPos; // this should be assigned only once
      hit->energyDeposit = aStep->GetTotalEnergyDeposit();

      // Add calo hit contributions
      //
      // Crete the first contribution associated to this hit
      Geant4Calorimeter::Hit::Contribution contrib;
      contrib.trackID = aStep->GetTrack()->GetTrackID();
      contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      contrib.deposit = aStep->GetTotalEnergyDeposit();
      contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
      contrib.x = HitCellPos.x();
      contrib.y = HitCellPos.y();
      contrib.z = HitCellPos.z();
      hit->truth.emplace_back(contrib);

      coll->add(VolID, hit); // add the hit to the hit collection
    } else {                 // if the hit exists already, increment its fields
      hit->energyDeposit += aStep->GetTotalEnergyDeposit();

      // Add calo hit contributions
      //
      // Add a new contribution associated to this hit
      Geant4Calorimeter::Hit::Contribution contrib;
      contrib.trackID = aStep->GetTrack()->GetTrackID();
      contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      contrib.deposit = aStep->GetTotalEnergyDeposit();
      contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
      Position HitCellPos(global.x(), global.y(), global.z());
      contrib.x = HitCellPos.x();
      contrib.y = HitCellPos.y();
      contrib.z = HitCellPos.z();
      hit->truth.emplace_back(contrib);
    }

    return true;
  } // end of Geant4SensitiveAction::process() method specialization

} // namespace sim
} // namespace dd4hep

//--- Factory declaration
namespace dd4hep {
namespace sim {
  typedef Geant4SensitiveAction<simplecaloSDData> SimpleCaloSDAction;
}
} // namespace dd4hep
DECLARE_GEANT4SENSITIVE(SimpleCaloSDAction)

//**************************************************************************
