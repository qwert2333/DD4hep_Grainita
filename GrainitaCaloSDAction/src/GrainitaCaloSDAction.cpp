#ifndef GRAINITACALOSDACTION_C
#define GRAINITACALOSDACTION_C 1
#endif

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

#include "GrainitaCaloSDAction.h"
#include "detectorSegmentations/FCCSWModularGridRhoPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "DD4hep/Segmentations.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>


// #define DEBUG

namespace dd4hep {
namespace sim {

  template <>
  Geant4SensitiveAction<GrainitaCaloSDData>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam,
                                                                    DetElement det, Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc), m_collectionName(), m_collectionID(0) {
    declareProperty("ReadoutName", m_readoutName);
    declareProperty("CollectionName", m_collectionName);
    declareProperty("responseFuncSlope", m_userData.slope);
    declareProperty("responseFuncIntersect", m_userData.intersect);
    declareProperty("responseFuncX0", m_userData.x0);
    declareProperty("responseFuncAttLength", m_userData.AttLength);
    declareProperty("responseFuncNorm", m_userData.norm);
    declareProperty("neighborCellSize", m_userData.neighborCellSize);
    declareProperty("fiberAttenuationLength", m_userData.fiberAttenuationLength);
    declareProperty("outerRadius", m_userData.outerRadius);


    InstanceCount::increment(this);

    m_hitCreationMode = HitCreationFlags::DETAILED_MODE;
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Define actions
  template <>
  void Geant4SensitiveAction<GrainitaCaloSDData>::initialize() {
    m_userData.sensitive = this;
    m_hitCreationMode = HitCreationFlags::DETAILED_MODE;
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Define collections created by this sensitivie action object
  template <>
  void Geant4SensitiveAction<GrainitaCaloSDData>::defineCollections() {
    std::string ROname = m_sensitive.readout().name();
    std::cout<<"  Readout name: "<<ROname<<std::endl;
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>(ROname);
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Method that accesses the G4Step object at each track step.
  template <>
  bool Geant4SensitiveAction<GrainitaCaloSDData>::process(const G4Step* aStep, G4TouchableHistory* /*history*/) {

#ifdef DEBUG
    std::cout << "-------------------------------" << std::endl;
    std::cout << "--> GrainitaCalo: track info: " << std::endl;
    std::cout << "----> Track #: " << aStep->GetTrack()->GetTrackID() << " "
              << "Step #: " << aStep->GetTrack()->GetCurrentStepNumber() << " "
              << "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " "
              << std::endl;
    std::cout << "--> GrainitaCalo:: position info(mm): " << std::endl;
    std::cout << "----> x: " << aStep->GetPreStepPoint()->GetPosition().x()
              << " y: " << aStep->GetPreStepPoint()->GetPosition().y()
              << " z: " << aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
    std::cout << "--> GrainitaCalo: particle info: " << std::endl;
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
    dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                          global.y() * dd4hep::millimeter / CLHEP::millimeter,
                          global.z() * dd4hep::millimeter / CLHEP::millimeter);

    auto cellID = m_segmentation->cellID(glob, glob, VolID);
    auto hitpos_dd4hep = m_segmentation->position(cellID); // in cm
    G4ThreeVector HitCellPos(hitpos_dd4hep.x()*dd4hep::centimeter/dd4hep::millimeter, hitpos_dd4hep.y()*dd4hep::centimeter/dd4hep::millimeter, hitpos_dd4hep.z()*dd4hep::centimeter/dd4hep::millimeter );

    auto modularSeg = dynamic_cast<const dd4hep::DDSegmentation::FCCSWModularGridRhoPhiTheta_k4geo*>(m_segmentation->segmentation);
    auto phiThetaSeg = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(m_segmentation->segmentation);
    const int phiIndex = decoder->index("phi");
    const int thetaIndex = decoder->index("theta");
    const int phiBins = phiThetaSeg ? phiThetaSeg->phiBins() : 0;
    const int currentPhiID = static_cast<int>(decoder->get(cellID, phiIndex));
    const int currentThetaID = static_cast<int>(decoder->get(cellID, thetaIndex));
    const int neighborSize = std::max(1, m_userData.neighborCellSize);
    const int neighborRadius = neighborSize / 2;

    std::vector<CellID> cellIDvec;
    std::vector<G4ThreeVector> cellPosVec;
    std::vector<G4double> responseVec;

    auto dd4hepPositionToG4 = [](const dd4hep::DDSegmentation::Vector3D& pos) {
      return G4ThreeVector(pos.X * dd4hep::centimeter / dd4hep::millimeter,
                           pos.Y * dd4hep::centimeter / dd4hep::millimeter,
                           pos.Z * dd4hep::centimeter / dd4hep::millimeter);
    };

    auto addCell = [&](CellID id) {
      auto pos = m_segmentation->position(id);
      cellIDvec.push_back(id);
      cellPosVec.push_back(dd4hepPositionToG4(pos));
    };

    addCell(cellID);
    for (int dPhi = -neighborRadius; dPhi <= neighborRadius; ++dPhi) {
      for (int dTheta = -neighborRadius; dTheta <= neighborRadius; ++dTheta) {
        if (dPhi == 0 && dTheta == 0) {
          continue;
        }

        int neighborPhiID = currentPhiID + dPhi;
        if (phiBins > 0) {
          while (neighborPhiID >= phiBins) {
            neighborPhiID -= phiBins;
          }
          while (neighborPhiID < 0) {
            neighborPhiID += phiBins;
          }
        } else {
          std::cout << "Error: phiBins is not well defined: " << phiBins << ". Cannot apply periodic boundary conditions." << std::endl;
          continue;
        }

        const int neighborThetaID = currentThetaID + dTheta;

        CellID neighborCellID = cellID;
        decoder->set(neighborCellID, phiIndex, neighborPhiID);
        decoder->set(neighborCellID, thetaIndex, neighborThetaID);
        addCell(neighborCellID);
      }
    }

    auto transverseDistance = [&](CellID id, const G4ThreeVector& cellPos) {
      G4ThreeVector axis = cellPos.unit();
      if (modularSeg) {
        auto fiberDir = modularSeg->fiberDirection(id);
        axis = G4ThreeVector(fiberDir.X, fiberDir.Y, fiberDir.Z).unit();
      }
      G4ThreeVector rel = global - cellPos;
      G4double transverse2 = rel.mag2() - rel.dot(axis) * rel.dot(axis);
      return std::sqrt(std::max(0., transverse2));
    };

    responseVec.reserve(cellIDvec.size());
    for (std::size_t i = 0; i < cellIDvec.size(); ++i) {
      responseVec.push_back(m_userData.lightResponse(transverseDistance(cellIDvec[i], cellPosVec[i])));
    }


#ifdef DEBUG
    auto phiID = m_segmentation->decoder()->get(cellID, "phi");
    auto thetaID = m_segmentation->decoder()->get(cellID, "theta");
    auto rhoID = m_segmentation->decoder()->get(cellID, "rho");
    std::cout<<"--> Step global position: ("<<global.x()<<", "<<global.y()<<", "<<global.z()<<") ";
    std::cout<<" (theta, phi, rho) = "<<"("<<global.theta()<<", "<<global.phi()<<", "<<global.mag()<<") "<<std::endl;
    std::cout<<"  phiID: "<<phiID<<", thetaID "<<thetaID<<", rhoID "<<rhoID<<", cellID "<<cellID<<std::endl;
    std::cout<<"  Cell position: ("<<HitCellPos.x()<<", "<<HitCellPos.y()<<", "<<HitCellPos.z()<<std::endl;
    std::cout<<" (theta, phi, rho) = "<<"("<<HitCellPos.theta()<<", "<<HitCellPos.phi()<<", "<<HitCellPos.mag()<<") "<<std::endl;
    std::cout<<"  Neighbor cell count: "<<cellIDvec.size()<<std::endl;
    responseSum = 0;
    for(size_t i=0; i<cellIDvec.size(); ++i) {
      phiID = m_segmentation->decoder()->get(cellIDvec[i], "phi");
      thetaID = m_segmentation->decoder()->get(cellIDvec[i], "theta");
      rhoID = m_segmentation->decoder()->get(cellIDvec[i], "rho");
      std::cout<<"  Neighbor cell "<<i<<": phiID "<<phiID<<", thetaID "<<thetaID<<", rhoID "<<rhoID<<", cellID "<<cellIDvec[i]<<", response "<<responseVec[i]<<std::endl;
      responseSum += responseVec[i];
    }
    std::cout<<"  Response sum: "<<responseSum<<std::endl;
#endif

    // Create the hits and accumulate contributions from multiple steps
    //
    Geant4HitCollection* coll = collection(m_collectionID);
    for (std::size_t i = 0; i < cellIDvec.size(); ++i) {
      const CellID hitCellID = cellIDvec[i];
      G4double longitudinalAttenuation = 1.;
      if (m_userData.fiberAttenuationLength > 0.) {
        longitudinalAttenuation = std::exp(-(m_userData.outerRadius-global.perp()) / m_userData.fiberAttenuationLength);
      }
      const G4double step_E = responseVec[i] * longitudinalAttenuation * aStep->GetTotalEnergyDeposit();
      Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(hitCellID); // the hit

      if (!hit) { // if the hit does not exist yet, create it
        hit = new Geant4Calorimeter::Hit();
        hit->cellID = hitCellID;
        hit->position = cellPosVec[i]; // this should be assigned only once
        hit->energyDeposit = step_E;
        coll->add(hitCellID, hit); // add the hit to the hit collection
      } else {                 // if the hit exists already, increment its fields
        hit->energyDeposit += step_E;
      }

      // Add calo hit contributions
      if(i == 0) { // only save the truth info for the central cell to avoid duplication
        Geant4Calorimeter::Hit::Contribution contrib;
        contrib.trackID = aStep->GetTrack()->GetTrackID();
        contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        contrib.deposit = aStep->GetTotalEnergyDeposit();
        contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
        contrib.x = global.x();
        contrib.y = global.y();
        contrib.z = global.z();
        hit->truth.emplace_back(contrib);
      }
    }

    return true;
  } // end of Geant4SensitiveAction::process() method specialization

} // namespace sim
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(GrainitaCaloSDAction)

//**************************************************************************

