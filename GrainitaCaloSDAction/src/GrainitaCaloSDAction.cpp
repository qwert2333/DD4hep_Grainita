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
#include "DD4hep/Segmentations.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include <cmath>


//#define DEBUG

namespace dd4hep {
namespace sim {

  template <>
  Geant4SensitiveAction<GrainitaCaloSDData>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam,
                                                                    DetElement det, Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc), m_collectionName(), m_collectionID(0) {
    declareProperty("ReadoutName", m_readoutName);
    declareProperty("CollectionName", m_collectionName);
    declareProperty("responseInnerCell", m_userData.responseInnerCell);
    declareProperty("responseMapHalfWidth", m_userData.responseMapHalfWidth);
    declareProperty("responseMapBins", m_userData.responseMapBins);
    declareProperty("responseIntercellNeighbor", m_userData.responseIntercellNeighbor);
    declareProperty("responseIntercellDiag", m_userData.responseIntercellDiag);
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


    //Get step local position in cell frame for the response map.
    // Prefer the modular tilted virtual fiber direction when this readout uses it.
    G4ThreeVector ez = HitCellPos.unit();
    auto modularSeg = dynamic_cast<const dd4hep::DDSegmentation::FCCSWModularGridRhoPhiTheta_k4geo*>(m_segmentation->segmentation);
    if (modularSeg && modularSeg->tiltedFiberEnabled()) {
      auto fiberDir = modularSeg->fiberDirection(cellID);
      ez = G4ThreeVector(fiberDir.X, fiberDir.Y, fiberDir.Z).unit();
    }
    G4ThreeVector ref(0,0,1);
    if (std::abs(ez.dot(ref)) > 0.99) ref = G4ThreeVector(1,0,0);

    G4ThreeVector ex = (ref.cross(ez)).unit();
    G4ThreeVector ey = (ez.cross(ex)).unit();    
    G4ThreeVector rel = global - HitCellPos;
    G4ThreeVector local(rel.dot(ex), rel.dot(ey), rel.dot(ez));

    const int responseBins = m_userData.responseMapBins > 0 ? m_userData.responseMapBins : 7;
    const double responseHalfWidth = m_userData.responseMapHalfWidth > 0. ? m_userData.responseMapHalfWidth : 4.;
    const double responseWidth = 2. * responseHalfWidth;
    int step_idx = floor((local.x() + responseHalfWidth) / responseWidth * responseBins);
    int step_idy = floor((local.y() + responseHalfWidth) / responseWidth * responseBins);
    if(step_idx<0) step_idx = 0;
    if(step_idx>=responseBins) step_idx = responseBins - 1;
    if(step_idy<0) step_idy = 0;
    if(step_idy>=responseBins) step_idy = responseBins - 1;
    const int responseIndex = step_idx * responseBins + step_idy;
    const double response = responseIndex < static_cast<int>(m_userData.responseInnerCell.size()) ? m_userData.responseInnerCell[responseIndex] : 0.;
    double step_E = (1 + response/100.) * aStep->GetTotalEnergyDeposit();


#ifdef DEBUG
    auto phiID = m_segmentation->decoder()->get(cellID, "phi");
    auto thetaID = m_segmentation->decoder()->get(cellID, "theta");
    auto rhoID = m_segmentation->decoder()->get(cellID, "rho");
    std::cout<<"--> Step global position: ("<<global.x()<<", "<<global.y()<<", "<<global.z()<<") ";
    std::cout<<" (theta, phi, rho) = "<<"("<<global.theta()<<", "<<global.phi()<<", "<<global.mag()<<") "<<std::endl;
    std::cout<<"    local position: ("<<local.x()<<", "<<local.y()<<", "<<local.z()<<") "<<std::endl;
    std::cout<<"  phiID: "<<phiID<<", thetaID "<<thetaID<<", rhoID "<<rhoID<<", cellID "<<cellID<<std::endl;
    std::cout<<"  Cell position: ("<<HitCellPos.x()<<", "<<HitCellPos.y()<<", "<<HitCellPos.z()<<std::endl;
    std::cout<<" (theta, phi, rho) = "<<"("<<HitCellPos.theta()<<", "<<HitCellPos.phi()<<", "<<HitCellPos.mag()<<") "<<std::endl;
#endif

    // Create the hits and accumulate contributions from multiple steps
    //
    Geant4HitCollection* coll = collection(m_collectionID);
    Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(cellID); // the hit

    if (!hit) { // if the hit does not exist yet, create it
      hit = new Geant4Calorimeter::Hit();

      hit->cellID = cellID;
      hit->position = HitCellPos; // this should be assigned only once
      hit->energyDeposit = step_E;

      // Add calo hit contributions
      //
      // Crete the first contribution associated to this hit
      Geant4Calorimeter::Hit::Contribution contrib;
      contrib.trackID = aStep->GetTrack()->GetTrackID();
      contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      contrib.deposit = step_E;
      contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
      contrib.x = global.x();
      contrib.y = global.y();
      contrib.z = global.z();
      hit->truth.emplace_back(contrib);

      coll->add(cellID, hit); // add the hit to the hit collection
    } else {                 // if the hit exists already, increment its fields
      hit->energyDeposit += step_E;

      // Add calo hit contributions
      //
      // Add a new contribution associated to this hit
      Geant4Calorimeter::Hit::Contribution contrib;
      contrib.trackID = aStep->GetTrack()->GetTrackID();
      contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      contrib.deposit = step_E;
      contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
      contrib.x = global.x();
      contrib.y = global.y();
      contrib.z = global.z();
      hit->truth.emplace_back(contrib);
    }

    return true;
  } // end of Geant4SensitiveAction::process() method specialization

} // namespace sim
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(GrainitaCaloSDAction)

//**************************************************************************

