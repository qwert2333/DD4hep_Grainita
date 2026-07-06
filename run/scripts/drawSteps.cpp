#include "podio/Reader.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

void drawSteps(){

  auto reader = podio::makeReader("CaloSim_gamma_5GeV_theta85-95_phi80-100.root");


  TH2D *h2_stepmap = new TH2D("h2_stepmap", "", 400, -400, 400, 400, 2000, 2800);
  TH2D *h2_stepmap_RZ = new TH2D("h2_stepmap_RZ", "", 200, -400, 400, 200, 2000, 2800);
  TH2D *h2_hitmap = new TH2D("h2_hitmap", "", 200, -400, 400, 200, 2000, 2800);
  TH2D *h2_hitmap_RZ = new TH2D("h2_hitmap_RZ", "", 200, -400, 400, 200, 2000, 2800);
  int Nevnts = reader.getEvents();
  cout<<"Total events: "<<Nevnts<<endl;

  for (size_t i = 0; i < Nevnts; ++i) {
    auto event = reader.readNextEvent();
    auto& calo_hits = event.get<edm4hep::SimCalorimeterHitCollection>("GrainitaCalorimeterHitsRaw");
    for(const auto& hit : calo_hits){
      //if(hit.getPosition().y>2600 ) continue;
      //if(fabs(hit.getPosition().z)>100) continue;
      h2_hitmap->Fill(hit.getPosition().x, hit.getPosition().y, hit.getEnergy());
      h2_hitmap_RZ->Fill( hit.getPosition().z, sqrt(hit.getPosition().x*hit.getPosition().x+hit.getPosition().y*hit.getPosition().y ), hit.getEnergy() );
      auto calo_steps = hit.getContributions();
      for(const auto& step : calo_steps){
        h2_stepmap->Fill(step.getStepPosition().x, step.getStepPosition().y, step.getEnergy());
        h2_stepmap_RZ->Fill(step.getStepPosition().z, sqrt(step.getStepPosition().x*step.getStepPosition().x+step.getStepPosition().y*step.getStepPosition().y), step.getEnergy());
      }
    }    
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 800);
  c1->Divide(2,1);
  c1->cd(1);
  h2_stepmap->GetXaxis()->SetTitle("x / mm");
  h2_stepmap->GetYaxis()->SetTitle("y / mm");
  h2_stepmap->GetYaxis()->SetTitleOffset(1.8);
  h2_stepmap->Draw("colz");
  h2_hitmap->Draw("same");

  c1->cd(2);
  h2_stepmap_RZ->GetXaxis()->SetTitle("z / mm");
  h2_stepmap_RZ->GetYaxis()->SetTitle("R / mm");
  h2_stepmap_RZ->GetYaxis()->SetTitleOffset(1.8);
  h2_stepmap_RZ->Draw("colz");
  h2_hitmap_RZ->Draw("same");

}

