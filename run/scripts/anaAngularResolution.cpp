#include "podio/Reader.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

class Result {
  public:
    double theta;
    double phi;
    double totalE;
    bool valid;
};

Result GetWeightedPosition(
    vector<double>* hit_x,
    vector<double>* hit_y,
    vector<double>* hit_z,
    vector<double>* hit_E,
    int weightMethod = 0; 
)
{
    // =========================================
    // apply cuts
    // =========================================

    double sumE = 0;
    TVector3 dir(0,0,0);
    double ave_x = 0.;
    double ave_y = 0.;
    double ave_z = 0.;
    double sumW = 0;

    // Get the total energy for log-weight calculation.
    for(size_t i=0;i<hit_E->size();i++)
        sumE += hit_E->at(i);

    // Get shower position 
    for(size_t i=0;i<hit_E->size();i++)
    {
        double x = hit_x->at(i);
        double y = hit_y->at(i);
        double z = hit_z->at(i);
        double E = hit_E->at(i);

        dir += E * TVector3(x, y, z);
        //Log weighted
        double weight = std::max(0., log(hit_E->at(i)) - log(sumE) + 4);
        sumW += weight;
        ave_x += x * weight;
        ave_y += y * weight;
        ave_z += z * weight;
    }

    if(weightMethod == 0)
      dir = dir.Unit();
    else(weightMethod == 1)
      dir.SetXYZ(ave_x/sumW, ave_y/sumW, ave_z/sumW);
    else{
      std::cout<<"Unknown weight method: 0 for linear weight, 1 for log weight. "<<std::endl;
    }

    Result r;

    if(sumE <= 0)
    {
        r.valid = false;
        return r;
    }

    r.theta = dir.Theta();
    r.phi = dir.Phi();
    r.totalE = sumE;
    r.valid = true;

    return r;
}  


vector<double> getAngularRes(TString s_filename = "CaloSim_gamma_5GeV_theta85-95_phi80-100.root"){

  gStyle->SetOptFit(1011);

  auto reader = podio::makeReader(s_filename);
  int Nevnts = reader.getEvents();
  cout<<"File name: "<<s_filename<<", total events: "<<Nevnts<<endl;

  TH1F *h_En = new TH1F("h_En", "", 100, 4.2, 5.2);
  TH1F *h_theta = new TH1F("h_theta", "", 100, -3, 3);
  TH1F *h_phi = new TH1F("h_phi", "", 100, -3, 3);
  TH2F *h2_hitmap = new TH2F("h2_hitmap", "h2_hitmap", 100, -2700, 2700, 100, -2700, 2700);
  TH2F *h2_theta = new TH2F("h2_theta", "", 100, 1.567, 1.574, 100, -3, 3);
  TH2F *h2_phi = new TH2F("h2_phi", "", 100, 1.567, 1.574, 100, -3, 3);
  TH2F *h2_hit_E_delta = new TH2F("h2_hit_E_delta", "", 100, -10, 2, 100, -6, 1);
  int aveNhit = 0;

  for (size_t i = 0; i < Nevnts; ++i) {
    auto event = reader.readNextEvent();
    auto& mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
    auto& calo_hits = event.get<edm4hep::SimCalorimeterHitCollection>("GrainitaCalorimeterHits");

    TVector3 mcP3;
    for(const auto&mcp : mcparticles){
      if(mcp.parents_size()==0) { mcP3.SetXYZ(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z); break; }
    }

    // get hit information
    vector<double> hit_x;
    vector<double> hit_y;
    vector<double> hit_z;
    vector<double> hit_E;
    for(const auto& hit : calo_hits){
      hit_x.push_back(hit.getPosition().x);
      hit_y.push_back(hit.getPosition().y);
      hit_z.push_back(hit.getPosition().z);
      hit_E.push_back(hit.getEnergy());
    }

    Result res = GetWeightedPosition(&hit_x, &hit_y, &hit_z, &hit_E, 0);
    if(!res.valid) continue;

    h_En->Fill(res.totalE);
    h_theta->Fill( (res.theta - mcP3.Theta()) * 1000); // in mrad
    h_phi->Fill( (res.phi - mcP3.Phi()) * 1000 ); // in mrad

    h2_theta->Fill(res.theta, (res.theta - mcP3.Theta())*1000.);
    h2_phi->Fill(res.phi, (res.phi - mcP3.Phi())*1000.);
  }


  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  c1->Divide(2,1);
  c1->cd(1);
  h_theta->GetXaxis()->SetTitle("#Delta#theta / mrad");
  h_theta->Fit("gaus", "Q");
  h_theta->Draw();
  c1->cd(2);
  h_phi->GetXaxis()->SetTitle("#Delta#phi / mrad");
  h_phi->Fit("gaus", "Q");
  h_phi->Draw();
  c1->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  c2->Divide(2,1);
  c2->cd(1);
  h2_theta->GetXaxis()->SetTitle("#theta_{rec} / rad");
  h2_theta->GetYaxis()->SetTitle("#Delta #theta / mrad");
  h2_theta->Draw("colz");
  c2->cd(2);
  h2_phi->GetXaxis()->SetTitle("#phi_{rec} / rad");
  h2_phi->GetYaxis()->SetTitle("#Delta #phi / mrad");
  h2_phi->Draw("colz");
  c2->Draw();
  //h2_hitmap->Draw("colz");

  vector<double> results;
  results.push_back( h_theta->GetFunction("gaus")->GetParameter(2) );
  results.push_back( h_phi->GetFunction("gaus")->GetParameter(2) );
  results.push_back( h_theta->GetFunction("gaus")->GetParError(2) );
  results.push_back( h_phi->GetFunction("gaus")->GetParError(2) );
  return results; 
}


void anaAngularResolution(){

  gStyle->SetPaintTextFormat("4.3f");
  const int Npoint = 8;
  double pitch[Npoint] = {3., 5., 7., 10., 15., 20., 30., 40.};
  double thetaRes[Npoint] = {0.};
  double phiRes[Npoint] = {0.};
  double thetaResErr[Npoint] = {0.};
  double phiResErr[Npoint] = {0.};
  Selection method[5];
  method[1].useWindow = true;
  method[2] = method[1];
  method[2].Ecut = 0.01;
  method[3] = method[2];
  method[3].firstLayerOnly = true;
  method[4] = method[3];
  method[4].Ecut = 0.05;

  //for(int ii=0; ii<Npoint; ii++){
  //  //int ii=2; 
  //  vector<double> tmp_res = getAngularRes( (int)pitch[ii], method[2] );
  //  thetaRes[ii] = tmp_res[0];
  //  phiRes[ii] = tmp_res[1];
  //  thetaResErr[ii] = tmp_res[2];
  //  phiResErr[ii] = tmp_res[3];
  //}
  //for(int ii=0; ii<Npoint; ii++) cout<<pitch[ii]<<"\t"<<(thetaRes[ii]+phiRes[ii])/2*1000.<<'\t'<<sqrt(pow(thetaResErr[ii],2)+pow(phiResErr[ii],2))/2*1000.<<endl;


  //for(int ii=0; ii<5; ii++){
   int ii = 2; 
   vector<double> tmp_res = getAngularRes( 7, method[ii] );
   thetaRes[ii] = tmp_res[0];
   phiRes[ii] = tmp_res[1];
   thetaResErr[ii] = tmp_res[2];
   phiResErr[ii] = tmp_res[3];
  //}

  for(int ii=0; ii<Npoint; ii++) cout<<(thetaRes[ii]+phiRes[ii])/2*1000.<<'\t'<<sqrt(pow(thetaResErr[ii],2)+pow(phiResErr[ii],2))/2*1000.<<endl;

  //TGraph *gr_AngRes = new TGraph(Npoint, pitch, thetaRes);
  //gr_AngRes->SetLineWidth(2);
  //gr_AngRes->Draw("ALP");

  TCanvas *c3 = new TCanvas();
  TH1F *gr_AngRes_method = new TH1F("gr_AngRes_method", "", 5, 0, 5);
  for(int i=0; i<5; i++){ 
    gr_AngRes_method->SetBinContent(i+1, (thetaRes[i]+phiRes[i])/2*1000.);
    gr_AngRes_method->SetBinError(i+1, sqrt(pow(thetaResErr[i],2)+pow(phiResErr[i],2))/2*1000.);
  }
  gr_AngRes_method->GetXaxis()->SetBinLabel(1, "No cut");
  gr_AngRes_method->GetXaxis()->SetBinLabel(2, " + 20 cm");
  gr_AngRes_method->GetXaxis()->SetBinLabel(3, " + 10 MeV");
  gr_AngRes_method->GetXaxis()->SetBinLabel(4, " + 1st layer");
  gr_AngRes_method->GetXaxis()->SetBinLabel(5, " + 50 MeV");
  gr_AngRes_method->GetYaxis()->SetTitle("Angular resolution / mrad");
  gr_AngRes_method->Draw("LP E0 TEXT0");

}




