// usage
// root -b
// .L plot_cell_energy_sum.C
// plot_cell_energy_sum(path_to_file)
// NB: in this case, a simple ROOT command in the interpreter (events->Draw("Sum$(simplecaloRO.energy)")) would give you the same information but the goal here is to teach how to properly analyze podio files (needed for more advanced analyses)

#include <string>

#include "podio/Reader.h"

#include "edm4hep/SimCalorimeterHitCollection.h"

int plot_cell_energy_sum(std::string input_file) {

  // Use the podio reader to open the file (it could also be opened with plain ROOT but that would not preserve the podio functionalities)
  auto reader = podio::makeReader(input_file);

  // Prepare a ROOT TH1 to store the cell energy sum
  TH1F* cell_energy_sum_th1 = new TH1F("cell_energy_sum_th1", ";Cell energy sum [GeV]; Number of events", 100, 0, 4.0);

  // Fill the TH1 with the cell energy sum
  for (size_t i = 0; i < reader.getEvents(); ++i) {
    auto event = reader.readNextEvent();
    auto& calo_cells = event.get<edm4hep::SimCalorimeterHitCollection>("simplecaloRO");
    float total_energy = 0.;
    for (const auto& calo_cell : calo_cells) {
      total_energy += calo_cell.getEnergy();
    }
    cell_energy_sum_th1->Fill(total_energy);
  }

  // Fit and draw the histogram
  TCanvas* cell_energy_sum_canvas = new TCanvas("cell_energy_sum_canvas");
  float fit_range_min = cell_energy_sum_th1->GetXaxis()->GetBinCenter(cell_energy_sum_th1->GetMaximumBin()) - 3 * cell_energy_sum_th1->GetRMS();
  float fit_range_max = cell_energy_sum_th1->GetXaxis()->GetBinCenter(cell_energy_sum_th1->GetMaximumBin()) + 3 * cell_energy_sum_th1->GetRMS();
  cell_energy_sum_th1->Fit("gaus", "SQ", "", fit_range_min, fit_range_max);
  gStyle->SetOptFit(1111); // to get the fit result displayed
  cell_energy_sum_th1->Draw();
  // Print the canvas in png with a name based on the input file name
  size_t pos = input_file.find(".root");
  if (pos != std::string::npos) {
    input_file.replace(pos, 5, "_cell_energy_sum.png");
  }
  cell_energy_sum_canvas->Print(input_file.c_str());

  return 0;
}
