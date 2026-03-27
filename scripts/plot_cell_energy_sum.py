import argparse
import ROOT
# prevent ROOT to display things on the flight
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from podio import root_io

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="Path to the input file")
args = parser.parse_args()

# Use the podio reader to open the file (it could also be opened with plain ROOT but that would not preserve the podio functionalities)
reader = root_io.Reader(args.input_file)

# Prepare a ROOT TH1 to store the cell energy sum
cell_energy_sum_th1 = ROOT.TH1F("cell_energy_sum_th1", ";Cell energy sum [GeV]; Number of events", 100, 0, 4.0) 

# Fill the TH1 with the cell energy sum
for event in reader.get("events"):
    calo_cells = event.get("simplecaloRO")
    total_energy = sum([calo_cell.getEnergy() for calo_cell in calo_cells])
    cell_energy_sum_th1.Fill(total_energy)

# Fit and draw the histogram
cell_energy_sum_canvas = ROOT.TCanvas("cell_energy_sum_canvas")
fit_range_min = cell_energy_sum_th1.GetXaxis().GetBinCenter(cell_energy_sum_th1.GetMaximumBin()) - 3 * cell_energy_sum_th1.GetRMS()
fit_range_max = cell_energy_sum_th1.GetXaxis().GetBinCenter(cell_energy_sum_th1.GetMaximumBin()) + 3 * cell_energy_sum_th1.GetRMS()
cell_energy_sum_th1.Fit("gaus", "SQ", "", fit_range_min, fit_range_max)
ROOT.gStyle.SetOptFit(1111) # to get the fit result displayed
cell_energy_sum_th1.Draw()
cell_energy_sum_canvas.Print(args.input_file.replace(".root", "_cell_energy_sum.png"))





