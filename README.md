# Polarization_analysis

This repository provides tools for plotting and comparing Standard Model (SM) vs. Effective Field Theory (EFT) predictions across different polarization states in WZ production.

## Features

- Reads ROOT files containing physics data (currently with Pythia shower)
- Computes histograms normalized by absolute cross-section
- Plots:
  - Leading variable distributions (e.g. `el_pt`, `jet_pt`, `MET`)
  - Invariant mass distributions (`M_ee`, `M_mumu`, `M_epmup`)
  - Ratio plots (EFT / SM) including per-polarization breakdowns for better visualisation
## Installation

dependencies can be found in requirements.txt file:

pip install -r requirements.txt

## Usage

Activate your virtual environment and run:

```bash
python3 main.py --var el_pt        # Plot a specific variable
python3 main.py --all              # Plot all standard ROOT variables
