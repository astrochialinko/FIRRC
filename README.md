# FIRRC: COSMOS Radio & Sub-mm Data Analysis

## Goals

This repository contains Python scripts, modules, and utilities for processing and analyzing astronomical survey data in the COSMOS field, focusing on radio (VLA 1.4 GHz and 3 GHz) and sub-millimeter (JCMT 450 µm) observations to investigate whether or not the infrared-radio correlation varies with redshift

## Installation

1. Clone the repository

   ```
   git clone <repository_url>
   cd FIRRC
   ```

2. Create and activate the conda environment

   ```
   conda create -n irrc38 python=3.8.11
   conda activate irrc38
   pip install -r requirements.txt
   ```

## Dependencies

Defined in `requirements.txt`, includes:

- Python 3.8.11
- numpy, scipy, astropy, photutils, gvar, pandas, matplotlib
- bdsf (PyBDSF) for radio source finding

## Data Requirements

Most scripts depend on external input data that isn't included in this repository.
Before running the analysis, make sure you have:

- Raw FITS images for VLA at 1.4 GHz and 3 GHz
  - `Data/COSMOS/Image/VLA/cosmos_vla_3GHz_2017image_uJy.fits`
  - `Data/COSMOS/Image/VLA/cosmos_vla_3GHz_2017rms_uJy.fits`
  - `Data/COSMOS/Image/VLA/cosmos_vla_1d4GHz_XS_2021image_uJy.fits`
  - `Data/COSMOS/Image/VLA/cosmos_vla_1d4GHz_XS_2021rms_uJy.fits`
- Catalog files in FITS or text format (for JCMT, IRAC, MIPS, VLA) in:
  - `Data/COSMOS/Catalog/` (e.g. `catalog_vla/*.fits`, `catalog_jcmt/*.fits`)
  - A master list of catalog filenames should be listed in `Script/catalog.txt`
- Primary-beamed corrected images and any additional RMS or weight maps required by `script_02_make_radec.py`
- Cross-match reference catalogs for matching in `script_01_make_crossmatch.py`

Once you place the data in these paths, the scripts will read the files according to the paths defined in `Script/path.py`.

## Data Organization

- Input Data: `Data/COSMOS/Image/VLA/`, `Data/COSMOS/Catalog/`
- Cross-matched Tables: `Data/COSMOS/Catalog/CrossMatch/`
- Stacking Outputs: `Data/COSMOS/Image/Stacking/`
- Figures: `Figures/`

## Usage

Before running scripts, change into the `Script/` directory:

```
cd Script
```

1. Catalog Preparation and Cross-matching

   ```
   python script_00_remodel_catalog.py
   python script_01_make_crossmatch.py
   ```

2. Coordinate and Flux Corrections

   ```
   python script_02_make_radec.py
   ```

3. Stacking Analyses

   ```
   python script_03_make_stacking_studies.py
   python script_04_bdsf_stacking_Coord_3GHz.py
   ```

4. Spectral Index Computation

   ```
   python script_05_cal_spectral_index_irac.py
   python script_cal_spectral_index.py
   ```

5. Plotting Final Results

   ```
   python script_06_plot_stacking_result.py
   ```

6. Notebooks for Plotting

   Typically, you can launch a Jupyter server at the project root and then open any notebook via the browser file browser:

   ```
   jupyter notebook
   ```

   Alternatively, to open a specific notebook directly, you can pass its path:

   ```
   jupyter notebook script_plot_alpha.ipynb
   jupyter notebook script_plot_diffR.ipynb
   ```

## Included files

```
FIRRC/
├── Figures/
│   └── *.pdf                     # diagnostic figures
├── README.md                     # Project overview and usage guide
├── requirements.txt              # Python environment and dependencies (python=3.8.11)
└── Script/
    ├── ClassFluxDensity.py       # Measure peak, integrated, and aperture flux
    ├── ClassSpecIndex.py         # Compute spectral indices from multi-frequency catalogs
    ├── ClassStacking.py          # Perform image stacking analyses
    ├── ClassStatistic.py         # Statistical utilities: bootstrapping, binning, errors
    ├── FuncTableCrossmatch.py    # Cross-match catalogs via KD-tree and angular matching
    ├── path.py                   # Dataset and output path definitions
    ├── catalog.txt               # Master list of catalog files for batch processing
    ├── script_00_remodel_catalog.py           # Convert and rename raw catalogs to CSV/FITS
    ├── script_01_make_crossmatch.py           # Cross-match multi-wavelength catalogs
    ├── script_02_make_radec.py                # Apply primary-beam correction, prepare RA/Dec catalogs
    ├── script_03_make_stacking_studies.py     # Generate stacked images for different source subsets
    ├── script_04_bdsf_stacking_Coord_3GHz.py  # Run PyBDSF on stacks, extract source parameters
    ├── script_05_cal_spectral_index_irac.py   # Calculate spectral indices on stacked datasets
    ├── script_cal_spectral_index.py           # Additional spectral index routines and batch runner
    ├── script_06_plot_stacking_result.py      # Produce final diagnostic plots
    ├── script_casa_subimage_VLA.py            # Create CASA sub-images and export FITS
    ├── script_make_matchtable_for_latex.ipynb # Notebook to generate LaTeX-ready match tables
    ├── script_image.ipynb                     # Exploratory notebook for image processing
    ├── script_plot_alpha.ipynb                # Notebook for spectral index visualization
    └── script_plot_diffR.ipynb                # Notebook for cross-match radius and flux comparisons
```
