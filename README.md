# es-crustal-thermal-state

Radiogenic heat production, magnetic source depths and Monte Carlo conductive geotherms for Espírito Santo (SE Brazil, Araçuaí orogen).

Code and derived data for the manuscript:

> do Vale Pena, C. *A cold crust beneath a once-hot orogen: radiogenic heat production, magnetic source depths and conductive geotherms in Espírito Santo, Araçuaí orogen, SE Brazil.* Submitted to the Journal of South American Earth Sciences.

## Contents

```
run_all.py runs scripts 01–13 in order (outputs go to results/ and figures/)
scripts/   processing and analysis
results/   derived grids and tables used in the paper
figures/   Figures 1–5 of the paper (300 dpi)
data/      place the input data here (not redistributed, see below)
```

| Script | What it does | Main outputs |
|---|---|---|
| `01_heat_production_50m.py` | Heat production A = ρ(9.52 U + 2.56 Th + 3.48 K)·10⁻⁵ (Rybach, 1988) from the Project 1093 grids, clipped to ES | `calor_radiogenico_uWm3_ES_v2.tif` (50 m) |
| `02_colatina_corridor.py` | 30-km corridor on the Ecoporanga–Vitória axis; NNW structures inside it | `lineamento_colatina_*.gpkg` |
| `03_heat_production_125m_and_reference_geotherm.py` | Area-average to 125 m; deterministic reference geotherm | `A0_125m.tif` |
| `04_heat_production_by_lithology.py` | A₀ by lithological group, inside/outside the corridor, block permutation test | `A0_por_litologia.csv` (Table 1) |
| `05_emag2v3_crop.py` | Crop of EMAG2v3 (4 km upward continued) around ES | `emag2v3_recorte.nc` |
| `bouligand.py` | Radially averaged spectra and the fractal model of Bouligand et al. (2009) | – |
| `06_aeromagnetic_fractal_fit.py` | β and depth to top from 150 km aeromagnetic windows | printed values |
| `07_centroid_estimates.py` | Centroid estimates with and without fractal correction (EMAG2v3) | printed values |
| `08_ghfdb_extract.py` | GHFDB entries between 15.5–23°S and 38.5–44°W | `ghfdb_regiao.csv` |
| `09_monte_carlo.py` | 4000-draw Monte Carlo conditioned on the heat flow at Cachoeiro de Itapemirim | `mc_parametros.csv`, `mc_*_500m.tif` |
| `10_sensitivity.py` | Sensitivity of the posterior to the heat-flow constraint | `sensitivity_heatflow.csv` (Table 3) |
| `11–13_*.py` | Figures 1–5 | `figures/*.png` |

Comments and printed messages in some scripts are in Portuguese.

## Input data (not included)

| Data | Source |
|---|---|
| Project 1093 gamma-ray grids (K, eTh, eU) and magnetic anomaly | Serviço Geológico do Brasil (SGB-CPRM), GeoSGB. Grids regridded at 50 m by the author with Oasis montaj. |
| Geological, structural and state-boundary shapefiles | SGB-CPRM, Geologia e recursos minerais do estado do Espírito Santo, 1:400,000 (Vieira and Menezes, 2015) |
| EMAG2v3, upward continued GeoTIFF | NOAA/NCEI, https://doi.org/10.7289/V5H70CVX |
| IHFC Global Heat Flow Database, Release 2024 (v2026.03) | GFZ Data Services, https://doi.org/10.5880/fidgeo.2024.014 |

Set the paths at the top of `run_all.py` or with the environment variables `ES_DATA`, `EMAG2_TIF` and `GHFDB_TXT` (or `GHFDB_ZIP`; the .txt is then extracted automatically).

## Environment

Python 3.13 with the packages in `requirements.txt` (versions used for the paper).

```
pip install -r requirements.txt
python run_all.py              # everything
python run_all.py --from 9     # resume from script 09
python run_all.py --only 12    # a single script
```

A log is written to `run_all.log`.

The Monte Carlo uses fixed random seeds, so the tables are reproduced exactly.

## Results files

- `A0_125m.tif`: surface radiogenic heat production (µW m⁻³), 125 m, UTM 24S (EPSG:32724).
- `mc_q0_media_500m.tif`, `mc_q0_dp_500m.tif`: posterior mean and standard deviation of surface heat flow (mW m⁻²).
- `mc_TMoho_media_500m.tif`: posterior mean Moho temperature (°C).
- `mc_z580_P50_500m.tif`: posterior median depth of the 580 °C isotherm (km).
- `mc_parametros.csv`: Monte Carlo draws and likelihood weights.
- `ghfdb_regiao.csv`: GHFDB subset (CC BY 4.0, Global Heat Flow Data Assessment Group et al., 2024).

## Citation

If you use this code or data, please cite the paper and this repository (see `CITATION.cff`).

## Licence

Code: MIT (see `LICENSE`). Derived data: CC BY 4.0.

## Acknowledgements

Part of the doctoral thesis of C. do Vale Pena at INPE, with a CAPES fellowship (Finance Code 001). Developed with the assistance of Claude (Anthropic); see the AI declaration in the paper.
