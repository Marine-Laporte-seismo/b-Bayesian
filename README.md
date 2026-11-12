# b_Bayesian_binned_McMC – 2026
 
## Probabilistic Estimate of b-value Temporal Variations
 
> ✅ **This is the latest and recommended version of b-Bayesian.**
> It is the fastest version and the only one that combines a full McMC exploration while considering discrete magnitudes. 
**Use this version for any new study.**
 
**Last update:** 28/10/2026  
**Author:** Marine Laporte  
**Contact:** marinelaporte@icm.csic.es / marinelaporte.ensg@gmail.com 
**Version:** 2026 – full McMC with binned magnitudes (version free of rights without numerical recipes)
 
---
 
## 📚 Citation
 
If you use this code, please cite:
 
> Laporte, M., Durand, S., Bodin, T., Gardonio, B., & Marsan, D. (2025). *b-bayesian: The full Probabilistic estimate of b-value Temporal variations for non-truncated catalogs.* Journal of Geophysical Research: Solid Earth, 130, e2024JB029973. https://doi.org/10.1029/2024JB029973
 
---
 
## 🗂️ Which version should I use?
 
| Version | Exploration of β, μ, σ | Magnitudes | Status |
|---|---|---|---|
| **`b-Bayesian_binned_McMC`** (this folder) | full McMC | **binned** (rounded catalogs) | ✅ **latest, recommended** | FAST
| [`b-Bayesian_McMC`](../b-Bayesian_McMC) | full McMC | continuous | for catalogs with non-rounded magnitudes | FAST
| [`b-Bayesian_binned`](../b-Bayesian_binned) | Monte Carlo integration in each segment | binned | previous version (2025) | SLOW
| [`b-Bayesian`](../b-Bayesian) | as in Laporte et al. (2025) | continuous | original version | VERY SLOW
 
---
 
## 🧭 Description
 
**b-Bayesian** performs a **transdimensional Bayesian inversion** of the temporal evolution of:
 
- the **b-value** of the Gutenberg-Richter law (the code works with β = b·ln(10))
- **μ (mu)** and **σ (sigma)**, which describe the detectability of earthquakes with a cumulative normal (error function) detection model (Ogata & Katsura, 1993)
The magnitude distribution in each time segment is
 
p(m) ∝ Φ((m − μ)/σ) · β · exp(−β (m − M₀)),  m ≥ M₀
 
where Φ is the standard normal cumulative distribution and M₀ is the lower edge of the first magnitude bin.
 
A single reversible-jump McMC explores the number and positions of the temporal discontinuities **and** the values of β, μ and σ in each segment.
 
### Why this version?
 
- **Magnitude discretization.** Catalog magnitudes are rounded (usually to 0.1). Instead of treating them as exact values, the likelihood integrates the magnitude distribution over each bin [m − Δm/2, m + Δm/2]. The bins are centred on the rounded magnitudes and built automatically from the data. This removes the bias that a continuous likelihood introduces on μ and σ, especially when σ is small.
- **Computation time.** The likelihood of a segment is computed on its magnitude histogram: its cost depends on the number of magnitude bins (a few tens), not on the number of events. Combined with the full McMC (no Monte Carlo integration over the priors), this makes it the fastest version for catalogs of thousands to hundreds of thousands of events.
- **Robustness.** The likelihood is numerically stable over the whole prior range (no `-Infinity` or `NaN`), and a minimum number of events per segment can be imposed (`N_min_seg`) to avoid spurious discontinuities isolating a handful of events.
---
 
## ⚙️ Installation
 
### Dependencies
 
- **To run b-Bayesian:** a Fortran 90 compiler (set it in the `Makefile`)
- **To generate figures:** Python 3 with
  - `numpy`
  - `scipy`
  - `pandas`
  - `matplotlib`
  - `mpmath`
  - `plotly`
  - `cmcrameri`
```bash
pip install numpy scipy pandas matplotlib mpmath plotly cmcrameri
```
 
---
 
## 📥 Input File
 
A plain text file with **two columns**, separated by spaces:
 
```
time   magnitude
```
 
- Events **must be sorted in time** (the program stops otherwise).
- Time can be in any unit (days, decimal years, …). Absolute times (e.g. 2006.5) are accepted: the program subtracts the first event time internally and writes this offset to `hypercube.txt`.
- Magnitudes are considered to be **rounded to `delta_bin`** (e.g. 0.1). If some magnitudes are not multiples of `delta_bin`, the program prints a warning: check that `delta_bin` matches the resolution of your catalog. 
---
 
## 🔄 Running the Inversion
 
### 1. Set the input file and parameters
 
Edit `/inversion/b-Bayesian_binned_McMC.f90`:
 
- **Data**
  - `Magnitude_file = "your/data/path"`
  - `delta_bin`: magnitude resolution of the catalog (e.g. `0.1d0`). The magnitude bins (`min_bin`, `N_bins`) are then computed automatically.
  - `T_bins`: number of time bins of the outputs
- **McMC:** `it_max`, `it_burnin` (must be smaller than `it_max`, typically 40–50 %), `it_thin`, `it_std`
- **Priors:** bounds of the uniform priors of **β**, **μ** and **σ** (start with broad priors)
- **Proposals:** step sizes `beta_stdi`, `mu_stdi`, `sigma_stdi` (aim for acceptance rates of 25–50 % in `Acceptances.txt`) and `*_std_birth` for births and deaths
- **Segments:** `max_layers` (maximum number of discontinuities) and `N_min_seg`, the minimum number of events per segment (50 by default). Note that this constraint makes the effective prior on the number of discontinuities slightly decreasing.
### 2. Compile and run
 
In the `/inversion/` folder:
 
```bash
make
./run
```
 
To check the code (array bounds, NaN, uninitialized variables) on a short run:
 
```bash
make debug
./run
```
 
Recompile with `make clean && make` before production runs: the debug version is 5–20 times slower.
 
### 3. Check convergence
 
Run **3–4 chains with different random seeds** (`ra`) and check that they give the same distribution of the number of discontinuities and the same profiles. If they do not, increase `it_max`.
 
---
 
## 📂 Output Files
 
b-Bayesian writes the following text files to `/outputs/`:
 
| File | Description |
|---|---|
| `input.txt` | Copy of the input catalog (original times) |
| `hypercube.txt` | Run parameters, one line each: μ, β, σ (prior bounds, number of bins, initial value, initial step); iterations; time offset and duration; bin widths of the β, μ, σ histograms; magnitude binning (`delta_bin`, `min_bin`, `N_bins`) |
| `prop.txt` | Number of post-burn-in models with a discontinuity in each time bin (size: `T_bins`) |
| `NC_tot.txt` | Number of discontinuities of the stored models (one value every `it_thin` iterations, burn-in included) |
| `Likelihood_b.txt` | Histogram of β for each time bin (rows: β bins, columns: time bins) |
| `Likelihood_Mu.txt` | Histogram of μ for each time bin (rows: μ bins, columns: time bins) |
| `Likelihood_Sigma.txt` | Histogram of σ for each time bin (rows: σ bins, columns: time bins) |
| `Acceptances.txt` | Number of accepted and proposed models and acceptance rate (%) for each move type: geometry (birth, death, move) and GR parameters (β, μ, σ) |
| `AR_evolution.txt` | Acceptance rates of each move type over successive windows of iterations |
| `Beta_tot.txt`, `Mu_tot.txt`, `Sigma_tot.txt`, `Times_tot.txt` | Parameters and discontinuity times of the current model at regular intervals (rows: segments, columns: storage times; `99` = no segment) |
 
> The `Likelihood_*` files are **counts**: divide each column by its sum to get a probability density. `prop.txt` divided by the number of stored post-burn-in models gives the probability of a discontinuity in each time bin.
 
---
 
## 📈 Plotting Figures
 
The script `output_make_figures.py` (in `/figures`) reads the inversion outputs and generates the figures. It needs `imports_outfiles_.py` in the same folder.
 
### Steps
 
1. In the main section of the script, set the run name, check the paths to the output files.
2. Run the script:
```bash
python3 output_make_figures.py
```
 
### Output
 
Creates a new folder named `Run_day_month_year_hourh` in the folder where the script is run, containing:
 
- **Fig0:** Magnitudes over time with the probability of discontinuities
- **Fig1:** 4-panel summary:
  - (Top-left) Probability of discontinuities over time bins
  - (Top-right) Accepted vs. proposed models for each move type
  - (Bottom-left) Number of discontinuities along the chain
  - (Bottom-right) Histogram of the number of discontinuities
- **Fig2:** Probability density of **μ** over time
- **Fig3:** Probability density of the **b-value** over time
- **Fig4:** Probability density of **σ** over time
---
 
## 🧪 Synthetic Tests
 
- `generation_synthetic_etas_binned.py`: piecewise-constant catalog (b, μ, σ) with ETAS clustering
It writes the catalog for the Fortran code and a text file with the true parameters. To test this version in realistic conditions, it rounds the synthetic magnitudes to `delta_bin`.
 
---
 
## 📝 License
 
Academic use only. Please contact the author for permission to use, modify, or distribute this code.

