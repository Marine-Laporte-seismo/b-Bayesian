# b_Bayesian_binned – v2 – 2025
 
## Probabilistic Estimate of b-value Temporal Variations
 
**Last update:** 27/10/2026  
**Author:** Marine Laporte  
**Contact:** [marinelaporte@icm.csic.es](mailto:marinelaporte@icm.csic.es)  
**Version:** 2 – magnitude binning and faster computation
 
-> ⚠️ **This is not the latest version of b-Bayesian.**
> This version is kept for reproducibility of previous studies (2025).
> For any new study, use the latest and recommended version:
> 👉 [`b-Bayesian_binned_McMC`](../b-Bayesian_binned_McMC_2026) – full McMC, faster, accounts for magnitude discretization.--
 
## 📚 Citation
 
If you use this code, please cite:
 
> Laporte, M., Durand, S., Bodin, T., Gardonio, B., & Marsan, D. (2025). *b-bayesian: The full Probabilistic estimate of b-value Temporal variations for non-truncated catalogs.* Journal of Geophysical Research: Solid Earth, 130, e2024JB029973. https://doi.org/10.1029/2024JB029973
 
---
 
## 🧭 Description
 
**b-Bayesian** performs a **transdimensional Bayesian inversion** of the temporal evolution of:
 
- the **b-value** of the Gutenberg-Richter law (the code works with β = b·ln(10))
- **μ (mu)** and **σ (sigma)**, which describe the detectability of earthquakes with a cumulative normal (error function) detection model (Ogata & Katsura, 1993)
The magnitude distribution in each time segment is
 
p(m) ∝ Φ((m − μ)/σ) · β · exp(−β (m − M_min)),  m ≥ M_min
 
where Φ is the standard normal cumulative distribution.
 
A reversible-jump McMC explores the number and positions of the temporal discontinuities. For each proposed segmentation, β, μ and σ are integrated out in each segment by Monte Carlo sampling of their uniform priors (`n_draws` random draws), which gives the marginal likelihood of the segmentation and the marginal distributions of β, μ and σ.
 
> **What is new in this *binned* version:**
> the likelihood is computed on the magnitude histogram of each segment (bins of width `delta_bin`, usually 0.1) instead of event by event. Its cost depends on the number of magnitude bins, not on the number of events, which strongly reduces the computation time for catalogs of several thousand events. It is also the exact likelihood for catalogs whose magnitudes are rounded.
 
---
 
## ⚙️ Installation
 
### Dependencies
 
- **To run b-Bayesian:** a Fortran 90 compiler (set it in the `Makefile`)
- **To generate figures:** Python 3 with
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `cmcrameri`
```bash
pip install numpy scipy matplotlib cmcrameri
```
 
---
 
## 📥 Input File
 
A plain text file with **two columns**, separated by spaces:
 
```
time   magnitude
```
 
- Events **must be sorted in time** (the program stops otherwise).
- Time can be in any unit (days, decimal years, …). The program subtracts the time of the first event, so the inversion runs on [0, duration].
- Magnitudes are expected to be rounded to `delta_bin` (e.g. 0.1).
---
 
## 🔄 Running the Inversion
 
### 1. Set the input file and parameters
 
Edit `/inversion/b_Bayesian.f90`:
 
- **Data**
  - `Magnitude_file = "your/data/path"`
  - `delta_bin`: magnitude rounding of the catalog (e.g. `0.1d0`)
  - `min_bin`: lower edge of the first magnitude bin, also the lower integration bound. It must be ≤ the smallest magnitude. For a catalog rounded to 0.1, use (smallest magnitude) − 0.05 so that the bins are centred on the rounded values.
  - `max_bin`: must be > the largest magnitude
  - `T_bins`: number of time bins of the outputs
- **McMC:** `it_max`, `it_burnin`, `it_thin`, `it_std`
- **Priors:** bounds of the uniform priors of **β**, **μ** and **σ** (start with broad priors) and the number of bins of their marginal distributions (`Beta_bins`, `Mu_bins`, `Sigma_bins`)
- **Monte Carlo integration:** `n_draws`. If the program prints `ATTENTION : seulement ... tirages efficaces`, `n_draws` is too small for at least one segment: increase it or narrow the priors.
### 2. Compile and run
 
In the `/inversion/` folder:
 
```bash
make
./run
```
 
### 3. Check convergence
 
Run **3–4 chains with different random seeds** (`ra`) and check that they give the same distribution of the number of discontinuities and the same profiles. If they do not, increase `it_max`.
 
---
 
## 📂 Output Files
 
b-Bayesian writes 7 text files to `/outputs/`:
 
| File | Description |
|---|---|
| `input.txt` | Copy of the input catalog, with times starting at 0 |
| `prior.txt` | Run parameters: one line per parameter for μ, β, σ (`min max n_bins`), then `it_max n_draws T_bins computation_time` |
| `Probability_discontinuities.txt` | Number of post-burn-in models with a discontinuity in each time bin (size: `T_bins`). Divide by the number of stored models to get a probability. |
| `Number_of_discontinuities.txt` | Number of discontinuities at each iteration (size: `it_max`, burn-in included) |
| `Marginal_b.txt` | Marginal distribution of **β** for each time bin (rows: β bins, columns: time bins) |
| `Marginal_Mu.txt` | Marginal distribution of μ for each time bin (rows: μ bins, columns: time bins) |
| `Marginal_Sigma.txt` | Marginal distribution of σ for each time bin (rows: σ bins, columns: time bins) |
| `Acceptances_Rates.txt` | 3 × 3 matrix. Columns: birth, death, move. Rows: accepted, proposed, acceptance rate (fraction) |
 
> The `Marginal_*` files are sums over the stored models: normalise each column to get a probability density. Bin `k` of a marginal covers [min + (k−1)·L/N, min + k·L/N), where L is the width of the prior and N the number of bins.
 
---
 
## 📈 Plotting Figures
 
The script `output_figures_bMuSigma.py` (in `/figures`) reads the inversion outputs and generates the figures.
 
### Steps
 
1. In the main section of the script, set:
   - the run name and the paths to the output files
   - `year0`: date of the first event of the original catalog (0 to keep relative times)
   - `it_burnin` and `it_thin`: same values as in the Fortran code
2. Run the script:
```bash
python3 output_figures_bMuSigma.py
```
 
### Output
 
Creates a new folder named `Run_day_month_year_hourh` in the folder where the script is run, containing:
 
- **Fig0:** Magnitudes over time with the probability of discontinuities
- **Fig1:** 4-panel summary:
  - (Top-left) Probability of discontinuities over time bins
  - (Top-right) Accepted vs. proposed models (birth/death/move)
  - (Bottom-left) Number of discontinuities along the chain
  - (Bottom-right) Histogram of the number of discontinuities after burn-in
- **Fig2:** Probability density of **μ** over time
- **Fig3:** Probability density of the **b-value** over time
- **Fig4:** Probability density of **σ** over time
---
 
## 📝 License
 
Academic use only. Please contact the author for permission to use, modify, or distribute this code.

