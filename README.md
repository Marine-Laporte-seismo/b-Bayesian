# b_Bayesian_McMC – 2026

## Probabilistic Estimate of b-value Temporal Variations

**Last update:** 27/10/2026  
**Author:** Marine Laporte  
**Contact:** marinelaporte@icm.csic.es  
**Version:** 2026 – faster code, fully McMC

---
> ⚠️ **This is not the latest version of b-Bayesian.**
> This version uses a continuous likelihood: use it only if your catalog has **non-rounded magnitudes**.
> For rounded catalogs and any new study, use the latest and recommended version:
> 👉 [`b-Bayesian_binned_McMC`](../b-Bayesian_binned_McMC_2026) – full McMC, faster, accounts for magnitude discretization.

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

where Φ is the standard normal cumulative distribution and M_min is the minimum magnitude of the whole catalog.

> **What is new in this McMC version:**
> *b-Bayesian_McMC* explores, within a single reversible-jump McMC, the number and positions of the temporal discontinuities **and** the values of β, μ and σ in each segment.

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

---

## 🔄 Running the Inversion

### 1. Set the input file and parameters

Edit `/inversion/b_Bayesian.f90`:

- Set your dataset path:
  `Magnitude_file = "your/data/path"`
- Adjust the McMC parameters (`it_max`, `it_burnin`, `it_thin`, `max_layers`, `T_bins`, …)
- Set the prior bounds for **β**, **μ** and **σ** (start with broad priors)

### 2. Compile and run

In the `/inversion/` folder:

```bash
make
./run
```

### 3. Check convergence

Run **3–4 chains with different random seeds** and check that they give the same distribution of the number of discontinuities and the same profiles. If they do not, increase `it_max`.

---

## 📂 Output Files

b-Bayesian writes the following text files to `/outputs/`:

| File | Description |
|---|---|
| `input.txt` | Copy of the input catalog (original times) |
| `hypercube.txt` | Run parameters: prior bounds and number of bins for μ, β, σ; iterations; time offset and duration; bin widths |
| `prop.txt` | Number of post-burn-in models with a discontinuity in each time bin (size: `T_bins`) |
| `NC_tot.txt` | Number of discontinuities of the stored models (one value every `it_thin` iterations, burn-in included) |
| `Likelihood_b.txt` | Histogram of β for each time bin (rows: β bins, columns: time bins) |
| `Likelihood_Mu.txt` | Histogram of μ for each time bin (rows: μ bins, columns: time bins) |
| `Likelihood_Sigma.txt` | Histogram of σ for each time bin (rows: σ bins, columns: time bins) |
| `Acceptances.txt` | Number of accepted and proposed models and acceptance rate (%) for each move type: geometry (birth, death, move) and GR parameters (β, μ, σ) |
| `AR_evolution.txt` | Acceptance rates of each move type over successive windows of iterations |
| `Beta_tot.txt`, `Mu_tot.txt`, `Sigma_tot.txt`, `Times_tot.txt` | Parameters and discontinuity times of the current model at regular intervals (rows: segments, columns: storage times; `99` = no segment) |

> The `Likelihood_*` files are **counts**: divide each column by its sum to get a probability density.

---

## 📈 Plotting Figures

The script `output_make_figures.py` (in `/figures`) reads the inversion outputs and generates the figures. It needs `imports_outfiles_.py` in the same folder.

### Steps

1. In the main section of the script, set the run name and check the paths to the output files.
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
- **Fig5:** b-value of each segment along the chain (convergence check)

---

## 🧪 Synthetic Tests

- `generation_synthetic_etas.py`: piecewise-constant catalog (b, μ, σ) with ETAS clustering

Writes the catalog for the Fortran code and a text file with the true parameters.

---

## 📝 License

Academic use only. Please contact the author for permission to use, modify, or distribute this code.
