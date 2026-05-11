# Forecasting Dynamic Panel Models with Shrinkage of Fixed Effects

Replication files for the paper **"Forecasting Dynamic Panel Models with Shrinkage of Fixed Effects"** by Magdalena Cornejo and Walter Sosa-Escudero.

## Repository structure

- `code/simulations/01_monte_carlo_simulation.R`: Monte Carlo simulation for one design.
- `code/simulations/02_make_mse_figures.R`: figures based on the Monte Carlo output tables.
- `code/empirical/03_empirical_application_compustat.R`: empirical application using Compustat data from WRDS.
- `data/`: local data folder. Raw Compustat data should not be committed.
- `output/figures/`: generated figures.
- `output/tables/`: generated tables.

## Data availability

The Monte Carlo results are reproducible using the simulation scripts in this repository.

The empirical application uses Compustat Fundamentals Annual data accessed through WRDS. These data are proprietary and cannot be redistributed. Researchers with WRDS access can reproduce the empirical application by placing the WRDS extract locally at:

```text
data/raw/df0_WRDS.csv
```

The script `code/empirical/03_empirical_application_compustat.R` documents the variable construction and sample restrictions.
