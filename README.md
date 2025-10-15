# Monte Carlo Analysis of Carbon Emissions from Household Showers

## Author: Muhammad Umar IJAZ

This repository contains the complete set of materials for the Monte Carlo Analysis project, conducted in October 2025. The project simulates the carbon footprint of household showers using stochastic modeling to account for variability in user behavior and infrastructure.

## Project Structure

- **Report/**: Contains the final academic project report (`Monte_Carlo_Analysis_Report.pdf`).
- **Presentation/**: Contains the project presentation deck (`Monte_Carlo_Analysis_Presentation.pptx`).
- **Notebook/**: Contains the original Google Colab Jupyter Notebook (`Monte_Carlo_Simulation_Notebook.ipynb`).
- **Data/**:
  - `Input_Data_Template.xlsx`: The original data source and template.
  - `Simulation_Results.xlsx`: Raw results from the 1,000 Monte Carlo simulation runs.
  - Visualizations: Key charts including the CO2 distribution histogram and comparative bar charts.
- **Monte_Carlo_Simulation_Script.py**: The standalone Python script used to execute the simulations.

## Key Insights

1. **Median Carbon Impact**: ~0.87 kg CO₂eq per shower.
2. **Bimodal Distribution**: Identified two distinct emission profiles based on regional energy infrastructure.
3. **Primary Driver**: Regional energy grids (especially low-carbon Swiss heating) significantly impact the total carbon risk.

## Methodology

The study employs a Monte Carlo approach (N=1000) with the following distributions:
- **Shower Length**: Lognormal Distribution (Right-skewed).
- **Heat & Water Source**: Multinomial Distribution based on regional mix.

---
*Completed in October 2025 as part of an environmental impact assessment study.*
