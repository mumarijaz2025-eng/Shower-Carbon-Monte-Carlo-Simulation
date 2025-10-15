import numpy as np
import pandas as pd
from scipy.stats import lognorm
import matplotlib.pyplot as plt

# --- Parameters from the Excel document ---

# Fixed parameters
FLOW_RATE_GALLONS_PER_MINUTE = 2.5
HEAT_ENERGY_PER_LITER_MJ = 0.144

# Conversion factor: gallons to liters (1 gallon = 3.78541 liters)
GALLONS_TO_LITERS = 3.78541

# Shower length (Lognormal distribution: mean=10, sigma=3 of the lognormal distribution itself)
SHOWER_LENGTH_MEAN_LOGNORM = 10  # minutes
SHOWER_LENGTH_STD_LOGNORM = 3

# Heat Energy Emission Factor (Multinomial distribution)
heat_emission_factors = {
    "Heat, central or small-scale, natural gas {CH}": {"proportion": 0.3855365835, "CO2": 0.0740342955991673},
    "Heat, central or small-scale, other than natural gas {CH}": {"proportion": 0.2101522058, "CO2": 0.104351332558709},
    "Heat, district or industrial, other than natural gas {CH}": {"proportion": 0.1220737078, "CO2": 0.000259857240833919},
    "Heat, borehole heat pump {CH}": {"proportion": 0.04620258055, "CO2": 0.0293168163797332},
    "Heat, central or small-scale, other than natural gas {CH}_wood": {"proportion": 0.007339874836, "CO2": 0.00601753582460092},
    "Heat, central or small-scale, other than natural gas {CH}_solar": {"proportion": 0.03708568338, "CO2": 0.00279358531172397},
    "Electricity, low voltage, Zurich, Switzerland": {"proportion": 0.1916093641, "CO2": 0.0332840451827167}
}

heat_types = list(heat_emission_factors.keys())
heat_proportions = [data["proportion"] for data in heat_emission_factors.values()]
heat_co2_factors = [data["CO2"] for data in heat_emission_factors.values()]

# Water Emission Factor (Multinomial distribution)
water_emission_factors = {
    "tap water production, underground water with disinfection": {"proportion": 0.1428571429, "CO2": 0.000136390779948015},
    "tap water production, underground water with chemical treatment": {"proportion": 0.1428571429, "CO2": 0.000297428660031339},
    "tap water production, underground water without treatment": {"proportion": 0.1428571429, "CO2": 0.000014},
    "tap water production, conventional treatment": {"proportion": 0.1428571429, "CO2": 0.000218},
    "tap water production, microstrainer treatment": {"proportion": 0.1428571429, "CO2": 0.000376},
    "tap water production, conventional with biological treatment": {"proportion": 0.1428571429, "CO2": 0.000288},
    "tap water production, seawater reverse osmosis, conventional pretreatment, baseline module, single stage": {"proportion": 0.1428571429, "CO2": 0.00577}
}

water_types = list(water_emission_factors.keys())
water_proportions = [data["proportion"] for data in water_emission_factors.values()]
water_co2_factors = [data["CO2"] for data in water_emission_factors.values()]

# Number of simulations
NUM_SIMULATIONS = 1000

# --- Monte Carlo Simulation ---

results = []

# Calculate parameters for scipy.stats.lognorm from mean and std dev of the lognormal distribution
m = SHOWER_LENGTH_MEAN_LOGNORM
sd = SHOWER_LENGTH_STD_LOGNORM

# sigma_normal (s_param) is the standard deviation of the underlying normal distribution
s_param = np.sqrt(np.log((sd**2 / m**2) + 1))
# mu (mu_param) is the mean of the underlying normal distribution
mu_param = np.log(m) - s_param**2 / 2
# scale parameter for scipy.stats.lognorm is exp(mu)
scale_param = np.exp(mu_param)

for i in range(NUM_SIMULATIONS):
    # 1. Sample Shower Length (Lognormal)
    shower_length = lognorm.rvs(s=s_param, scale=scale_param, size=1)[0]

    # 2. Calculate Total Water
    total_water_gallons = shower_length * FLOW_RATE_GALLONS_PER_MINUTE
    total_water_liters = total_water_gallons * GALLONS_TO_LITERS
    total_water_kg = total_water_liters # Assuming 1 liter of water is 1 kg

    # 3. Sample Water Type and its Emission Factor
    water_choice_index = np.random.choice(len(water_types), p=water_proportions)
    water_type = water_types[water_choice_index]
    ef_water = water_co2_factors[water_choice_index]

    # 4. Sample Heat Type and its Emission Factor
    heat_choice_index = np.random.choice(len(heat_types), p=heat_proportions)
    heat_type = heat_types[heat_choice_index]
    ef_heat_energy = heat_co2_factors[heat_choice_index]

    # 5. Calculate Total Heat Energy (in MJ)
    total_heat_mj = total_water_liters * HEAT_ENERGY_PER_LITER_MJ

    # 6. Calculate Total CO2eq
    emissions_from_water = total_water_kg * ef_water
    emissions_from_heating = total_heat_mj * ef_heat_energy

    total_co2eq = emissions_from_water + emissions_from_heating

    results.append({
        "simulation_number": i + 1,
        "shower_length": shower_length,
        "total_water": total_water_liters,
        "water_type": water_type,
        "heat_type": heat_type,
        "total_heat": total_heat_mj,
        "total_CO2eq": total_co2eq
    })

df_results = pd.DataFrame(results)

# --- Analysis and Reporting ---

# 1. Report on the lowest, highest and median shower in terms of carbon impact, as well as percentiles (Interquartile range)
lowest_carbon = df_results["total_CO2eq"].min()
highest_carbon = df_results["total_CO2eq"].max()
median_carbon = df_results["total_CO2eq"].median()
q1_carbon = df_results["total_CO2eq"].quantile(0.25)
q3_carbon = df_results["total_CO2eq"].quantile(0.75)

print("\n--- Monte Carlo Simulation Results ---")
print(f"Number of simulations: {NUM_SIMULATIONS}")
print("\nCarbon Impact Statistics (kg CO2eq):")
print(f"  Lowest: {lowest_carbon:.4f}")
print(f"  Highest: {highest_carbon:.4f}")
print(f"  Median: {median_carbon:.4f}")
print(f"  25th Percentile (Q1): {q1_carbon:.4f}")
print(f"  75th Percentile (Q3): {q3_carbon:.4f}")
print(f"  Interquartile Range (IQR): {(q3_carbon - q1_carbon):.4f}")

# Save results to CSV for further analysis/reporting
df_results.to_csv("/home/ubuntu/monte_carlo_results.csv", index=False)
print("\nSimulation results saved to /home/ubuntu/monte_carlo_results.csv")

# Also save the raw data for the report
df_results.to_excel("/home/ubuntu/monte_carlo_results.xlsx", index=False)
print("Simulation results saved to /home/ubuntu/monte_carlo_results.xlsx")

# Basic visualization (histogram of CO2eq)
plt.figure(figsize=(10, 6))
plt.hist(df_results["total_CO2eq"], bins=50, edgecolor='black')
plt.title('Distribution of Shower Carbon Emissions (Monte Carlo Simulation)')
plt.xlabel('Total CO2eq (kg)')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.savefig("/home/ubuntu/co2_distribution_histogram.png")
print("Histogram saved to /home/ubuntu/co2_distribution_histogram.png")

# Further analysis: breakdown by water and heat type
# Average CO2eq by water type
water_type_avg_co2 = df_results.groupby('water_type')['total_CO2eq'].mean().sort_values(ascending=False)
print("\nAverage CO2eq by Water Type:")
print(water_type_avg_co2.to_string())
plt.figure(figsize=(10, 6))
water_type_avg_co2.plot(kind='barh')
plt.title('Average CO2eq by Water Type')
plt.xlabel('Average CO2eq (kg)')
plt.ylabel('Water Type')
plt.tight_layout()
plt.savefig("/home/ubuntu/co2_by_water_type.png")
print("Bar chart for CO2 by water type saved to /home/ubuntu/co2_by_water_type.png")

# Average CO2eq by heat type
heat_type_avg_co2 = df_results.groupby('heat_type')['total_CO2eq'].mean().sort_values(ascending=False)
print("\nAverage CO2eq by Heat Type:")
print(heat_type_avg_co2.to_string())
plt.figure(figsize=(10, 6))
heat_type_avg_co2.plot(kind='barh')
plt.title('Average CO2eq by Heat Type')
plt.xlabel('Average CO2eq (kg)')
plt.ylabel('Heat Type')
plt.tight_layout()
plt.savefig("/home/ubuntu/co2_by_heat_type.png")
print("Bar chart for CO2 by heat type saved to /home/ubuntu/co2_by_heat_type.png")
