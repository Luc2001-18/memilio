#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Kilian Volmer
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
import numpy as np
import pyabc
import os
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# OBSERVED DATA (2010-2018, 9 yearly values per region)
# =============================================================================

observed_data_dict = {
    "prev_north":  np.array([54.71, 50.02, 45.59, 42.68, 42.36, 44.16, 46.54, 44.84, 42.41]),
    "inc_north":   np.array([532.99, 510.55, 483.01, 454.34, 445.97, 460.37, 476.18, 464.68, 445.21]),
    "prev_center": np.array([42.69, 32.84, 30.12, 32.97, 36.42, 40.64, 44.79, 43.27, 39.82]),
    "inc_center":  np.array([456.79, 377.24, 354.06, 374.01, 395.41, 430.04, 463.12, 452.02, 422.72]),
    "prev_south":  np.array([26.03, 23.15, 22.29, 23.74, 25.56, 28.76, 32.76, 33.62, 33.94]),
    "inc_south":   np.array([318.27, 293.35, 284.87, 294.09, 305.23, 330.18, 362.03, 371.36, 375.60])
}

# =============================================================================
# SUMMARY STATISTICS
# Slices 2000-2018 simulation, keeps only 2010-2018 (burn-in = first 10 years)
# Computes prevalence PfPR (age 2-10) and incidence (total population)
# =============================================================================

def calculate_summary_stats(results_list):
    start_day = 3650  # Day 0 = year 2000; Day 3650 = year 2010
    stats = {}
    region_names = ["north", "center", "south"]

    # Compartment indices in the flattened state vector
    # Group 0 (kids <2):   S=0,  E=1,  IA=2,  IS=3,  R=4,  Sv=5,  Iv=6
    # Group 1 (2-10):      S=7,  E=8,  IA=9,  IS=10, R=11, Sv=12, Iv=13
    # Group 2 (adults):    S=14, E=15, IA=16, IS=17, R=18, Sv=19, Iv=20
    # Group 3 (mosquitos): S=21, E=22, IA=23, IS=24, R=25, Sv=26, Iv=27

    exposed_indices = [1, 8, 15]
    human_indices   = [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18]

    for i, name in enumerate(region_names):
        # results_list[i] is already a numpy array (shape: compartments x timepoints)
        data = results_list[i].as_ndarray().T[start_day:, :] # shape: (timepoints, compartments)

        prev_years = []
        inc_years  = []

        for year in range(9):
            y_data = data[year * 365 : (year + 1) * 365, :]

            # Prevalence PfPR 2-10: (IA + IS) / total_2_10
            pop_2_10 = y_data[:, 7:12].sum(axis=1)
            inf_2_10 = y_data[:, 9] + y_data[:, 10]
            prev_years.append(
                np.mean(inf_2_10 / np.where(pop_2_10 == 0, 1, pop_2_10)) * 100
            )

            # Incidence: new cases per 1000 people per year
            # Approximated as daily flow E -> I = E / TimeExposed, summed over year
            new_cases = np.sum(y_data[:, exposed_indices].sum(axis=1) / 15.0)
            pop_total = y_data[:, human_indices].sum(axis=1).mean()
            inc_years.append((new_cases / pop_total) * 1000)

        stats[f"prev_{name}"] = np.array(prev_years)
        stats[f"inc_{name}"]  = np.array(inc_years)

    return stats

# =============================================================================
# SIMULATION FUNCTION
# =============================================================================
def run_benin_simulation(params):
    """
    Runs the Benin Malaria model simulation (Stochastic Metapopulation Model) via Python bindings.
    """
  #  print("Starting Simulation...")
    print(f"  Running simulation with BitingRates: "
          f"N={params['BitingRateNorth']:.3f}, "
          f"C={params['BitingRateCenter']:.3f}, "
          f"S={params['BitingRateSouth']:.3f}")

    results = oseirvector.simulate(
        t0                           = 0.0,
        tmax                         = 6935.0,
        dt                           = 1.0,
        TimeExposed                  = 15.0,
        TimeInfectedAsymptomatic     = 80.0,
        TimeInfectedSymptomatic      = 7.0,
        TransmissionProbabilityOnContact = 0.1,
        AsymptomaticProbability      = 0.287,
        TimeWaningImmunity           = 180, #730.0,
        BitingRateNorth              = params["BitingRateNorth"],
        BitingRateCenter             = params["BitingRateCenter"],
        BitingRateSouth              = params["BitingRateSouth"]
    )

    return calculate_summary_stats(results)

# =============================================================================
#  DISTANCE FUNCTIONS (one per observable, sum of absolute yearly differences)
# Aggregated and auto-scaled by AdaptiveAggregatedDistance
# =============================================================================
def distance_prev_north(simulation, real_data):
    # Mean absolute error per year, normalized to 0-1
    return np.mean(np.abs(simulation["prev_north"] - real_data["prev_north"])) / 100.0

def distance_inc_north(simulation, real_data):
    return np.mean(np.abs(simulation["inc_north"] - real_data["inc_north"])) / 1000.0

def distance_prev_center(simulation, real_data):
    return np.mean(np.abs(simulation["prev_center"] - real_data["prev_center"])) / 100.0

def distance_inc_center(simulation, real_data):
    return np.mean(np.abs(simulation["inc_center"] - real_data["inc_center"])) / 1000.0

def distance_prev_south(simulation, real_data):
    return np.mean(np.abs(simulation["prev_south"] - real_data["prev_south"])) / 100.0

def distance_inc_south(simulation, real_data):
    return np.mean(np.abs(simulation["inc_south"] - real_data["inc_south"])) / 1000.0

# Combine all 6 distances — auto-scaling handles the prevalence vs incidence
# scale mismatch so they contribute more equally to the total distance
distance = pyabc.AdaptiveAggregatedDistance(
    [
        distance_prev_north,
        distance_inc_north,
        distance_prev_center,
        distance_inc_center,
        distance_prev_south,
        distance_inc_south,
    ],
    adaptive=True,
    scale_function=pyabc.distance.median
)

# =============================================================================
# VISUALIZATION HELPERS
# =============================================================================

def plot_region(history, region_name, observed_data_dict, output_dir="."):
    """
    Plots simulated trajectories vs observed data for one region.
    Saves as <region_name>.png
    """
    prev_key = f"prev_{region_name}"
    inc_key  = f"inc_{region_name}"
    years    = np.arange(2010, 2019)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Benin Malaria — {region_name.capitalize()} Region", fontsize=14)

    for ax, key, ylabel, title in zip(
        axes,
        [prev_key, inc_key],
        ["Prevalence PfPR 2-10 (%)", "Incidence (per 1000/year)"],
        ["Prevalence", "Incidence"]
    ):
    
        # Plot all accepted simulations from last population (grey)
        sims, weights = history.get_distribution(m=0, t=history.max_t)
        for _, row in sims.iterrows():
            sim_result = run_benin_simulation(dict(row))
            print(f"  {key}: {sim_result[key]}")  # ADD THIS
          #  ax.plot(years, sim_result[key], color="grey", alpha=0.1)
            ax.scatter(years, sim_result[key], color="grey", alpha=0.2, s=20)

        # Plot observed data
        ax.scatter(years, observed_data_dict[key],
                   color="C1", zorder=3, label="Observed data")

        # Plot weighted mean
        all_sim_data = []
        w_array = []
        for (_, row), w in zip(sims.iterrows(), weights):
            sim_result = run_benin_simulation(dict(row))
            all_sim_data.append(sim_result[key])
            w_array.append(w)
        w_array  = np.array(w_array)
        w_array /= w_array.sum()
        mean     = (np.array(all_sim_data) * w_array[:, None]).sum(axis=0)
        #ax.plot(years, mean, color="C2", label="Simulation mean", linewidth=2)
        ax.scatter(years, mean, color="C2", zorder=3, label="Simulation mean", s=20, marker="D")
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()

    plt.tight_layout()
    save_path = os.path.join(output_dir, f"{region_name}.png")
    plt.savefig(save_path)
    print(f"  Figure saved: {save_path}")
    plt.close()

# =============================================================================
# MAIN — PRIOR, ABC SETUP, RUN, RESULTS
# =============================================================================

if __name__ == "__main__":

    # --- Prior distributions ---
    prior = pyabc.Distribution(
    BitingRateNorth  = pyabc.RV("uniform", 0.3, 0.5),
    BitingRateCenter = pyabc.RV("uniform", 0.3, 0.5),
    BitingRateSouth  = pyabc.RV("uniform", 0.3, 0.5),
    )

    # --- Quick sanity check before launching ABC ---
    print("=== SANITY CHECK: running one simulation before ABC ===")
    test_params = prior.rvs()
    print(f"  Sampled params: {test_params}")
    test_result = run_benin_simulation(test_params)
    print(f"  prev_north (first year): {test_result['prev_north'][0]:.2f} %")
    print(f"  inc_north  (first year): {test_result['inc_north'][0]:.2f} per 1000")
    print("=== Sanity check passed ===\n")

    # --- ABC setup ---
    population_size = 100  # increase later for better accuracy

    abc = pyabc.ABCSMC(
        run_benin_simulation,
        prior,
        distance,
        population_size=population_size
    )

    db_path = "sqlite:///" + os.path.join(os.getcwd(), "benin_malaria.db")
    abc.new(db_path, observed_data_dict)

    # --- Run calibration ---
    print("=== STARTING CALIBRATION ===")
    history = abc.run(minimum_epsilon=0.5, max_nr_populations=30)
    print(f"\nCalibration finished. Results saved in {db_path}")

    # --- Posterior summary ---
    df, weights = history.get_distribution(m=0, t=history.max_t)
    print("\n=== ESTIMATED BITING RATES (posterior summary) ===")
    print(df[["BitingRateNorth", "BitingRateCenter", "BitingRateSouth"]].describe())

    # --- Plots: one figure per region ---
    print("\n=== GENERATING FIGURES ===")
    for region in ["north", "center", "south"]:
        print(f"  Plotting {region}...")
        plot_region(history, region, observed_data_dict, output_dir=".")

    print("\nDone. Open north.png, center.png, south.png in VS Code to view results.")