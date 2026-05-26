#############################################################################
# Malaria Full Calibration
# Estimates all 7 parameters simultaneously:
#   - alpha_north, alpha_center, alpha_south (initial condition scaling)
#   - reporting_rate
#   - BitingRateNorth, BitingRateCenter, BitingRateSouth
# Fits to observed prevalence and incidence 2010-2018 for 3 regions
#############################################################################

import numpy as np
import pyabc
import os
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# 1. OBSERVED DATA (2010-2018, 9 yearly values per region)
# =============================================================================

observed_data = {
    "prev_north":  np.array([54.71, 50.02, 45.59, 42.68, 42.36, 44.16, 46.54, 44.84, 42.41]),
    "inc_north":   np.array([532.99, 510.55, 483.01, 454.34, 445.97, 460.37, 476.18, 464.68, 445.21]),
    "prev_center": np.array([42.69, 32.84, 30.12, 32.97, 36.42, 40.64, 44.79, 43.27, 39.82]),
    "inc_center":  np.array([456.79, 377.24, 354.06, 374.01, 395.41, 430.04, 463.12, 452.02, 422.72]),
    "prev_south":  np.array([26.03, 23.15, 22.29, 23.74, 25.56, 28.76, 32.76, 33.62, 33.94]),
    "inc_south":   np.array([318.27, 293.35, 284.87, 294.09, 305.23, 330.18, 362.03, 371.36, 375.60])
}

# =============================================================================
# 2. FIXED PARAMETERS
# =============================================================================

TMAX                         = 3285.0   # 9 years: 2010-2018
DT                           = 1.0
TIME_EXPOSED                 = 15.0
TIME_INFECTED_ASYMPTOMATIC   = 100.0
TIME_INFECTED_SYMPTOMATIC    = 7.0
TRANSMISSION_PROB_ON_CONTACT = 0.1
ASYMPTOMATIC_PROBABILITY     = 0.287
TIME_WANING_IMMUNITY         = 180.0    # adults; kids overridden in C++
TRANSMISSION_VECTOR_TO_HUMAN = 0.27
TRANSMISSION_HUMAN_TO_VECTOR = 0.02

# =============================================================================
# 3. COMPARTMENT INDICES
# Time is column 0 — all indices shifted by +1
# Group 0 (kids <2):   S=1,  E=2,  IA=3,  IS=4,  R=5
# Group 1 (2-10):      S=8,  E=9,  IA=10, IS=11, R=12
# Group 2 (adults):    S=15, E=16, IA=17, IS=18, R=19
# =============================================================================

HUMAN_INDICES        = [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19]
SYMPTOMATIC_FRACTIONS = [0.95, 0.7, 0.3]   # 1 - AsymptomaticProbability per group
E_INDICES             = [2, 9, 16]

# =============================================================================
# 4. SUMMARY STATISTICS
# =============================================================================

def calculate_summary_stats(results_list, reporting_rate):
    stats = {}
    region_names = ["north", "center", "south"]

    for i, name in enumerate(region_names):
        data = results_list[i].as_ndarray().T  # shape: (timepoints, 29)

        prev_years = []
        inc_years  = []

        for year in range(9):
            y_data = data[year * 365 : (year + 1) * 365, :]

            # Prevalence PfPR 2-10: mean(IA_g1 + IS_g1) / N_g1
            pop_2_10 = y_data[:, 8:13].sum(axis=1)
            inf_2_10 = y_data[:, 10] + y_data[:, 11]
            prev_years.append(
                np.mean(inf_2_10 / np.where(pop_2_10 == 0, 1, pop_2_10)) * 100
            )

            # Incidence: symptomatic new cases per 1000 total population per year
            daily_new   = sum(
                SYMPTOMATIC_FRACTIONS[g] * y_data[:, E_INDICES[g]] / TIME_EXPOSED
                for g in range(3)
            )
            pop_at_risk = y_data[0, HUMAN_INDICES].sum()
            inc_years.append(
                (np.sum(daily_new) * reporting_rate / pop_at_risk) * 1000
            )

        stats[f"prev_{name}"] = np.array(prev_years)
        stats[f"inc_{name}"]  = np.array(inc_years)

    return stats

# =============================================================================
# 5. SIMULATION FUNCTION
# =============================================================================

def run_simulation(params):
    print(f"  BR=({params['BitingRateNorth']:.3f}, "
          f"{params['BitingRateCenter']:.3f}, {params['BitingRateSouth']:.3f}) "
          f"alpha=({params['alpha_north']:.3f}, "
          f"{params['alpha_center']:.3f}, {params['alpha_south']:.3f}) "
          f"r={params['reporting_rate']:.3f}")

    results = oseirvector.simulate(
        t0                           = 0.0,
        tmax                         = TMAX,
        dt                           = DT,
        TimeExposed                  = TIME_EXPOSED,
        TimeInfectedAsymptomatic     = TIME_INFECTED_ASYMPTOMATIC,
        TimeInfectedSymptomatic      = TIME_INFECTED_SYMPTOMATIC,
        TransmissionProbabilityOnContact = TRANSMISSION_PROB_ON_CONTACT,
        AsymptomaticProbability      = ASYMPTOMATIC_PROBABILITY,
        TimeWaningImmunity           = TIME_WANING_IMMUNITY,
        TransmissionVectorToHuman    = TRANSMISSION_VECTOR_TO_HUMAN,
        TransmissionHumanToVector    = TRANSMISSION_HUMAN_TO_VECTOR,
        BitingRateNorth              = params["BitingRateNorth"],
        BitingRateCenter             = params["BitingRateCenter"],
        BitingRateSouth              = params["BitingRateSouth"],
        ic_scale_north               = params["alpha_north"],
        ic_scale_center              = params["alpha_center"],
        ic_scale_south               = params["alpha_south"],
    )

    return calculate_summary_stats(results, params["reporting_rate"])

# =============================================================================
# 6. DISTANCE FUNCTIONS
# =============================================================================

def distance_prev_north(sim, obs):
    return np.mean(np.abs(sim["prev_north"] - obs["prev_north"])) / 100.0

def distance_inc_north(sim, obs):
    return np.mean(np.abs(sim["inc_north"] - obs["inc_north"])) / 1000.0

def distance_prev_center(sim, obs):
    return np.mean(np.abs(sim["prev_center"] - obs["prev_center"])) / 100.0

def distance_inc_center(sim, obs):
    return np.mean(np.abs(sim["inc_center"] - obs["inc_center"])) / 1000.0

def distance_prev_south(sim, obs):
    return np.mean(np.abs(sim["prev_south"] - obs["prev_south"])) / 100.0

def distance_inc_south(sim, obs):
    return np.mean(np.abs(sim["inc_south"] - obs["inc_south"])) / 1000.0

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
# 7. VISUALIZATION
# =============================================================================

def plot_region(history, region_name, output_dir="."):
    prev_key = f"prev_{region_name}"
    inc_key  = f"inc_{region_name}"
    years    = np.arange(2010, 2019)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Benin Malaria — {region_name.capitalize()} Region", fontsize=14)

    sims, weights = history.get_distribution(m=0, t=history.max_t)
    all_sim_data  = {prev_key: [], inc_key: []}
    w_array       = []

    for (_, row), w in zip(sims.iterrows(), weights):
        sim_result = run_simulation(dict(row))
        all_sim_data[prev_key].append(sim_result[prev_key])
        all_sim_data[inc_key].append(sim_result[inc_key])
        w_array.append(w)

    w_array  = np.array(w_array)
    w_array /= w_array.sum()

    for ax, key, ylabel, title in zip(
        axes,
        [prev_key, inc_key],
        ["Prevalence PfPR 2-10 (%)", "Incidence (per 1000/year)"],
        ["Prevalence", "Incidence"]
    ):
       # for sim_vals in all_sim_data[key]:
       #     ax.scatter(years, sim_vals, color="grey", alpha=0.2, s=20)
        ax.scatter(years, observed_data[key],
                   color="C1", zorder=3, label="Observed data")
        mean = (np.array(all_sim_data[key]) * w_array[:, None]).sum(axis=0)
        ax.scatter(years, mean, color="C2", zorder=3,
                   label="Simulation mean", s=20, marker="D")
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()

    plt.tight_layout()
    save_path = os.path.join(output_dir, f"{region_name}_full.png")
    plt.savefig(save_path, dpi=150)
    print(f"  Saved: {save_path}")
    plt.close()


def plot_posteriors(df, weights, best_vals, output_dir="."):
    params_info = [
        ("BitingRateNorth",  best_vals["br_north"],   "steelblue",  "Biting Rate North"),
        ("BitingRateCenter", best_vals["br_center"],  "darkorange", "Biting Rate Center"),
        ("BitingRateSouth",  best_vals["br_south"],   "seagreen",   "Biting Rate South"),
        ("alpha_north",      best_vals["alpha_north"], "steelblue", "Alpha North"),
        ("alpha_center",     best_vals["alpha_center"],"darkorange","Alpha Center"),
        ("alpha_south",      best_vals["alpha_south"], "seagreen",  "Alpha South"),
        ("reporting_rate",   best_vals["r"],           "purple",    "Reporting Rate"),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(20, 8))
    fig.suptitle("Posterior Distributions — All Parameters", fontsize=14)

    for ax, (param, best_val, color, label) in zip(axes.flatten(), params_info):
        ax.hist(df[param], weights=weights, bins=20,
                color=color, edgecolor="white", alpha=0.8)
        ax.axvline(best_val, color="red", linewidth=2,
                   label=f"Mean = {best_val:.4f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Weighted Count")
        ax.set_title(label)
        ax.legend(fontsize=8)

    axes.flatten()[-1].set_visible(False)
    plt.tight_layout()
    save_path = os.path.join(output_dir, "all_parameters_posterior.png")
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"  Saved: {save_path}")

# =============================================================================
# 8. MAIN
# =============================================================================

if __name__ == "__main__":

    # --- Prior ---
    prior = pyabc.Distribution(
        alpha_north      = pyabc.RV("uniform", 0.01, 0.99),
        alpha_center     = pyabc.RV("uniform", 0.01, 0.99),
        alpha_south      = pyabc.RV("uniform", 0.01, 0.99),
        reporting_rate   = pyabc.RV("uniform", 0.01, 0.99),
        BitingRateNorth  = pyabc.RV("uniform", 0.1, 0.8),
        BitingRateCenter = pyabc.RV("uniform", 0.1, 0.8),
        BitingRateSouth  = pyabc.RV("uniform", 0.1, 0.8),
    )

    # --- Sanity check ---
    print("=== SANITY CHECK ===")
    test_params = {
        "alpha_north": 0.5, "alpha_center": 0.5, "alpha_south": 0.5,
        "reporting_rate": 0.45,
        "BitingRateNorth": 0.5, "BitingRateCenter": 0.3, "BitingRateSouth": 0.25,
    }
    test_result = run_simulation(test_params)
    print(f"  prev_north  year1 = {test_result['prev_north'][0]:.2f}%  (target: 54.71%)")
    print(f"  inc_north   year1 = {test_result['inc_north'][0]:.2f}    (target: 532.99)")
    print(f"  prev_south  year1 = {test_result['prev_south'][0]:.2f}%  (target: 26.03%)")
    print(f"  inc_south   year1 = {test_result['inc_south'][0]:.2f}    (target: 318.27)")
    print("=== Sanity check done ===\n")

    # --- ABC setup ---
    population_size = 100

    abc = pyabc.ABCSMC(
        run_simulation,
        prior,
        distance,
        population_size=population_size
    )

    db_path = "sqlite:///" + os.path.join(os.getcwd(), "benin_malaria_full.db")
    abc.new(db_path, observed_data)

    # --- Run ---
    print("=== STARTING FULL CALIBRATION ===")
    history = abc.run(minimum_epsilon=0.3, max_nr_populations=20)
    print(f"\nCalibration finished. Results saved in {db_path}")

    # --- Posterior summary ---
    df, weights = history.get_distribution(m=0, t=history.max_t)

    best_vals = {
        "br_north":    (df["BitingRateNorth"]  * weights).sum(),
        "br_center":   (df["BitingRateCenter"] * weights).sum(),
        "br_south":    (df["BitingRateSouth"]  * weights).sum(),
        "alpha_north":  (df["alpha_north"]  * weights).sum(),
        "alpha_center": (df["alpha_center"] * weights).sum(),
        "alpha_south":  (df["alpha_south"]  * weights).sum(),
        "r":            (df["reporting_rate"] * weights).sum(),
    }

    # --- Plots ---
    print("\n=== GENERATING PLOTS ===")
    plot_posteriors(df, weights, best_vals)

    for region in ["north", "center", "south"]:
        print(f"  Plotting {region}...")
        plot_region(history, region)

    print("\n" + "=" * 55)
    print(f"{'ESTIMATED PARAMETERS (weighted means)':^55}")
    print("=" * 55)
    print(f"  Biting Rate North:   {best_vals['br_north']:.4f}")
    print(f"  Biting Rate Center:  {best_vals['br_center']:.4f}")
    print(f"  Biting Rate South:   {best_vals['br_south']:.4f}")
    print(f"  Alpha North:         {best_vals['alpha_north']:.4f}")
    print(f"  Alpha Center:        {best_vals['alpha_center']:.4f}")
    print(f"  Alpha South:         {best_vals['alpha_south']:.4f}")
    print(f"  Reporting Rate:      {best_vals['r']:.4f}")
    print("=" * 55)

    print("\nDone.")