#############################################################################
# Malaria IC Calibration
# Uses pyABC to estimate a single scaling factor alpha that modifies
# initial conditions so that after 1 year, the model matches 2010 data
# Biting rates fixed at 0.3 for all regions
#############################################################################

import numpy as np
import pyabc
import os
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# 1. OBSERVED DATA — 2010 targets (one value per region)
# =============================================================================

# observed_2010_2012 = {
#     "prev_north":  np.array([54.71]),
#     "inc_north":   np.array([532.99]),
#     "prev_center": np.array([42.69]),
#     "inc_center":  np.array([456.79]),
#     "prev_south":  np.array([26.03]),
#     "inc_south":   np.array([318.27]),
#  }
observed_2010 = {
    "prev_north":  np.array([54.71]),
    "inc_north":   np.array([532.99]),
    "prev_center": np.array([42.69]),
    "inc_center":  np.array([456.79]),
    "prev_south":  np.array([26.03]),
    "inc_south":   np.array([318.27]),
}


# =============================================================================
# 2. FIXED PARAMETERS
# =============================================================================

BITING_RATE          = 0.4   # same for all regions
T0                   = 0.0
TMAX                 = 365.0
DT                   = 1.0
TIME_EXPOSED         = 15.0
TIME_INFECTED_ASYMPTOMATIC = 80.0
TIME_INFECTED_SYMPTOMATIC  = 7.0
TRANSMISSION_PROB_ON_CONTACT = 0.1
ASYMPTOMATIC_PROBABILITY     = 0.287
TIME_WANING_IMMUNITY         = 180 #1825.0
#SYMPTOMATIC_FRACTION         = 1.0 - ASYMPTOMATIC_PROBABILITY

# =============================================================================
# 3. COMPARTMENT INDICES
# Group 0 (kids <2):   S=0,  E=1,  IA=2,  IS=3,  R=4
# Group 1 (2-10):      S=7,  E=8,  IA=9,  IS=10, R=11
# Group 2 (adults):    S=14, E=15, IA=16, IS=17, R=18
# =============================================================================

COMPARTMENTS = {
    0: {"S": 1,  "E": 2,  "IA": 3,  "IS": 4,  "R": 5},
    1: {"S": 8,  "E": 9,  "IA": 10, "IS": 11, "R": 12},
    2: {"S": 15, "E": 16, "IA": 17, "IS": 18, "R": 19},
}
EXPOSED_INDICES = [2, 9, 16]
HUMAN_INDICES   = [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19]
# =============================================================================
# 4. SUMMARY STATISTICS
# Computes prevalence and incidence after 1 year for each region
# =============================================================================

def calculate_summary_stats(results_list, reporting_rate):
    stats = {}
    region_names = ["north", "center", "south"]

    for i, name in enumerate(region_names):
        data = results_list[i].as_ndarray().T  # shape: (365, 28)

        # Prevalence PfPR 2-10: mean over year of (IA + IS) / N for group 1
        idx1 = COMPARTMENTS[1]
        N_g1 = (data[:, idx1["S"]] + data[:, idx1["E"]] +
                data[:, idx1["IA"]] + data[:, idx1["IS"]] +
                data[:, idx1["R"]])
        inf_g1     = data[:, idx1["IA"]] + data[:, idx1["IS"]]
        prevalence = np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100

        # Age-group specific symptomatic fractions (matches C++ AsymptomaticProbability values)
        # group 0: p_asymp=0.0  → symptomatic fraction = 1.0
        # group 1: p_asymp=0.2  → symptomatic fraction = 0.8
        # group 2: p_asymp=0.35 → symptomatic fraction = 0.65
        SYMPTOMATIC_FRACTIONS = [1.0, 0.8, 0.65]
        E_INDICES = [2, 9, 16]
        idx1 = COMPARTMENTS[1]

        prev_years = []
        inc_years  = []

        for yr in range(1):
            start  = yr * 365
            end    = start + 365
            y_data = data[start:end, :]

            # Prevalence PfPR 2-10
            N_g1 = (y_data[:, idx1["S"]] + y_data[:, idx1["E"]] +
                    y_data[:, idx1["IA"]] + y_data[:, idx1["IS"]] +
                    y_data[:, idx1["R"]])
            inf_g1 = y_data[:, idx1["IA"]] + y_data[:, idx1["IS"]]
            prev_years.append(np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100)

            # Incidence
            daily_clinical = sum(
                SYMPTOMATIC_FRACTIONS[g] * y_data[:, E_INDICES[g]] / TIME_EXPOSED
                for g in range(1)
            )
            pop_at_risk = y_data[0, HUMAN_INDICES].sum()
            inc_years.append((np.sum(daily_clinical) * reporting_rate / pop_at_risk) * 1000)

        # Return as 3-element arrays matching observed_2010_2012 format
        stats[f"prev_{name}"] = np.array(prev_years)
        stats[f"inc_{name}"]  = np.array(inc_years)

    return stats

# =============================================================================
# 5. SIMULATION FUNCTION
# =============================================================================

def run_simulation(params):
    alpha = params["alpha"]
    r = params["reporting_rate"]
    print(f"  Running with alpha={alpha:.4f}")

    results = oseirvector.simulate(
        t0                           = T0,
        tmax                         = TMAX,
        dt                           = DT,
        TimeExposed                  = TIME_EXPOSED,
        TimeInfectedAsymptomatic     = TIME_INFECTED_ASYMPTOMATIC,
        TimeInfectedSymptomatic      = TIME_INFECTED_SYMPTOMATIC,
        TransmissionProbabilityOnContact = TRANSMISSION_PROB_ON_CONTACT,
        AsymptomaticProbability      = ASYMPTOMATIC_PROBABILITY,
        TimeWaningImmunity           = TIME_WANING_IMMUNITY,
        BitingRateNorth              = BITING_RATE,
        BitingRateCenter             = BITING_RATE,
        BitingRateSouth              = BITING_RATE,
        TransmissionVectorToHuman    =0.27, #0.24,
        TransmissionHumanToVector    =0.02, #0.02,
        ic_scale                     = alpha,
    )

    return calculate_summary_stats(results,r )

# =============================================================================
# 6. DISTANCE FUNCTIONS — normalized to 0-1
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
# 7. SANITY CHECK
# =============================================================================

test = run_simulation({"alpha": 1.0, "reporting_rate": 1.0})
print(f"  alpha=1.0 → prev_north year1={test['prev_north'][0]:.2f}%")
print(f"  Target    → 54.71%")
print(f"  alpha=1.0 → inc_north year1={test['inc_north'][0]:.2f}")
print(f"  Target    → 532.99")
# =============================================================================
# 8. MAIN — ABC SETUP AND RUN
# =============================================================================

if __name__ == "__main__":

    prior = pyabc.Distribution(
        alpha = pyabc.RV("uniform", 0.05, 1),  # range: 0.1 to 1.0
        reporting_rate = pyabc.RV("uniform", 0.01, 1)  # Range: 0.2 to 0.8
    )

    population_size = 500  # small for first test, increase later

    abc = pyabc.ABCSMC(
        run_simulation,
        prior,
        distance,
        population_size=population_size
    )

    db_path = "sqlite:///" + os.path.join(os.getcwd(), "ic_calibration.db")
    abc.new(db_path, observed_2010)

    print("=== STARTING IC CALIBRATION ===")
    history = abc.run(minimum_epsilon=0.01, max_nr_populations=10)
    print("\nCalibration finished.")

    # --- Posterior summary ---
    df, weights = history.get_distribution(m=0, t=history.max_t)
    
    # Calculate weighted means
    best_alpha = (df["alpha"] * weights).sum()
    best_r = (df["reporting_rate"] * weights).sum()

    print("\n=== ESTIMATED PARAMETERS (Weighted Means) ===")
    print(f"Alpha:          {best_alpha:.4f}")
    print(f"Reporting Rate: {best_r:.4f}")

    # --- What do the initial conditions look like at best alpha? ---
    print("\n=== INITIAL CONDITIONS AT BEST ALPHA ===")
    prop_E  = [0.01,  0.01,  0.005]
    prop_IA = [0.02,  0.05,  0.05]
    prop_IS = [0.01,  0.01,  0.005]
    prop_R  = [0.10,  0.20,  0.30]
    groups  = ["Kids <2y", "Children 2-10y", "Adults"]
    for i, g in enumerate(groups):

        s = best_alpha * (prop_E[i] + prop_IA[i] + prop_IS[i] + prop_R[i])

        print(f"  {g}: E={best_alpha*prop_E[i]:.4f}  "

              f"IA={best_alpha*prop_IA[i]:.4f}  "

              f"IS={best_alpha*prop_IS[i]:.4f}  "

              f"R={best_alpha*prop_R[i]:.4f}  "

              f"S={1-s:.4f}")

    # =============================================================================
    # 9. COMPREHENSIVE PLOTTING
    # =============================================================================
    
    # Plot 1: Alpha Posterior
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["alpha"], weights=weights, bins=20, color="steelblue", edgecolor="white")
    ax.axvline(best_alpha, color="red", linewidth=2, label=f"Mean = {best_alpha:.4f}")
    ax.set_xlabel("Alpha (IC Scale)")
    ax.set_ylabel("Weighted Count")
    ax.set_title("Posterior Distribution of Alpha")
    ax.legend()
    plt.tight_layout()
    plt.savefig("alpha_posterior.png", dpi=150)
    plt.close()

    # Plot 2: Reporting Rate Posterior
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["reporting_rate"], weights=weights, bins=20, color="seagreen", edgecolor="white")
    ax.axvline(best_r, color="red", linewidth=2, label=f"Mean = {best_r:.4f}")
    ax.set_xlabel("Reporting Rate (r)")
    ax.set_ylabel("Weighted Count")
    ax.set_title("Posterior Distribution of Reporting Rate")
    ax.legend()
    plt.tight_layout()
    plt.savefig("reporting_rate_posterior.png", dpi=150)
    plt.close()

    # Plot 3: Parameter Correlation (2D)
    # This is crucial to see if one parameter is "fighting" the other
    fig, ax = plt.subplots(figsize=(7, 6))
    scatter = ax.scatter(df["alpha"], df["reporting_rate"], c=weights, cmap="viridis", alpha=0.6)
    ax.set_xlabel("Alpha")
    ax.set_ylabel("Reporting Rate")
    ax.set_title("Parameter Correlation: Alpha vs Reporting Rate")
    plt.colorbar(scatter, label='Weight')
    plt.tight_layout()
    plt.savefig("correlation_plot.png", dpi=150)
    plt.close()

    print("\nAll plots saved: alpha_posterior.png, reporting_rate_posterior.png, correlation_plot.png")