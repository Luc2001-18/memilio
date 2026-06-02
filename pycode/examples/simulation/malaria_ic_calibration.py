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

BITING_RATE_NORTH          = 0.4
BITING_RATE_CENTER         = 0.4
BITING_RATE_SOUTH          = 0.4
T0                   = 0.0
TMAX                 = 365.0
DT                   = 1.0
TIME_EXPOSED         = 15.0
TIME_INFECTED_ASYMPTOMATIC = 100
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
        
        # Age-group specific symptomatic fractions (matches C++ AsymptomaticProbability values)
        # group 0: p_asymp=0.0  → symptomatic fraction = 1.0
        # group 1: p_asymp=0.2  → symptomatic fraction = 0.8
        # group 2: p_asymp=0.35 → symptomatic fraction = 0.65
        SYMPTOMATIC_FRACTIONS = [0.95, 0.7, 0.3]
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
                for g in range(3)
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
    r = params["reporting_rate"]
    print(f"  Running with prop_E=({params['prop_E_north']:.4f}, "
          f"{params['prop_E_center']:.4f}, {params['prop_E_south']:.4f}) "
          f"r={r:.3f}")

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
        BitingRateNorth              = BITING_RATE_NORTH,
        BitingRateCenter             = BITING_RATE_CENTER,
        BitingRateSouth              = BITING_RATE_SOUTH,
        TransmissionVectorToHuman    = 0.27,
        TransmissionHumanToVector    = 0.02,
        prop_E_north                 = params["prop_E_north"],
        prop_E_center                = params["prop_E_center"],
        prop_E_south                 = params["prop_E_south"],
    )

    return calculate_summary_stats(results, r)

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

test = run_simulation({"prop_E_north": 0.10, "prop_E_center": 0.08,
                       "prop_E_south": 0.06, "reporting_rate": 1.0})
print(f"  prop_E=0.10 → prev_north year1={test['prev_north'][0]:.2f}%")
print(f"  Target      → 54.71%")
print(f"  prop_E=0.10 → inc_north  year1={test['inc_north'][0]:.2f}")
print(f"  Target      → 532.99")
# =============================================================================
# 8. MAIN — ABC SETUP AND RUN
# =============================================================================

if __name__ == "__main__":

    prior = pyabc.Distribution(
        prop_E_north   = pyabc.RV("uniform", 0.01, 0.29),
        prop_E_center  = pyabc.RV("uniform", 0.01, 0.28),
        prop_E_south   = pyabc.RV("uniform", 0.01, 0.28),
        reporting_rate = pyabc.RV("uniform", 0.01, 0.99)
    )
    population_size = 1000  # small for first test, increase later

    abc = pyabc.ABCSMC(
        run_simulation,
        prior,
        distance,
        population_size=population_size
    )

    db_path = "sqlite:///" + os.path.join(os.getcwd(), "ic_calibration.db")
    abc.new(db_path, observed_2010)

    print("=== STARTING IC CALIBRATION ===")
    history = abc.run(minimum_epsilon=0.005, max_nr_populations=10)
    print("\nCalibration finished.")

    # --- Posterior summary ---
    df, weights = history.get_distribution(m=0, t=history.max_t)
    
    # Calculate weighted means
    best_prop_E_north  = (df["prop_E_north"]  * weights).sum()
    best_prop_E_center = (df["prop_E_center"] * weights).sum()
    best_prop_E_south  = (df["prop_E_south"]  * weights).sum()
    best_r             = (df["reporting_rate"] * weights).sum()

    print("\n=== ESTIMATED PARAMETERS (Weighted Means) ===")
    print(f"prop_E North:   {best_prop_E_north:.4f}")
    print(f"prop_E Center:  {best_prop_E_center:.4f}")
    print(f"prop_E South:   {best_prop_E_south:.4f}")
    print(f"Reporting Rate: {best_r:.4f}")

    # Show implied PfPR from estimated prop_E
    p_asymp_g1 = 0.3
    T_E = 15.0; T_IA = 100.0; T_IS = 7.0
    pfpr_multiplier = p_asymp_g1*(T_IA/T_E) + (1-p_asymp_g1)*(T_IS/T_E)
    print(f"\n=== IMPLIED PfPR 2-10 FROM prop_E ===")
    print(f"  North:  {best_prop_E_north  * pfpr_multiplier * 100:.2f}%  (target 54.71%)")
    print(f"  Center: {best_prop_E_center * pfpr_multiplier * 100:.2f}%  (target 42.69%)")
    print(f"  South:  {best_prop_E_south  * pfpr_multiplier * 100:.2f}%  (target 26.03%)")

   
    # =============================================================================
    # 9. COMPREHENSIVE PLOTTING
    # =============================================================================

    # Figure 1 — Posterior distributions of the three alphas
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("Posterior Distributions of Alpha Parameters", fontsize=14)

    for ax, param, best_val, color, label in zip(
        axes,
        ["prop_E_north", "prop_E_center", "prop_E_south"],
        [best_prop_E_north, best_prop_E_center, best_prop_E_south],
        ["steelblue", "darkorange", "seagreen"],
        ["prop_E North", "prop_E Center", "prop_E South"]
    ):
        ax.hist(df[param], weights=weights, bins=20,
                color=color, edgecolor="white", alpha=0.8)
        ax.axvline(best_val, color="red", linewidth=2,
                label=f"Mean = {best_val:.4f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Weighted Count")
        ax.set_title(label)
        ax.legend()

    plt.tight_layout()
    plt.savefig("alpha_posteriors.png", dpi=150)
    plt.close()
    print("  Saved: alpha_posteriors.png")

    
    # Figure 2 — Correlation plots: each alpha vs reporting rate
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Correlation: prop_E Parameters vs Reporting Rate", fontsize=14)

    for ax, param, label, color in zip(
        axes,
        ["prop_E_north", "prop_E_center", "prop_E_south"],
        ["prop_E North", "prop_E Center", "prop_E South"],
        ["steelblue", "darkorange", "seagreen"]
    ):
        scatter = ax.scatter(df[param], df["reporting_rate"],
                            c=weights, cmap="viridis", alpha=0.7, s=30)
        ax.set_xlabel(label)
        ax.set_ylabel("Reporting Rate")
        ax.set_title(f"{label} vs Reporting Rate")
        plt.colorbar(scatter, ax=ax, label="Weight")

    plt.tight_layout()
    plt.savefig("prop_E_correlation.png", dpi=150)
    plt.close()
    print("  Saved: prop_E_correlation.png")

    # Figure 3 — Reporting rate posterior
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["reporting_rate"], weights=weights, bins=20,
            color="purple", edgecolor="white", alpha=0.8)
    ax.axvline(best_r, color="red", linewidth=2,
            label=f"Mean = {best_r:.4f}")
    ax.set_xlabel("Reporting Rate")
    ax.set_ylabel("Weighted Count")
    ax.set_title("Posterior Distribution of Reporting Rate")
    ax.legend()
    plt.tight_layout()
    plt.savefig("reporting_rate_posterior.png", dpi=150)
    plt.close()
    print("  Saved: reporting_rate_posterior.png")

    print("\nAll plots saved.")