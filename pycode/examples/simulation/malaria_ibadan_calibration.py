#############################################################################
# Malaria Ibadan Calibration
# Uses pyABC to estimate prop_E and BitingRate for Ibadan, Nigeria
# using 2015-2016 monthly prevalence data (24 data points)
# We use the South region as proxy for Ibadan
# No ITN effect — ITN block commented out in C++
#############################################################################

import numpy as np
import pandas as pd
import pyabc
import os
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# 1. LOAD OBSERVED DATA — 2015-2016 monthly prevalence (24 points)
# =============================================================================

df = pd.read_excel(
    "ibadan_malaria-prevalence_dataset_1996-2017.xlsx",
    sheet_name="vars-preP_val"
)

ibadan_2015_2016 = df[(df["year"] >= 2015) & (df["year"] <= 2016)].copy().reset_index(drop=True)
observed_prev = ibadan_2015_2016["preP"].values * 100  # convert to percentage

print("=== OBSERVED DATA 2015-2016 ===")
for i, v in enumerate(observed_prev):
    print(f"  Month {i+1:2d}: {v:.2f}%")
print(f"  Mean: {observed_prev.mean():.2f}%")

# pyABC observed dictionary
observed_data = {
    "monthly_prev": observed_prev  # shape (24,)
}

# =============================================================================
# 2. FIXED PARAMETERS
# =============================================================================

T0   = 0.0
TMAX = 730.0  # 2 years: 2015-2016
DT   = 1.0

TIME_EXPOSED                 = 15.0
TIME_INFECTED_ASYMPTOMATIC   = 100.0
TIME_INFECTED_SYMPTOMATIC    = 7.0
TRANSMISSION_PROB_ON_CONTACT = 0.1
ASYMPTOMATIC_PROBABILITY     = 0.287
TIME_WANING_IMMUNITY         = 180.0
TRANSMISSION_VECTOR_TO_HUMAN = 0.27
TRANSMISSION_HUMAN_TO_VECTOR = 0.02

# =============================================================================
# 3. COMPARTMENT INDICES
# Time is column 0 — all indices shifted by +1
# Group 1 (2-10): S=8, E=9, IA=10, IS=11, R=12
# =============================================================================

S_g1  = 8
E_g1  = 9
IA_g1 = 10
IS_g1 = 11
R_g1  = 12

# =============================================================================
# 4. SUMMARY STATISTICS
# Compute monthly average prevalence for group 1 (children 2-10)
# over 24 months using South region output
# =============================================================================

def calculate_summary_stats(results_list):
    # Use South region (index 2) as Ibadan proxy
    data = results_list[2].as_ndarray().T  # shape: (731, 29)

    monthly_prev = []

    for month_idx in range(24):
        day_start = month_idx * 30
        day_end   = day_start + 30
        y_data    = data[day_start:day_end, :]

        N_g1   = (y_data[:, S_g1] + y_data[:, E_g1] +
                  y_data[:, IA_g1] + y_data[:, IS_g1] +
                  y_data[:, R_g1])
        inf_g1 = y_data[:, IA_g1] + y_data[:, IS_g1]
        prev   = np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100
        monthly_prev.append(prev)

    return {"monthly_prev": np.array(monthly_prev)}

# =============================================================================
# 5. SIMULATION FUNCTION
# =============================================================================

def run_simulation(params):
    biting_rate = params["biting_rate"]
    prop_E      = params["prop_E"]

    print(f"  biting_rate={biting_rate:.4f}, prop_E={prop_E:.4f}")

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
        TransmissionVectorToHuman    = TRANSMISSION_VECTOR_TO_HUMAN,
        TransmissionHumanToVector    = TRANSMISSION_HUMAN_TO_VECTOR,
        BitingRateNorth              = biting_rate,
        BitingRateCenter             = biting_rate,
        BitingRateSouth              = biting_rate,
        prop_E_north                 = prop_E,
        prop_E_center                = prop_E,
        prop_E_south                 = prop_E,
    )

    return calculate_summary_stats(results)

# =============================================================================
# 6. DISTANCE FUNCTION
# Mean absolute error over 24 monthly prevalence values, normalized
# =============================================================================

def distance(sim, obs):
    return np.mean(np.abs(sim["monthly_prev"] - obs["monthly_prev"])) / 100.0

# =============================================================================
# 7. SANITY CHECK
# =============================================================================

print("\n=== SANITY CHECK ===")
test = run_simulation({"biting_rate": 0.15, "prop_E": 0.06})
print(f"  Model monthly prev: {test['monthly_prev'].round(2)}")
print(f"  Model mean:         {test['monthly_prev'].mean():.2f}%")
print(f"  Observed mean:      {observed_prev.mean():.2f}%")
print("=== Sanity check done ===\n")

# =============================================================================
# 8. MAIN — ABC SETUP AND RUN
# =============================================================================

if __name__ == "__main__":

    prior = pyabc.Distribution(
        biting_rate = pyabc.RV("uniform", 0.05, 0.45),  # range [0.05, 0.50]
        prop_E      = pyabc.RV("uniform", 0.01, 0.14),  # range [0.01, 0.15]
    )

    population_size = 500

    abc = pyabc.ABCSMC(
        run_simulation,
        prior,
        distance,
        population_size=population_size
    )

    db_path = "sqlite:///" + os.path.join(os.getcwd(), "ibadan_calibration.db")
    abc.new(db_path, observed_data)

    print("=== STARTING IBADAN CALIBRATION ===")
    history = abc.run(minimum_epsilon=0.005, max_nr_populations=15)
    print("\nCalibration finished.")

    # --- Posterior summary ---
    df_post, weights = history.get_distribution(m=0, t=history.max_t)

    best_br     = (df_post["biting_rate"] * weights).sum()
    best_prop_E = (df_post["prop_E"]      * weights).sum()

    print("\n" + "=" * 45)
    print(f"{'ESTIMATED PARAMETERS':^45}")
    print("=" * 45)
    print(f"  Biting Rate:  {best_br:.4f}")
    print(f"  prop_E:       {best_prop_E:.4f}")
    print(f"  Implied PfPR: {best_prop_E * 2.327 * 100:.2f}%")
    print("=" * 45)

    # --- Verify with best parameters ---
    best_result = run_simulation({"biting_rate": best_br, "prop_E": best_prop_E})

    month_names_ext = [
        "Jan\n2015","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec\n2015",
        "Jan\n2016","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec\n2016"
    ]

    print(f"\n=== MODEL VS OBSERVED (2015-2016) ===")
    print(f"  {'Month':<12} {'Model (%)':>10} {'Observed (%)':>14}")
    print("  " + "-" * 40)
    for i in range(24):
        print(f"  {month_names_ext[i]:<12} {best_result['monthly_prev'][i]:>10.2f} "
              f"{observed_prev[i]:>14.2f}")
    print(f"\n  Model mean:    {best_result['monthly_prev'].mean():.2f}%")
    print(f"  Observed mean: {observed_prev.mean():.2f}%")

    # --- Posterior plots ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle("Posterior Distributions — Ibadan Calibration (2015-2016)", fontsize=14)

    for ax, param, best_val, color, label in zip(
        axes,
        ["biting_rate", "prop_E"],
        [best_br, best_prop_E],
        ["steelblue", "darkorange"],
        ["Biting Rate", "prop_E"]
    ):
        ax.hist(df_post[param], weights=weights, bins=20,
                color=color, edgecolor="white", alpha=0.8)
        ax.axvline(best_val, color="red", linewidth=2,
                   label=f"Mean = {best_val:.4f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Weighted Count")
        ax.set_title(f"Posterior: {label}")
        ax.legend()

    plt.tight_layout()
    plt.savefig("ibadan_posteriors.png", dpi=150)
    plt.close()
    print("\n  Saved: ibadan_posteriors.png")

    # --- Model vs observed plot ---
    months = np.arange(1, 25)

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(months, best_result["monthly_prev"],
            color="steelblue", linewidth=2, marker="D",
            markersize=6, label="Model (South/Ibadan)")
    ax.scatter(months, observed_prev,
               color="C1", s=80, zorder=3, label="Observed Ibadan 2015-2016")
    ax.set_xlabel("Month")
    ax.set_ylabel("Prevalence PfPR (%)")
    ax.set_title("Model vs Observed Monthly Prevalence — Ibadan 2015-2016")
    ax.set_xticks(months)
    ax.set_xticklabels(month_names_ext, fontsize=7)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("ibadan_fit_2015_2016.png", dpi=150)
    plt.close()
    print("  Saved: ibadan_fit_2015_2016.png")
    print("\nDone.")
