#############################################################################
# Malaria Diagnostics — Year-by-year exploration
# Simulates the model for 1 year with fixed biting rates
# Plots S, IA, IS, R per age group per region
# Prints prevalence (PfPR 2-10) and incidence per region
#############################################################################

import numpy as np
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# CONFIGURATION 
# =============================================================================

BITING_RATE_NORTH  = 0.7
BITING_RATE_CENTER = 0.5
BITING_RATE_SOUTH  = 0.5
REPORTING_RATE = 0.1848
T0   = 0.0
TMAX = 1095.0   # one year
DT   = 1.0     # daily resolution
#alpha = 0.3776
PROP_E_NORTH  = 0.2166
PROP_E_CENTER = 0.2089
PROP_E_SOUTH  = 0.0749
# Fixed parameters (not being calibrated here)
TIME_EXPOSED                    = 15.0
TIME_INFECTED_ASYMPTOMATIC      = 100 #80.0
TIME_INFECTED_SYMPTOMATIC       = 7.0
TIME_WANING_IMMUNITY            = 180.0 #180 #730.0
TRANSMISSION_VECTOR_TO_HUMAN    = 0.27 #0.27
TRANSMISSION_HUMAN_TO_VECTOR    = 0.02

# Time is column 0 — all indices shifted by +1
# Group 0 (kids <2):   S=1,  E=2,  IA=3,  IS=4,  R=5,  Sv=6,  Iv=7
# Group 1 (2-10):      S=8,  E=9,  IA=10, IS=11, R=12, Sv=13, Iv=14
# Group 2 (adults):    S=15, E=16, IA=17, IS=18, R=19, Sv=20, Iv=21
# Group 3 (mosquitos): S=22, E=23, IA=24, IS=25, R=26, Sv=27, Iv=28

COMPARTMENTS = {
    0: {"S": 1,  "E": 2,  "IA": 3,  "IS": 4,  "R": 5},
    1: {"S": 8,  "E": 9,  "IA": 10, "IS": 11, "R": 12},
    2: {"S": 15, "E": 16, "IA": 17, "IS": 18, "R": 19},
}

HUMAN_INDICES   = [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19]
SYMPTOMATIC_FRACTIONS = [0.95, 0.7, 0.3]
E_INDICES             = [2, 9, 16]

AGE_GROUP_NAMES = ["Kids < 2y", "Children 2-10y", "Adults 10+"]
REGION_NAMES    = ["north", "center", "south"]

# OBSERVED = {
#     "north":  {"prev": [54.71, 50.02], "inc": [532.99, 510.55]},
#     "center": {"prev": [42.69, 32.84], "inc": [456.79, 377.24]},
#     "south":  {"prev": [26.03, 23.15], "inc": [318.27, 293.35]},
# }
OBSERVED = {
    "north":  {"prev": [54.71, 50.02, 45.59], "inc": [532.99, 510.55, 483.01]},
    "center": {"prev": [42.69, 32.84, 30.12], "inc": [456.79, 377.24, 354.06]},
    "south":  {"prev": [26.03, 23.15, 22.29], "inc": [318.27, 293.35, 284.87]},
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

print("Running simulation for 3 years (2010-2012)...")
print(f"  BitingRates: North={BITING_RATE_NORTH}, "
      f"Center={BITING_RATE_CENTER}, South={BITING_RATE_SOUTH}")
print(f"  prop_E: North={PROP_E_NORTH}, Center={PROP_E_CENTER}, South={PROP_E_SOUTH}")
print(f"  Reporting Rate: {REPORTING_RATE}")

results = oseirvector.simulate(
    t0                           = T0,
    tmax                         = TMAX,
    dt                           = DT,
    TimeExposed                  = TIME_EXPOSED,
    TimeInfectedAsymptomatic     = TIME_INFECTED_ASYMPTOMATIC,
    TimeInfectedSymptomatic      = TIME_INFECTED_SYMPTOMATIC,
    TimeWaningImmunity           = TIME_WANING_IMMUNITY,
    TransmissionVectorToHuman    = TRANSMISSION_VECTOR_TO_HUMAN,
    TransmissionHumanToVector    = TRANSMISSION_HUMAN_TO_VECTOR,
    BitingRateNorth              = BITING_RATE_NORTH,
    BitingRateCenter             = BITING_RATE_CENTER,
    BitingRateSouth              = BITING_RATE_SOUTH,
    prop_E_north                 = PROP_E_NORTH,
    prop_E_center                = PROP_E_CENTER,
    prop_E_south                 = PROP_E_SOUTH,
)



# Convert to numpy — shape becomes (timepoints, compartments) after transpose
region_data = {}
for i, name in enumerate(REGION_NAMES):
    region_data[name] = results[i].as_ndarray().T  # shape: (365, 28)

print(f"  Simulation done. Output shape: {region_data['north'].shape}")

# Quick day-0 check
d = region_data["north"]
N_g1 = d[0,8]+d[0,9]+d[0,10]+d[0,11]+d[0,12]
inf_g1 = d[0,10]+d[0,11]
print(f"North day-0 PfPR: {inf_g1/N_g1*100:.2f}%  (should be ~54%)")

d = region_data["south"]
N_g1 = d[0,8]+d[0,9]+d[0,10]+d[0,11]+d[0,12]
inf_g1 = d[0,10]+d[0,11]
print(f"South day-0 PfPR: {inf_g1/N_g1*100:.2f}%  (should be ~26%)")
# =============================================================================
# MONTHS AXIS (approximate: 30 days per month)
# =============================================================================

days   = np.arange(region_data["north"].shape[0])
months = days / 30.437  # convert days to months (average month length)

# =============================================================================
# PLOTTING — one figure per region, 3 subplots (one per age group)
# =============================================================================

for region in REGION_NAMES:
    data = region_data[region]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    fig.suptitle(f"Benin Malaria — {region.capitalize()} Region (2010-2012)", fontsize=14)
    for g, ax in enumerate(axes):
        idx = COMPARTMENTS[g]

        # Total population of this age group (for normalizing to proportions)
        N_g = (data[:, idx["S"]] + data[:, idx["E"]] +
               data[:, idx["IA"]] + data[:, idx["IS"]] + data[:, idx["R"]])

        # Plot each compartment as proportion of age group population
        ax.plot(months, data[:, idx["S"]]  / N_g * 100, label="S",  color="steelblue",  linewidth=1.5)
        ax.plot(months, data[:, idx["E"]]  / N_g * 100, label="E",  color="purple",     linewidth=1.5)
        ax.plot(months, data[:, idx["IA"]] / N_g * 100, label="IA", color="orange",     linewidth=1.5)
        ax.plot(months, data[:, idx["IS"]] / N_g * 100, label="IS", color="red",        linewidth=1.5)
        ax.plot(months, data[:, idx["R"]]  / N_g * 100, label="R",  color="green",      linewidth=1.5)

        ax.set_title(AGE_GROUP_NAMES[g])
        ax.set_xlabel("Month")
        ax.set_ylabel("% of age group")
        ax.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24])
        ax.set_xticklabels(["Jan\n2010", "Apr", "Jul", "Oct",
                    "Jan\n2011", "Apr", "Jul", "Oct",
                    "Jan\n2012"])
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fname = f"{region}_compartments.png" 
    plt.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close()

# =============================================================================
# DISTANCE PLOTS — model vs observed per year, per region
# =============================================================================

years_labels = [2010, 2011, 2012]

for region in REGION_NAMES:
    data = region_data[region]

    
    idx_g1                = COMPARTMENTS[1]

    model_prev = []
    model_inc  = []

    for year_idx in range(3):
        start  = year_idx * 365
        end    = start + 365
        y_data = data[start:end, :]

        # Prevalence
        N_g1   = (y_data[:, idx_g1["S"]] + y_data[:, idx_g1["E"]] +
                  y_data[:, idx_g1["IA"]] + y_data[:, idx_g1["IS"]] +
                  y_data[:, idx_g1["R"]])
        inf_g1 = y_data[:, idx_g1["IA"]] + y_data[:, idx_g1["IS"]]
        model_prev.append(np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100)

        # Incidence
        daily_clinical = sum(
            SYMPTOMATIC_FRACTIONS[g] * y_data[:, E_INDICES[g]] / TIME_EXPOSED
            for g in range(3)
        )
        pop_at_risk = y_data[0, HUMAN_INDICES].sum()
        model_inc.append((np.sum(daily_clinical) * REPORTING_RATE / pop_at_risk) * 1000)

    obs_prev = OBSERVED[region]["prev"]
    obs_inc  = OBSERVED[region]["inc"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Model vs Observed — {region.capitalize()} Region (2010-2012)", fontsize=14)

    for ax, model_vals, obs_vals, ylabel, title in zip(
        axes,
        [model_prev, model_inc],
        [obs_prev,   obs_inc],
        ["Prevalence PfPR 2-10 (%)", "Incidence (per 1000/year)"],
        ["Prevalence", "Incidence"]
    ):
        for yr_idx, year in enumerate(years_labels):
            mv = model_vals[yr_idx]
            ov = obs_vals[yr_idx]

            # Vertical line connecting the two points
            ax.plot([year, year], [mv, ov],
                    color="grey", linewidth=1.5, zorder=1)

            # Model point
            ax.scatter(year, mv,
                       color="steelblue", s=80, zorder=3,
                       label="Model" if yr_idx == 0 else "")

            # Observed point
            ax.scatter(year, ov,
                       color="C1", s=80, zorder=3,
                       label="Observed" if yr_idx == 0 else "")

        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xticks(years_labels)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fname = f"distance_{region}.png"
    plt.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close()
# =============================================================================
# SUMMARY TABLE — prevalence and incidence per region
# =============================================================================

print("\n" + "=" * 75)
print(f"{'SUMMARY TABLE — 2010 to 2012':^75}")
print("=" * 75)
print(f"{'Region':<12} {'Year':<8} {'Prevalence PfPR 2-10 (%)':>26} {'Incidence (/1000/yr)':>22}")
print("-" * 75)

for region in REGION_NAMES:
    data = region_data[region]

    idx_g1                = COMPARTMENTS[1]

    yearly_prev = []
    yearly_inc  = []

    for year_idx in range(3):
        start = year_idx * 365
        end   = start + 365
        y_data = data[start:end, :]

        # Prevalence PfPR 2-10
        N_g1 = (y_data[:, idx_g1["S"]] + y_data[:, idx_g1["E"]] +
                y_data[:, idx_g1["IA"]] + y_data[:, idx_g1["IS"]] +
                y_data[:, idx_g1["R"]])
        inf_g1     = y_data[:, idx_g1["IA"]] + y_data[:, idx_g1["IS"]]
        prevalence = np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100

        # Incidence
        daily_clinical = sum(
            SYMPTOMATIC_FRACTIONS[g] * y_data[:, E_INDICES[g]] / TIME_EXPOSED
            for g in range(3)
        )
        total_reported = np.sum(daily_clinical) * REPORTING_RATE
        pop_at_risk    = y_data[0, HUMAN_INDICES].sum()
        incidence      = (total_reported / pop_at_risk) * 1000

        yearly_prev.append(prevalence)
        yearly_inc.append(incidence)

        year_label = str(2010 + year_idx)
        print(f"  {region.capitalize():<10} {year_label:<8} "
              f"{prevalence:>26.2f} {incidence:>22.2f}")

    print("-" * 75)

print("=" * 75)
print("\nDone. Open the .png files in VS Code to view the plots.")
