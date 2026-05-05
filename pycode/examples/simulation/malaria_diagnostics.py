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

BITING_RATE_NORTH  = 0.4
BITING_RATE_CENTER = 0.4
BITING_RATE_SOUTH  = 0.4
REPORTING_RATE = 0.0433
T0   = 0.0
TMAX = 730.0   # one year
DT   = 1.0     # daily resolution

# Fixed parameters (not being calibrated here)
TIME_EXPOSED                    = 15.0
TIME_INFECTED_ASYMPTOMATIC      = 80.0
TIME_INFECTED_SYMPTOMATIC       = 7.0
TRANSMISSION_PROB_ON_CONTACT    = 0.1
ASYMPTOMATIC_PROBABILITY        = 0.287
TIME_WANING_IMMUNITY            = 180 #730.0
TRANSMISSION_VECTOR_TO_HUMAN    = 0.27 #0.27
TRANSMISSION_HUMAN_TO_VECTOR    = 0.02 #0.02

# =============================================================================
# COMPARTMENT INDEX MAP
# Group 0 (kids <2):   S=0,  E=1,  IA=2,  IS=3,  R=4,  Sv=5,  Iv=6
# Group 1 (2-10):      S=7,  E=8,  IA=9,  IS=10, R=11, Sv=12, Iv=13
# Group 2 (adults):    S=14, E=15, IA=16, IS=17, R=18, Sv=19, Iv=20
# Group 3 (mosquitos): S=21, E=22, IA=23, IS=24, R=25, Sv=26, Iv=27
# =============================================================================

COMPARTMENTS = {
    0: {"S": 0,  "E": 1,  "IA": 2,  "IS": 3,  "R": 4},
    1: {"S": 7,  "E": 8,  "IA": 9,  "IS": 10, "R": 11},
    2: {"S": 14, "E": 15, "IA": 16, "IS": 17, "R": 18},
}

EXPOSED_INDICES = [1, 8, 15]
HUMAN_INDICES   = [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18]

AGE_GROUP_NAMES = ["Kids < 2y", "Children 2-10y", "Adults 10+"]
REGION_NAMES    = ["north", "center", "south"]

# =============================================================================
# RUN SIMULATION
# =============================================================================

print("Running simulation for 1 year...")
print(f"  BitingRates: North={BITING_RATE_NORTH}, "
      f"Center={BITING_RATE_CENTER}, South={BITING_RATE_SOUTH}")

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
    BitingRateNorth              = BITING_RATE_NORTH,
    BitingRateCenter             = BITING_RATE_CENTER,
    BitingRateSouth              = BITING_RATE_SOUTH
)

# Convert to numpy — shape becomes (timepoints, compartments) after transpose
region_data = {}
for i, name in enumerate(REGION_NAMES):
    region_data[name] = results[i].as_ndarray().T  # shape: (365, 28)

print(f"  Simulation done. Output shape: {region_data['north'].shape}")

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
    fig.suptitle(f"Benin Malaria — {region.capitalize()} Region (2010-2011)", fontsize=14)
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
    fname = f"{region}_compartmentss.png" 
    plt.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close()

# =============================================================================
# SUMMARY TABLE — prevalence and incidence per region
# =============================================================================

print("\n" + "=" * 75)
print(f"{'SUMMARY TABLE — 2010 to 2011':^75}")
print("=" * 75)
print(f"{'Region':<12} {'Year':<8} {'Prevalence PfPR 2-10 (%)':>26} {'Incidence (/1000/yr)':>22}")
print("-" * 75)

for region in REGION_NAMES:
    data = region_data[region]

    SYMPTOMATIC_FRACTIONS = [1.0, 0.8, 0.65]
    E_INDICES             = [1, 8, 15]
    idx_g1                = COMPARTMENTS[1]

    yearly_prev = []
    yearly_inc  = []

    for year_idx in range(2):
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

    # Average row
    avg_prev = np.mean(yearly_prev)
    avg_inc  = np.mean(yearly_inc)
    print(f"  {region.capitalize():<10} {'avg':<8} "
          f"{avg_prev:>26.2f} {avg_inc:>22.2f}")
    print("-" * 75)

print("=" * 75)
print("\nDone. Open the .png files in VS Code to view the plots.")