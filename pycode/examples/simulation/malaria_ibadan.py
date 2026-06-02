#############################################################################
# Malaria Ibadan Validation
# Uses calibrated parameters from 2015 to run model for 2015-2017
# Compares monthly and annual prevalence against observed Ibadan data
# Calibrated: biting_rate=0.2820, prop_E=0.1003
#############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from memilio.simulation import oseirvector

# =============================================================================
# 1. LOAD IBADAN DATA (2015-2017)
# =============================================================================

df = pd.read_excel(
    "ibadan_malaria-prevalence_dataset_1996-2017.xlsx",
    sheet_name="vars-preP_val"
)

ibadan = df[(df["year"] >= 2015) & (df["year"] <= 2017)].copy().reset_index(drop=True)
ibadan["prevalence_pct"] = ibadan["preP"] * 100
ibadan["month_index"]    = (ibadan["year"] - 2015) * 12 + ibadan["month"]

print("=== IBADAN DATA LOADED (2015-2017) ===")
print(f"  Rows: {len(ibadan)}")

# =============================================================================
# 2. CONFIGURATION — calibrated parameters from 2015
# =============================================================================

BITING_RATE  =  0.1820 #0.2820
PROP_E       = 0.1114

T0   = 0.0
TMAX = 1095.0   # 3 years = 2015-2017
DT   = 1.0

TIME_EXPOSED                 = 15.0
TIME_INFECTED_ASYMPTOMATIC   = 100.0
TIME_INFECTED_SYMPTOMATIC    = 7.0
TRANSMISSION_PROB_ON_CONTACT = 0.1
ASYMPTOMATIC_PROBABILITY     = 0.287
TIME_WANING_IMMUNITY         = 180.0
TRANSMISSION_VECTOR_TO_HUMAN = 0.27
TRANSMISSION_HUMAN_TO_VECTOR = 0.02

# Compartment indices — time is column 0
S_g1  = 8
E_g1  = 9
IA_g1 = 10
IS_g1 = 11
R_g1  = 12

# =============================================================================
# 3. RUN SIMULATION
# =============================================================================

print(f"\n=== RUNNING MODEL (2015-2017) ===")
print(f"  Biting rate: {BITING_RATE}")
print(f"  prop_E:      {PROP_E}")

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
    BitingRateNorth              = BITING_RATE,
    BitingRateCenter             = BITING_RATE,
    BitingRateSouth              = BITING_RATE,
    prop_E_north                 = PROP_E,
    prop_E_center                = PROP_E,
    prop_E_south                 = PROP_E,
)

# Use South region as Ibadan proxy
data = results[2].as_ndarray().T  # shape: (1096, 29)
print(f"  Simulation done. Shape: {data.shape}")

# =============================================================================
# 4. COMPUTE MONTHLY PREVALENCE (36 months)
# =============================================================================

model_monthly = []
for month_idx in range(36):
    day_start = month_idx * 30
    day_end   = day_start + 30

    y_data    = data[day_start:day_end, :]

    N_g1   = (y_data[:, S_g1] + y_data[:, E_g1] +
              y_data[:, IA_g1] + y_data[:, IS_g1] +
              y_data[:, R_g1])
    inf_g1 = y_data[:, IA_g1] + y_data[:, IS_g1]
    prev   = np.mean(inf_g1 / np.where(N_g1 == 0, 1, N_g1)) * 100
    model_monthly.append(prev)

model_monthly = np.array(model_monthly)

# =============================================================================
# 5. COMPUTE ANNUAL AVERAGES
# =============================================================================

years = [2015, 2016, 2017]
model_annual   = []
observed_annual = []

for yr_idx, yr in enumerate(years):
    start = yr_idx * 12
    end   = start + 12
    model_annual.append(model_monthly[start:end].mean())
    observed_annual.append(
        ibadan[ibadan["year"] == yr]["prevalence_pct"].mean()
    )

model_annual    = np.array(model_annual)
observed_annual = np.array(observed_annual)

# =============================================================================
# 6. PRINT SUMMARY TABLE
# =============================================================================

print("\n" + "=" * 65)
print(f"{'MONTHLY PREVALENCE — MODEL VS OBSERVED':^65}")
print("=" * 65)
print(f"  {'Month':<12} {'Model (%)':>10} {'Observed (%)':>14} {'Diff':>8}")
print("  " + "-" * 48)

month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

for yr_idx, yr in enumerate(years):
    print(f"\n  --- {yr} ---")
    for mo in range(12):
        global_idx = yr_idx * 12 + mo
        mv = model_monthly[global_idx]
        ov = ibadan.iloc[global_idx]["prevalence_pct"]
        diff = mv - ov
        print(f"  {month_names[mo]:<12} {mv:>10.2f} {ov:>14.2f} {diff:>+8.2f}")

print("\n" + "=" * 65)
print(f"{'ANNUAL AVERAGES':^65}")
print("=" * 65)
print(f"  {'Year':<8} {'Model (%)':>10} {'Observed (%)':>14} {'Diff':>8}")
print("  " + "-" * 40)
for i, yr in enumerate(years):
    diff = model_annual[i] - observed_annual[i]
    print(f"  {yr:<8} {model_annual[i]:>10.2f} {observed_annual[i]:>14.2f} {diff:>+8.2f}")
print("=" * 65)

# =============================================================================
# 7. PLOT 1 — MONTHLY PREVALENCE (2015-2017)
# =============================================================================

month_labels = []
for yr in years:
    for mo in month_names:
        month_labels.append(f"{mo}\n{yr}")

x = np.arange(1, 37)

fig, ax = plt.subplots(figsize=(16, 5))
ax.plot(x, model_monthly, color="steelblue", linewidth=2,
        marker="D", markersize=5, label="Model (calibrated on 2015)")
ax.scatter(x, ibadan["prevalence_pct"].values,
           color="C1", s=80, zorder=3, label="Observed Ibadan")
ax.set_xlabel("Month")
ax.set_ylabel("Prevalence PfPR (%)")
ax.set_title("Model vs Observed Monthly Prevalence — Ibadan 2015-2017")
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=7)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("ibadan_monthly_2015_2017.png", dpi=150)
plt.close()
print("\n  Saved: ibadan_monthly_2015_2017.png")

# =============================================================================
# 8. PLOT 2 — ANNUAL AVERAGES
# =============================================================================

fig, ax = plt.subplots(figsize=(8, 5))
x_yr = np.arange(len(years))
width = 0.35

bars1 = ax.bar(x_yr - width/2, model_annual, width,
               color="steelblue", alpha=0.8, label="Model")
bars2 = ax.bar(x_yr + width/2, observed_annual, width,
               color="C1", alpha=0.8, label="Observed")

# Add value labels on bars
for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9)
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9)

ax.set_xlabel("Year")
ax.set_ylabel("Mean Annual Prevalence (%)")
ax.set_title("Annual Average Prevalence — Model vs Observed Ibadan")
ax.set_xticks(x_yr)
ax.set_xticklabels([str(y) for y in years])
ax.legend()
ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig("ibadan_annual_2015_2017.png", dpi=150)
plt.close()
print("  Saved: ibadan_annual_2015_2017.png")

print("\nDone.")
