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
from memilio.simulation import oseirvector


def run_benin_simulation():
    """
    Runs the Benin Malaria model simulation (Stochastic Metapopulation Model) via Python bindings.
    """
    print("Starting Simulation...")

    results = oseirvector.simulate(
        t0 = 0.0,
        tmax= 1200,
        dt=1.0,
        TimeExposed = 15.0,
        TimeInfectedAsymptomatic = 100.0,  
        TimeInfectedSymptomatic = 7.0,
        TransmissionProbabilityOnContact = 0.1,
        AsymptomaticProbability = 0.287,
        TimeWaningImmunity = 180.0,
        BitingRateNorth = 0.4,
        BitingRateCenter = 0.4,
        BitingRateSouth = 0.4
    )
    result_array = results.as_ndarray()
    results.print_table(column_labels = ["S_K", "E_K", "I_A_K", "I_S_K", "R_K", "Sv_K_Ghost", "Iv_K_Ghost",
                        "S_5_15", "E_5_15", "I_A_5_15", "I_S_5_15", "R_5_15", "Sv_5_Ghost", "Iv_5_Ghost",
                        "S_Ad", "E_Ad", "I_A_Ad", "I_S_Ad", "R_Ad", "Sv_Ad_Ghost", "Iv_Ad_Ghost",
                        "S_M_Ghost", "E_M_Ghost", "I_A_M_Ghost", "I_S_M_Ghost", "R_M_Ghost", "Sv_Mosq", "Iv_Mosq"])

if __name__ == "__main__":
    run_benin_simulation()
