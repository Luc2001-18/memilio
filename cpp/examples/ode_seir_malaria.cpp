/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "ode_seir_vector/model.h"
#include "memilio/data/analyze_result.h"
#include "ode_seir_vector/infection_state.h"
#include "ode_seir_vector/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

#include "memilio/utils/time_series.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 100.;
    ScalarType dt   = 1.0;

    mio::log_info("Simulating ODE SEIR-Vector (Malaria); t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::oseirvector::Model<ScalarType> model(1);
    // Human Populations
    ScalarType total_population                                                        = 100000;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Exposed}]   = 1000;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Infected}]  = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Recovered}] = 1000;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Recovered}];
    // Vector Population
    ScalarType total_vector_population                                                          = 50000;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Infected_vector}]    = 500;
    model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Susceptible_vector}] = 
        total_vector_population - 
        model.populations[{mio::AgeGroup(0), mio::oseirvector::InfectionState::Infected_vector}];

    // Parameters    
    model.parameters.set<mio::oseirvector::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirvector::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirvector::MosquitoBitingRate<ScalarType>>(0.3);        // Mosquito bites a human every ~3 days
    model.parameters.set<mio::oseirvector::TransmissionVectorToHuman<ScalarType>>(0.1); // 10% chance of transmission to human
    model.parameters.set<mio::oseirvector::TransmissionHumanToVector<ScalarType>>(0.5); // 50% chance of transmission to mosquito
    model.parameters.set<mio::oseirvector::MosquitoBirthRate<ScalarType>>(0.05);        // Even if unused in flows, good to set to avoid uninitialized memory
    model.parameters.set<mio::oseirvector::MosquitoDeathRate<ScalarType>>(0.05);
    
    
    model.parameters.set<mio::oseirvector::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::oseirvector::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(30.));

    model.check_constraints();

    auto seir = mio::simulate<ScalarType>(t0, tmax, dt, model);

    // Interpolate results to get clean daily values (Day 0, 1, 2...)
    auto interpolated_results = mio::interpolate_simulation_result<ScalarType>(seir);

    // Print the table using your specific Malaria compartments
    interpolated_results.print_table({"S_H", "E_H", "I_H", "R_H", "S_V", "I_V"});
    
    // Print the total population to verify it is still 60000
    std::cout << "\nPopulation total: " << interpolated_results.get_last_value().sum() << "\n";

  //  seir.print_table({"S_H", "E_H", "I_H", "R_H", "S_V", "I_V"});
   // std::cout << "\nnumber total: " << seir.get_last_value().sum() << "\n";
}
