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

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph.h"

auto simulate(ScalarType t0 = 0, ScalarType tmax = 6940, ScalarType dt = 0.1, ScalarType TimeExposed = 15.0,
              ScalarType TimeInfectedAsymptomatic = 100.0, ScalarType TimeInfectedSymptomatic = 7.0,
              ScalarType TransmissionProbabilityOnContact = 0.1, ScalarType AsymptomaticProbability = 0.287,
              ScalarType TimeWaningImmunity = 180.0)
{
    mio::set_log_level(mio::LogLevel::warn);

    // ScalarType t0   = 0; // Day 0 :  2010
    // ScalarType tmax = 6940.0; //6940.0; // End of 2018 //365;
    // ScalarType dt   = 0.1;

    mio::log_info("Simulating ODE SEIR-Vector (Malaria); t={} ... {} with dt = {}.", t0, tmax, dt);
    // Initialize 4 groups (0, 1, 2 = Humans [less than 2 years, from 2 to 10, adult]; 3 = Mosquitoes)
    mio::oseirvector::Model<ScalarType> model(4);
    // Human Populations
    // ScalarType total_population                                                        = 100000;
    std::array<ScalarType, 3> age_group_props = {0.2075, 0.2590, 0.5335}; // according to Benin's data

    std::array<ScalarType, 3> prop_E  = {0.020, 0.020, 0.010}; // 2%, 2%, 1% are currently incubating
    std::array<ScalarType, 3> prop_IA = {0.000, 0.287, 0.150}; // 0% (Per your assumption!), 28.7% (Ouidah study), 15%
    std::array<ScalarType, 3> prop_IS = {0.100, 0.020, 0.005}; // 10%, 2%, 0.5% currently sick with fever
    std::array<ScalarType, 3> prop_R  = {0.150, 0.300, 0.450}; // 15%, 30%, 45% protected by recent infection
    ScalarType prop_vector_infected   = 0.01;

    for (size_t i = 0; i < 3; ++i) {
        ScalarType prop_S = 1.0 - prop_E[i] - prop_IA[i] - prop_IS[i] - prop_R[i];
        assert(prop_S > 0.0 && "Proportions sum to >= 1.0 for this age group!");

        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Susceptible}]          = prop_S;
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Exposed}]              = prop_E[i];
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedAsymptomatic}] = prop_IA[i];
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedSymptomatic}]  = prop_IS[i];
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Recovered}]            = prop_R[i];
        // ZERO-LOCK the vector compartments
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Susceptible_vector}] = 0.0;
        model.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Infected_vector}]    = 0.0;
    }
    // Vector Population
    // ScalarType total_vector_population                                                          = 50000;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Infected_vector}] = prop_vector_infected;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Susceptible_vector}] =
        1.0 - prop_vector_infected;
    // ZERO - LOCK the human compartment
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Susceptible}]          = 0.0;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Exposed}]              = 0.0;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::InfectedAsymptomatic}] = 0.0;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::InfectedSymptomatic}]  = 0.0;
    model.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Recovered}]            = 0.0;
    // Parameters
    // Global parametera
    model.parameters.set<mio::oseirvector::MosquitoBitingRate<ScalarType>>(0.4); // Mosquito bites a human every ~4 days
    model.parameters.set<mio::oseirvector::TransmissionVectorToHuman<ScalarType>>(0.27); // Data estimated
    model.parameters.set<mio::oseirvector::TransmissionHumanToVector<ScalarType>>(0.64); // Data estimated
    model.parameters.set<mio::oseirvector::MosquitoBirthRate<ScalarType>>(
        0.05); // Even if unused in flows, good to set to avoid uninitialized memory
    model.parameters.set<mio::oseirvector::MosquitoDeathRate<ScalarType>>(0.05);
    // Set the seasonality values
    model.parameters.set<mio::oseirvector::SeasonalityAmp1<ScalarType>>(0.1);
    model.parameters.set<mio::oseirvector::SeasonalityAmp2<ScalarType>>(0.1);
    model.parameters.set<mio::oseirvector::SeasonalityPhi1<ScalarType>>(2.5);
    model.parameters.set<mio::oseirvector::SeasonalityPhi2<ScalarType>>(6.5);
    model.parameters.set<mio::oseirvector::SeasonalityPeak<ScalarType>>(3);
    // Age group dependent parameters
    for (size_t i = 0; i < 3; ++i) {
        // Set values for human groups 0, 1, and 2
        model.parameters.get<mio::oseirvector::TimeExposed<ScalarType>>()[mio::AgeGroup(i)] = TimeExposed; // Data based
        model.parameters.get<mio::oseirvector::TimeInfectedAsymptomatic<ScalarType>>()[mio::AgeGroup(i)] =
            TimeInfectedAsymptomatic; // Data based
        model.parameters.get<mio::oseirvector::TimeInfectedSymptomatic<ScalarType>>()[mio::AgeGroup(i)] =
            TimeInfectedSymptomatic; // Data based
        model.parameters.get<mio::oseirvector::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(i)] =
            TransmissionProbabilityOnContact;
        model.parameters.get<mio::oseirvector::AsymptomaticProbability<ScalarType>>()[mio::AgeGroup(i)] =
            AsymptomaticProbability; // e.g., 30% are asymptomatic
        model.parameters.get<mio::oseirvector::TimeWaningImmunity<ScalarType>>()[mio::AgeGroup(i)] =
            TimeWaningImmunity; // Data Estimated
    }
    // Age_group dependant parameters
    model.parameters.get<mio::oseirvector::AsymptomaticProbability<ScalarType>>()[mio::AgeGroup(0)] =
        0.0; // asumption of the model for under 2
    model.parameters.get<mio::oseirvector::AsymptomaticProbability<ScalarType>>()[mio::AgeGroup(1)] = 0.2; // Data based
    model.parameters.get<mio::oseirvector::AsymptomaticProbability<ScalarType>>()[mio::AgeGroup(2)] =
        0.35; // Data based

    // All this is Dummy values
    model.parameters.get<mio::oseirvector::TimeExposed<ScalarType>>()[mio::AgeGroup(3)] = 1.0; // Dummy value
    model.parameters.get<mio::oseirvector::TimeInfectedAsymptomatic<ScalarType>>()[mio::AgeGroup(3)]         = 1.0;
    model.parameters.get<mio::oseirvector::TimeInfectedSymptomatic<ScalarType>>()[mio::AgeGroup(3)]          = 1.0;
    model.parameters.get<mio::oseirvector::AsymptomaticProbability<ScalarType>>()[mio::AgeGroup(3)]          = 0.0;
    model.parameters.get<mio::oseirvector::TimeWaningImmunity<ScalarType>>()[mio::AgeGroup(3)]               = 1.0;
    model.parameters.get<mio::oseirvector::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(3)] = 0.0;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::oseirvector::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(0.0);
    //  contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(30.));

    // model.check_constraints();

    // Metapopulation Set up !
    auto model_north  = model;
    auto model_center = model;
    auto model_south  = model;

    //  PATCH-SPECIFIC OVERRIDES

    auto set_patch_population = [&](mio::oseirvector::Model<ScalarType>& patch, ScalarType total_humans,
                                    ScalarType total_vectors) {
        for (size_t i = 0; i < 3; ++i) {
            ScalarType group_size = total_humans * age_group_props[i];
            patch.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Susceptible}] *= group_size;
            patch.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Exposed}] *= group_size;
            patch.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedAsymptomatic}] *= group_size;
            patch.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedSymptomatic}] *= group_size;
            patch.populations[{mio::AgeGroup(i), mio::oseirvector::InfectionState::Recovered}] *= group_size;
        }
        patch.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Susceptible_vector}] *= total_vectors;
        patch.populations[{mio::AgeGroup(3), mio::oseirvector::InfectionState::Infected_vector}] *= total_vectors;
    };

    // NORTH BENIN
    set_patch_population(model_north, 4639811.0, 9279622.0); // assume 2 times mosquitoes than human

    // CENTER BENIN
    set_patch_population(model_center, 2457289.0, 4914578.0); // assume 2 times mosquitoes than human

    // SOUTH BENIN
    set_patch_population(model_south, 7548534.0, 22645602.0); // assume 3 times mosquitoes than human, more humid

    model_north.check_constraints();
    model_center.check_constraints();
    model_south.check_constraints();

    // SETTING UP THE METAPOPULATION GRAPH

    // track Asymptomatic (I_A) and Symptomatic (I_S) humans who travel.
    std::vector<std::vector<size_t>> indices_save_edges(2);
    for (auto& vec : indices_save_edges) {
        vec.reserve(3); // Reserve for our 3 human age groups
    }

    // Only loop through age groups 0, 1, and 2 (Humans).
    for (size_t i = 0; i < 3; ++i) {
        indices_save_edges[0].emplace_back(model.populations.get_flat_index(
            {mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedAsymptomatic}));
        indices_save_edges[1].emplace_back(model.populations.get_flat_index(
            {mio::AgeGroup(i), mio::oseirvector::InfectionState::InfectedSymptomatic}));
    }

    // Build the graph
    mio::Graph<mio::SimulationNode<ScalarType, mio::Simulation<ScalarType, mio::oseirvector::Model<ScalarType>>>,
               mio::MobilityEdge<ScalarType>>
        g;
    // Add the 3 patches as nodes. (ID, model, start_time)
    g.add_node(1001, model_north, t0); // Node 0 in the graph
    g.add_node(1002, model_center, t0); // Node 1 in the graph
    g.add_node(1003, model_south, t0); // Node 2 in the graph

    // Create a helper function to set human mobility but Stop mosquitoes
    auto make_mobility_vector = [&model](ScalarType human_rate) {
        // Get the total size safely from the model
        Eigen::VectorX<ScalarType> rates =
            Eigen::VectorX<ScalarType>::Constant(model.populations.get_num_compartments(), human_rate);

        // Safely zero-lock the vector compartments specifically for Age Group 3 (Mosquitoes)
        rates[model.populations.get_flat_index(
            {mio::AgeGroup(3), mio::oseirvector::InfectionState::Susceptible_vector})] = 0.0;
        rates[model.populations.get_flat_index({mio::AgeGroup(3), mio::oseirvector::InfectionState::Infected_vector})] =
            0.0;

        return rates;
    };

    // Add the specific routes : go for Symmetric assumption for now..may change later

    // North (Node 0) <--> Center (Node 1): Rate 0.3
    g.add_edge(0, 1, make_mobility_vector(0.03), indices_save_edges);
    g.add_edge(1, 0, make_mobility_vector(0.03), indices_save_edges);

    // Center (Node 1) <--> South (Node 2): Rate 0.5
    g.add_edge(1, 2, make_mobility_vector(0.05), indices_save_edges);
    g.add_edge(2, 1, make_mobility_vector(0.05), indices_save_edges);

    // North (Node 0) <--> South (Node 2): Rate 0.1
    g.add_edge(0, 2, make_mobility_vector(0.01), indices_save_edges);
    g.add_edge(2, 0, make_mobility_vector(0.01), indices_save_edges);

    // Run the simulation
    ScalarType dt_2 = 5.0;
    auto sim        = mio::make_mobility_sim<ScalarType>(t0, dt_2, std::move(g));

    sim.advance(tmax);

    // OUTPUT RESULTS (Example: North Benin / Node 0)

    // Extract the raw results specifically from Node 0
    auto& results_north = sim.get_graph().nodes()[0].property.get_result();

    // Interpolate the results
    return mio::interpolate_simulation_result<ScalarType>(results_north);
}

#ifndef OSEIRVECTOR_BINDINGS_SKIP_MAIN

int main()
{
    auto interpolated_results_north = simulate();

    std::cout << "\n=== NORTH BENIN RESULTS (Node 0) ===\n";

    // Print the table using specific Malaria compartments
    interpolated_results_north.print_table(
        {// Group 0: Kids < 5
         "S_K", "E_K", "I_A_K", "I_S_K", "R_K", "Sv_K_Ghost", "Iv_K_Ghost",

         // Group 1: Ages 5 to 15
         "S_5_15", "E_5_15", "I_A_5_15", "I_S_5_15", "R_5_15", "Sv_5_Ghost", "Iv_5_Ghost",

         // Group 2: Adults 15+
         "S_Ad", "E_Ad", "I_A_Ad", "I_S_Ad", "R_Ad", "Sv_Ad_Ghost", "Iv_Ad_Ghost",

         // Group 3: Mosquitoes
         "S_M_Ghost", "E_M_Ghost", "I_A_M_Ghost", "I_S_M_Ghost", "R_M_Ghost", "Sv_Mosq", "Iv_Mosq"});

    std::cout << "\nNorth Benin Final Population: " << interpolated_results_north.get_last_value().sum() << "\n";

    /*
    // Commuter results per edge
    std::cout << "\n=== COMMUTERS: North --> Center ===\n";
    sim.get_graph().edges()[0].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});

    std::cout << "\n=== COMMUTERS: Center --> North ===\n";
    sim.get_graph().edges()[1].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});

    std::cout << "\n=== COMMUTERS: Center --> South ===\n";
    sim.get_graph().edges()[2].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});

    std::cout << "\n=== COMMUTERS: South --> Center ===\n";
    sim.get_graph().edges()[3].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});

    std::cout << "\n=== COMMUTERS: North --> South ===\n";
    sim.get_graph().edges()[4].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});

    std::cout << "\n=== COMMUTERS: South --> North ===\n";
    sim.get_graph().edges()[5].property.get_mobility_results()
        .print_table({"Total_I_A", "Total_I_S"});
  //  seir.print_table({"S_H", "E_H", "I_H", "R_H", "S_V", "I_V"});
   // std::cout << "\nnumber total: " << seir.get_last_value().sum() << "\n"; 
   */
    return interpolated_results_north.get_last_value().sum();
}
#endif