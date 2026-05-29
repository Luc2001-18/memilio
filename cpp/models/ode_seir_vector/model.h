/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SEIR_VECTOR_MODEL_H
#define SEIR_VECTOR_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir_vector/infection_state.h"
#include "ode_seir_vector/parameters.h"

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")
#include <Eigen/Dense>
GCC_CLANG_DIAGNOSTIC(pop)

namespace mio
{
namespace oseirvector
{

/********************
 * define the model *
 ********************/

// clang-format off
using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Exposed>,
                       Flow<InfectionState::Exposed,              InfectionState::InfectedAsymptomatic>,
                       Flow<InfectionState::Exposed,              InfectionState::InfectedSymptomatic>,
                       Flow<InfectionState::InfectedAsymptomatic, InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSymptomatic,  InfectionState::Recovered>,
                       Flow<InfectionState::Recovered,            InfectionState::Susceptible>,
                       Flow<InfectionState::Susceptible_vector,   InfectionState::Infected_vector>>;
// **TODO**: Add/Adjust the flows as needed for the model.

// clang-format on
template <typename FP>
class Model
    : public FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    // **TODO**: Adjust the get_flows function according to the model.
    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> /*pop*/, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t ,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const Index<AgeGroup> age_groups = reduce_index<Index<AgeGroup>>(this->populations.size());
        const auto& params               = this->parameters;
        //  Fetch Malaria-specific parameters
        const int year_idx = std::min(static_cast<int>(std::floor(t / 365.0)), 18);
        const FP a = params.template get<EffectiveBitingRate<FP>>()[year_idx];
       // const FP a    = params.template get<MosquitoBitingRate<FP>>();


       // const FP hbr  = params.template get<HumanBitingRate<FP>>();
        const FP p_vh = params.template get<TransmissionVectorToHuman<FP>>();
        const FP p_hv = params.template get<TransmissionHumanToVector<FP>>();
        // Fetch Seasonality Parameters 
        const FP amp1 = params.template get<SeasonalityAmp1<FP>>();
        const FP amp2 = params.template get<SeasonalityAmp2<FP>>();
        const FP phi1 = params.template get<SeasonalityPhi1<FP>>();
        const FP phi2 = params.template get<SeasonalityPhi2<FP>>();
        const FP peak = params.template get<SeasonalityPeak<FP>>();
        // CALCULATE SEASONALITY FOR THIS TIME STEP (t)
        const FP term1 = amp1 * std::pow(std::abs(std::cos((2.0 * M_PI * t / 365.0) - phi1)), peak);
        const FP term2 = amp2 * std::pow(std::abs(std::cos((4.0 * M_PI * t / 365.0) - phi2)), peak);
        const FP seas  = 1.0 + term1 + term2;
        // Fetch demographic parameters for mosquitoes
      //  const FP mu_b = params.template get<MosquitoBirthRate<FP>>();
       // const FP mu_d = params.template get<MosquitoDeathRate<FP>>();
      
        // Identify the vector group (the very last group)
        const size_t num_groups = (size_t)params.get_num_groups();
        const size_t vector_idx = num_groups - 1; 

        // Pre-calculate the total human population and total infected humans
        FP total_human_population = 0.0;
        FP total_infected_humans  = 0.0;

        for (size_t h = 0; h < vector_idx; ++h) {
            // Get the indices for the human compartments in age group 'h'
            const size_t S_h  = this->populations.get_flat_index({mio::AgeGroup(h), InfectionState::Susceptible});
            const size_t E_h  = this->populations.get_flat_index({mio::AgeGroup(h), InfectionState::Exposed});
            const size_t IA_h = this->populations.get_flat_index({mio::AgeGroup(h), InfectionState::InfectedAsymptomatic});
            const size_t IS_h = this->populations.get_flat_index({mio::AgeGroup(h), InfectionState::InfectedSymptomatic});
            const size_t R_h  = this->populations.get_flat_index({mio::AgeGroup(h), InfectionState::Recovered});
            
            // Sum up the values
            total_human_population += y[S_h] + y[E_h] + y[IA_h] + y[IS_h] + y[R_h];
            total_infected_humans  += y[IA_h] + y[IS_h];
        }
        
        // Prevent division by zero later
        const FP div_total_human = (total_human_population < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / total_human_population);

        for (auto i : make_index_range(age_groups)) {
           // Indices for Human Compartments
            const size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            const size_t IAi = this->populations.get_flat_index({i, InfectionState::InfectedAsymptomatic});
            const size_t ISi = this->populations.get_flat_index({i, InfectionState::InfectedSymptomatic});
            const size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});
            // Indices for Vector Compartments
            const size_t Sv_i = this->populations.get_flat_index({i, InfectionState::Susceptible_vector});
          //  const size_t Iv_i = this->populations.get_flat_index({i, InfectionState::Infected_vector});
            // Calculate Populations
         //   const FP Nh_i    = y[Si] + y[Ei] + y[IAi] + y[ISi] + y[Ri];
          //  const FP divNh_i = (Nh_i < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nh_i);
            
          //  const FP Nv_i    = y[Sv_i] + y[Iv_i]; // Total vectors
            
            // Calculate Forces of Infection

            FP FOI_H = 0.0;
            FP FOI_V = 0.0;

            if ((size_t)i != vector_idx) {
                // FOR HUMAN GROUPS (0, 1, 2):
                // Get the global infected mosquito population from the Vector Group
                const size_t Iv_global = this->populations.get_flat_index({mio::AgeGroup(vector_idx), InfectionState::Infected_vector});
                
                // 2. Humans get infected by the global mosquito pool
                FOI_H = (a * seas) * p_vh * (y[Iv_global] * div_total_human);
            } 
            else {
                // Mosquitoes get infected by the combined human pool we calculated in Step 1
                FOI_V = (a * seas) * p_hv * (total_infected_humans * div_total_human);
            }
           
            
           // const FP FOI_H = (a*seas) * p_vh * (y[Iv_i] * divNh_i);
           // const FP FOI_V = (a*seas) * p_hv * ((y[IAi] + y[ISi]) * divNh_i);
            // Fetch parameters
            const FP p_asymp = params.template get<AsymptomaticProbability<FP>>()[i];
            const FP t_E     = params.template get<TimeExposed<FP>>()[i];
            const FP gamma_A = params.template get<TimeInfectedAsymptomatic<FP>>()[i];
            const FP gamma_S = params.template get<TimeInfectedSymptomatic<FP>>()[i];
            const FP t_W = params.template get<TimeWaningImmunity<FP>>()[i];
            // Assign the Flows
            // --- HUMAN FLOWS ---
            // S_H -> E_H
            flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(i)] =
                FOI_H * y[Si];

            // E_H -> I_A
            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::InfectedAsymptomatic>(i)] =
                p_asymp * (1.0 / t_E) * y[Ei];

            // E_H -> I_S
            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::InfectedSymptomatic>(i)] =
               (1.0- p_asymp) * (1.0 / t_E) * y[Ei];

            // I_A -> R_H
            flows[Base::template get_flat_flow_index<InfectionState::InfectedAsymptomatic, InfectionState::Recovered>(i)] =
                (1.0 / gamma_A) * y[IAi];

            // I_S -> R_H
            flows[Base::template get_flat_flow_index<InfectionState::InfectedSymptomatic, InfectionState::Recovered>(i)] =
                (1.0 / gamma_S) * y[ISi];

            //R_H -> S_H
            flows[Base::template get_flat_flow_index<InfectionState::Recovered, InfectionState::Susceptible>(i)] =
                (1.0 / t_W) * y[Ri];

            // --- VECTOR FLOWS ---
            // S_V -> I_V
            flows[Base::template get_flat_flow_index<InfectionState::Susceptible_vector, InfectionState::Infected_vector>(i)] =
                FOI_V * y[Sv_i];
           // I_V -> S_V (mosquito turnover: infected mosquitoes die and are replaced by susceptible ones)
           // flows[Base::template get_flat_flow_index<InfectionState::Infected_vector, InfectionState::Susceptible_vector>(i)] =
            //    mu_d * y[Iv_i];
        }
    }
    
    #if 0
    /**
    *@brief Computes the reproduction number at a given index time of the Model output obtained by the Simulation.
    *@param t_idx The index time at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided index time.
    */
    IOResult<FP> get_reproduction_number(size_t t_idx, const mio::TimeSeries<FP>& y)
    {
        if (!(t_idx < static_cast<size_t>(y.get_num_time_points()))) {
            return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
        }

        auto const& params = this->parameters;

        const size_t num_groups                    = (size_t)params.get_num_groups();
        constexpr size_t num_infected_compartments = 2;
        const size_t total_infected_compartments   = num_infected_compartments * num_groups;

        ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<ScalarType>>();

        Eigen::MatrixX<FP> F = Eigen::MatrixX<FP>::Zero(total_infected_compartments, total_infected_compartments);
        Eigen::MatrixX<FP> V = Eigen::MatrixX<FP>::Zero(total_infected_compartments, total_infected_compartments);

        for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            for (auto j = AgeGroup(0); j < AgeGroup(num_groups); j++) {

                const ScalarType Nj    = this->populations.get_group_total(j);
                const ScalarType divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

                FP coeffStoE = contact_matrix.get_matrix_at(SimulationTime<FP>(y.get_time(t_idx)))(i.get(), j.get()) *
                               params.template get<TransmissionProbabilityOnContact<FP>>()[i] * divNj;
                F((size_t)i, (size_t)j + num_groups) = coeffStoE * y.get_value(t_idx)[Si];
            }

            FP T_Ei                                           = params.template get<TimeExposed<FP>>()[i];
            FP T_Ii                                           = params.template get<TimeInfected<FP>>()[i];
            V((size_t)i, (size_t)i)                           = 1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i)              = -1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i + num_groups) = 1.0 / T_Ii;
        }

        V = V.inverse();

        Eigen::MatrixXd NextGenMatrix(total_infected_compartments, total_infected_compartments);
        NextGenMatrix.noalias() = F * V;

        // Compute the largest eigenvalue in absolute value
        Eigen::ComplexEigenSolver<Eigen::MatrixX<ScalarType>> ces;
        ces.compute(NextGenMatrix);
        FP rho = ces.eigenvalues().cwiseAbs().maxCoeff();

        return mio::success(rho);
    }

    /**
    *@brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns vector containing all reproduction numbers
    */
    Eigen::VectorX<FP> get_reproduction_numbers(const mio::TimeSeries<FP>& y)
    {
        auto num_time_points = y.get_num_time_points();
        Eigen::VectorX<FP> temp(num_time_points);
        for (size_t i = 0; i < static_cast<size_t>(num_time_points); i++) {
            temp[i] = get_reproduction_number(i, y).value();
        }
        return temp;
    }

    /**
    *@brief Computes the reproduction number at a given time point of the Model output obtained by the Simulation. If the particular time point is not inside the output, a linearly interpolated value is returned.
    *@param t_value The time point at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided time point, potentially using linear interpolation.
    */
    IOResult<FP> get_reproduction_number(FP t_value, const mio::TimeSeries<FP>& y)
    {
        if (t_value < y.get_time(0) || t_value > y.get_last_time()) {
            return mio::failure(mio::StatusCode::OutOfRange,
                                "Cannot interpolate reproduction number outside computed horizon of the TimeSeries");
        }

        if (t_value == y.get_time(0)) {
            return mio::success(get_reproduction_number((size_t)0, y).value());
        }

        auto times = std::vector<FP>(y.get_times().begin(), y.get_times().end());

        auto time_late = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), t_value));

        FP y1 = get_reproduction_number(static_cast<size_t>(time_late - 1), y).value();
        FP y2 = get_reproduction_number(static_cast<size_t>(time_late), y).value();

        auto result = linear_interpolation(t_value, y.get_time(time_late - 1), y.get_time(time_late), y1, y2);
        return mio::success(static_cast<FP>(result));
    }
    
    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("Parameters", this->parameters);
        obj.add_element("Populations", this->populations);
    }
    
    #endif
    

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Model");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return Model{pop_, par_};
            },
            par, pop);
    }
};

} // namespace oseirvector
} // namespace mio

#endif // SEIR_VECTOR_MODEL_H
